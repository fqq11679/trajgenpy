import math

import geojson
import pyproj
import shapely
import shapely.plotting as shplt
from shapely.geometry.polygon import orient
import random
import trajgenpy.bindings as bindings
from trajgenpy import Logging

import itertools
import networkx as nx

from ortools.constraint_solver import pywrapcp, routing_enums_pb2

log = Logging.get_logger()


class GeoData:
    def __init__(self, geometry, crs="WGS84"):
        self.geometry = geometry
        self.crs = crs

    def set_crs(self, crs):
        if not isinstance(crs, str):
            msg = "New CRS must be a string."
            raise ValueError(msg)

        if crs != self.crs:
            # Apply the transformer to the geometry
            self._convert_to_crs(crs)
        self.crs = crs
        return self

    def _convert_to_crs(self, crs):  # noqa: ARG002
        msg = "_convert_to_crs(crs) sould be implemented in the data classes!"
        raise NotImplementedError(msg)

    def is_geometry_of_type(self, geometry, expected_class):
        if expected_class and not isinstance(geometry, expected_class):
            msg = f"Geometry must be a {expected_class.__name__}."
            raise ValueError(msg)

    def get_geometry(self):
        return self.geometry

    def buffer(self, distance, quad_segs=1, cap_style="square", join_style="bevel"):
        self.geometry = self.geometry.buffer(distance, quad_segs, cap_style, join_style)
        return self

    def __str__(self):
        return f"Geometry in CRS: {self.crs}\nGeometry: {self.geometry}"

    def __geo_interface__(self):
        return self.geometry.__geo_interface__

    def to_geojson(self, id=None, name=None, properties=None):
        if id is None:
            id = random.randint(0, 1000000000)
        if name is None:
            name = str(id)
        if properties is None:
            properties = {}

        properties["crs"] = self.crs
        properties["name"] = name
        return geojson.Feature(id, self.geometry, properties=properties)


class GeoTrajectory(GeoData):
    def __init__(self, geometry, crs="WGS84"):
        self.is_geometry_of_type(geometry, shapely.LineString)
        super().__init__(geometry, crs)

    def plot(self, ax=None, add_points=True, color=None, linewidth=2, **kwargs):
        if self.crs == "WGS84":
            log.warning(
                "Plotting in WGS84 is not recomended as this distorts the geometry!"
            )
        shplt.plot_line(self.geometry, ax, add_points, color, linewidth, **kwargs)

    def _convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)

        converted_coords = [
            transformer.transform(x, y) for x, y in list(self.geometry.coords)
        ]
        self.geometry = shapely.LineString(converted_coords)

class DFSearcher:
    def __init__(self, trajs, polygon, sub_polys):
        self.trajs = trajs
        self.polygon = polygon
        self.sub_polys = sub_polys
        self.visited = set()
        self.final_path = []

    def connect(self, polygon, point1, point2):
        vertices = list(polygon.exterior.coords)[0:-1]
        all_points = vertices + [point1, point2]

        G = nx.Graph()

        for i in range(len(all_points)):
            for j in range(i+1, len(all_points)):
                line = shapely.geometry.LineString([all_points[i], all_points[j]])
                if polygon.buffer(0.001).contains(line):
                    G.add_edge(i, j, weight=line.length)

        point1_index = all_points.index(point1)
        point2_index = all_points.index(point2)

        path = nx.dijkstra_path(G, source=point1_index, target=point2_index, weight='weight')
        path_length = nx.dijkstra_path_length(G, source=point1_index, target=point2_index, weight='weight')

        inter_points = [all_points[i] for i in path]

        return inter_points, path_length

    def find_nearest_traj(self, start_point):
        nearest_traj = None
        nearest_distance = float('inf')
        is_rev = None

        for traj in self.trajs:
            distance = start_point.distance(shapely.Point(list(traj.coords)[0]))
            if distance < nearest_distance:
                nearest_distance = distance
                nearest_traj = traj
                is_rev = 0

            distance = start_point.distance(shapely.Point(list(traj.coords)[-1]))
            if distance < nearest_distance:
                nearest_distance = distance
                nearest_traj = traj
                is_rev = 1

        return nearest_traj, is_rev

    def dfs(self, current_traj, is_rev, start_point, end_point):
        print(self.final_path)

        # 标记当前几何形状为已访问
        self.visited.add(current_traj)

        # 获取当前几何的起终点
        is_start = 0

        # 初始化 MultiLineString 用于回溯线
        backtrack_lines = []
        forward_lines = []

        if is_rev == 0:
            head_point = list(current_traj.coords)[0]
            tail_point = list(current_traj.coords)[-1]
        else:
            head_point = list(current_traj.coords)[-1]
            tail_point = list(current_traj.coords)[0]            

        # 添加连接路径
        if self.final_path:
            prev_point = self.final_path[-1]
            path, _ = self.connect(self.polygon, prev_point, head_point)
            self.final_path.extend(path)
            if len(path) > 1:
                forward_lines.append(shapely.LineString(path))
        else:
            is_start = 1
            self.final_path.append(start_point)  # 添加起始点

        if is_rev == 0:
            self.final_path.extend(list(current_traj.coords))
        else:
            self.final_path.extend(list(current_traj.coords)[::-1])

        # 检查是否已访问所有几何形状
        if len(self.visited) == len(self.trajs):
            # 连接到终点
            self.final_path.append(end_point)
            return True # 已完成遍历


        # 深度优先遍历相邻几何形状
        for next_traj in self.trajs:
            if next_traj not in self.visited:
                # 检查相邻性
                cur_inx = list(self.trajs).index(current_traj)
                next_inx = list(self.trajs).index(next_traj)
                if self.sub_polys[cur_inx].touches(self.sub_polys[next_inx]):  ## 可以根据需要使用其他相邻检查
                    _, next_length0 = self.connect(self.polygon, tail_point, list(next_traj.coords)[0])
                    _, next_length1 = self.connect(self.polygon, tail_point, list(next_traj.coords)[-1])
                    if next_length0 < next_length1:
                        if self.dfs(next_traj, 0, start_point, end_point):
                            return True
                    else:
                        if self.dfs(next_traj, 1, start_point, end_point):
                            return True

        print('------------------------')

        # 回溯：从当前几何形状回到起始点
        if is_start == 0:
            ret_path, _ = self.connect(self.polygon, tail_point, head_point)
            self.final_path.extend(ret_path)
            if len(ret_path) > 1:
                backtrack_lines.append(shapely.LineString(ret_path))
            ret_path, _ = self.connect(self.polygon, head_point, prev_point)
            self.final_path.extend(ret_path)
            if len(ret_path) > 1:
                backtrack_lines.append(shapely.LineString(ret_path))


        # 创建 MultiLineString 包含所有回溯线段
        if backtrack_lines:
            backtrack_multilines = shapely.MultiLineString(backtrack_lines)
            self.plot_multilines(backtrack_multilines)

        if forward_lines:
            forward_multilines = shapely.MultiLineString(forward_lines)
            self.forward_plot_multilines(forward_multilines)

        return False

    def plot_multilines(self, multilines):
        import matplotlib.pyplot as plt
        from shapely.affinity import rotate
        for line in multilines.geoms:
            line = rotate(line, angle=-45, origin=(0, 0))
            x, y = line.xy
            plt.plot(x, y, linestyle='dashed', color='red', linewidth=2)

        plt.show()

    def forward_plot_multilines(self, multilines):
        import matplotlib.pyplot as plt
        from shapely.affinity import rotate
        for line in multilines.geoms:
            line = rotate(line, angle=-45, origin=(0, 0))
            x, y = line.xy
            plt.plot(x, y, linestyle='dashed', color='blue', linewidth=2)

        plt.show()

    def run(self, start_point, end_point):
        nearest_traj, is_rev = self.find_nearest_traj(start_point)
        if nearest_traj:
            self.dfs(nearest_traj, is_rev, start_point, end_point)
        return self.final_path


class GeoMultiTrajectory(GeoData):
    def __init__(
        self,
        geometry: (
            shapely.MultiLineString
            | list[shapely.LineString]
            | list[GeoTrajectory]
            | shapely.LineString
        ),
        crs="WGS84",
    ):
        super().__init__(geometry, crs)
        if isinstance(geometry, list):
            for line in geometry:
                self.is_geometry_of_type(line, shapely.LineString)

            super().__init__(shapely.MultiLineString(geometry), crs)
        elif isinstance(geometry, shapely.LineString):
            self.is_geometry_of_type(geometry, shapely.LineString)
            super().__init__(shapely.MultiLineString([geometry]), crs)
        elif isinstance(geometry, GeoTrajectory):
            self.is_geometry_of_type(geometry, GeoTrajectory)
            super().__init__(shapely.MultiLineString([geometry.geometry]), crs)
        else:
            self.is_geometry_of_type(geometry, shapely.MultiLineString)
            super().__init__(geometry, crs)

    def concatenate_trajectories(self, polygon, sub_polys, start_point, end_point):
        searcher = DFSearcher(self.geometry.geoms, polygon, sub_polys)
        searcher.run(start_point, end_point)
        self.geometry = shapely.LineString(searcher.final_path)

    def plot(self, ax=None, add_points=False, color=None, linewidth=2, **kwargs):
        if self.crs == "WGS84":
            log.warning(
                "Plotting in WGS84 is not recomended as this distorts the geometry!"
            )
        if isinstance(self.geometry, shapely.MultiLineString):
            for line in self.geometry.geoms:
                shplt.plot_line(line, ax, add_points, color, linewidth, **kwargs)
        elif isinstance(self.geometry, shapely.geometry.LineString):
            shplt.plot_line(self.geometry, ax, add_points, color, linewidth, **kwargs)
        else:
            raise ValueError("Geometry is neither MultiLineString or LineString.")

    def _convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)
        # Convert the coordinates of each line in the MultiLineString
        converted_geoms = [
            [
                transformer.transform(x, y) for x, y in list(line.coords)
            ]  # Convert the coordinates of each line in the MultiLineString
            for line in self.geometry.geoms
        ]
        self.geometry = shapely.MultiLineString(converted_geoms)


class GeoPoint(GeoData):
    def __init__(self, geometry: shapely.Point, crs="WGS84"):
        self.is_geometry_of_type(geometry, shapely.Point)
        super().__init__(geometry, crs)

    def plot(self, ax=None, add_points=True, color=None, linewidth=2, **kwargs):
        if self.crs == "WGS84":
            log.warning(
                "Plotting in WGS84 is not recomended as this distorts the geometry!"
            )
        shplt.plot_points(self.geometry, ax, add_points, color, linewidth, **kwargs)

    def _convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)
        x, y = transformer.transform(self.geometry.x, self.geometry.y)
        self.geometry = shapely.Point(x, y)


class GeoPolygon(GeoData):
    def __init__(self, geometry: shapely.Polygon | shapely.LineString, crs="WGS84"):
        self.is_geometry_of_type(geometry, shapely.Polygon | shapely.LineString)
        geometry = shapely.Polygon(geometry)
        super().__init__(geometry, crs)

    def _convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)
        # Convert each point in the polygon
        exterior = [
            transformer.transform(x, y) for x, y in self.geometry.exterior.coords
        ]
        interiors = [
            [transformer.transform(x, y) for x, y in interior.coords]
            for interior in self.geometry.interiors
        ]
        self.geometry = shapely.Polygon(exterior, interiors)

    def plot(
        self,
        ax=None,
        add_points=False,
        color=None,
        facecolor=None,
        edgecolor=None,
        linewidth=2,
        **kwargs,
    ):
        if self.crs == "WGS84":
            log.warning(
                "Plotting in WGS84 is not recomended as this distorts the geometry!"
            )
        shplt.plot_polygon(
            polygon=self.geometry,
            ax=ax,
            add_points=add_points,
            color=color,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            **kwargs,
        )


class GeoMultiPolygon(GeoData):
    def __init__(self, geometry, crs="WGS84"):
        if isinstance(geometry, list):
            for geom in geometry:
                self.is_geometry_of_type(geom, shapely.Polygon)
            geometry = shapely.MultiPolygon(geometry)
        else:
            self.is_geometry_of_type(geometry, shapely.MultiPolygon)
        super().__init__(geometry, crs)

    def _convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)
        polygon_list = []
        for polygon in list(self.geometry.geoms):
            # Convert each point in the polygon
            exterior = [transformer.transform(x, y) for x, y in polygon.exterior.coords]
            interiors = [
                [transformer.transform(x, y) for x, y in interior.coords]
                for interior in polygon.interiors
            ]
            polygon_list.append(shapely.Polygon(exterior, interiors))

        self.geometry = shapely.MultiPolygon(polygon_list)

    def plot(
        self,
        ax=None,
        add_points=False,
        color=None,
        facecolor=None,
        edgecolor=None,
        linewidth=2,
        **kwargs,
    ):
        if self.crs == "WGS84":
            log.warning(
                "Plotting in WGS84 is not recomended as this distorts the geometry!"
            )
        shplt.plot_polygon(
            polygon=self.geometry,
            ax=ax,
            add_points=add_points,
            color=color,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            **kwargs,
        )


def multi_polygon_to_polygon_with_holes(multi_polygon):
    # Create an empty CGAL PolygonWithHoles
    polygon_with_holes = bindings.PolygonWithHoles(bindings.Polygon_2())

    # Iterate through each polygon in the MultiPolygon
    for polygon in multi_polygon:
        # Convert the Shapely polygon to a list of CGAL points
        cgal_points = [
            bindings.Point_2(point.x, point.y) for point in polygon.exterior.coords
        ]

        # Create a CGAL Polygon_2 from the list of points
        cgal_polygon = bindings.Polygon_2(cgal_points)

        # Add the polygon to the PolygonWithHoles as a hole
        polygon_with_holes.add_hole(cgal_polygon)

    # Return the resulting PolygonWithHoles
    return polygon_with_holes


def is_convex(polygon: shapely.Polygon):
    coords = polygon.exterior.coords
    num_coords = len(coords)

    if num_coords < 4:
        # A polygon with less than 4 vertices cannot be convex
        return False

    # Calculate the orientation of the first three points
    orientation = 0
    for i in range(num_coords):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % num_coords]
        x3, y3 = coords[(i + 2) % num_coords]

        # Calculate the cross product of the vectors (x2-x1, y2-y1) and (x3-x2, y3-y2)
        cross_product = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)

        if cross_product != 0:
            orientation = cross_product
            break

    # Check the orientation of the remaining vertices
    for i in range(num_coords):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % num_coords]
        x3, y3 = coords[(i + 2) % num_coords]

        cross_product = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)

        if cross_product * orientation < 0:
            return False

    return True


def shapely_polygon_to_cgal(polygon: shapely.Polygon):
    # Exctract all the points except the last, as this is the same as the first
    cgal_points = [
        bindings.Point_2(point[0], point[1]) for point in polygon.exterior.coords[:-1]
    ]
    # Create a CGAL Polygon_2 from the list of points
    return bindings.Polygon_2(cgal_points)


def get_sweep_offset(overlap=0.1, height=10, field_of_view=90):
    if overlap < 0 or overlap > 1:
        msg = "Overlap percentage has to be a float between 0 and 1!"
        raise ValueError(msg)

    return abs(
        2 * height * math.tan((field_of_view * math.pi / 180.0) / 2) * (1 - overlap)
    )


def generate_sweep_pattern(
    polygon: shapely.Polygon,
    sweep_offset,
    clockwise=True,
    connect_sweeps=False,
):
    # Make sure that the orientation of the polygon is counterclockwise and the interior is clockwise
    cgal_poly = shapely_polygon_to_cgal(orient(polygon=polygon))
    segments = bindings.generate_sweeps(
        cgal_poly, sweep_offset, clockwise, connect_sweeps
    )

    if connect_sweeps:
        # Combine all segments into a single LineString
        combined_line = []
        for seg in segments:
            combined_line.append([seg.source.x, seg.source.y])
            combined_line.append([seg.target.x, seg.target.y])
        result = [shapely.LineString(combined_line)]
    else:
        lines = [
            shapely.LineString(
                [[seg.source.x, seg.source.y], [seg.target.x, seg.target.y]]
            )
            for seg in segments
        ]
        result = lines

    return result


def decompose_polygon(
    boundary: shapely.Polygon, obstacles: shapely.MultiPolygon | shapely.Polygon = None
):
    if obstacles is not None:
        if isinstance(obstacles, shapely.Polygon):
            obstacles = shapely.MultiPolygon([obstacles])
        elif not isinstance(obstacles, shapely.MultiPolygon):
            msg = "Obstacles must be a Shapely MultiPolygon."
            raise ValueError(msg)

        # If the obstacles intersect with the boundary, take the union of the two and remove it from the obstacles list
        updated_obstacles = []
        for obstacle in obstacles.geoms:
            if obstacle.intersects(boundary.boundary):
                log.debug(
                    "Obstacles intersect with the boundary, the geometries will be merged."
                )
                boundary = obstacles.union(boundary)
            else:
                updated_obstacles.append(obstacle)

        obstacles = shapely.MultiPolygon(updated_obstacles)
    pwh = bindings.Polygon_with_holes_2(shapely_polygon_to_cgal(boundary))
    if obstacles is not None:
        for poly in obstacles.geoms:
            pwh.add_hole(shapely_polygon_to_cgal(poly))
    decompose_polygons = bindings.decompose(pwh)
    return [
        shapely.Polygon([(vertex.x, vertex.y) for vertex in polygon])
        for polygon in decompose_polygons
    ]
