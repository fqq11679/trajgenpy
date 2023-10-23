import math

import pyproj
import shapely
from shapely.geometry.polygon import orient

import trajgenpy.bindings as bindings


class GeoData:
    def __init__(self, geometry, crs="WGS84"):
        self.set_geometry(geometry)
        # self.set_crs(crs)
        self.crs = crs

    def set_geometry(self, geometry):
        if not isinstance(
            geometry,
            shapely.LineString | shapely.Point | shapely.Polygon | shapely.MultiPolygon,
        ):
            msg = "Geometry must be a Shapely LineString, Point, or Polygon."
            raise ValueError(msg)
        self.geometry = geometry

    def set_crs(self, crs):
        if not isinstance(crs, str):
            raise ValueError("New CRS must be a string.")

        if crs != self.crs:
            # Apply the transformer to the geometry
            self.convert_to_crs(crs)
        self.crs = crs

    def convert_to_crs(self, crs):
        raise NotImplementedError(
            "This method sould be implemented in the data classes!"
        )

    def get_geometry(self):
        return self.geometry


class Trajectory(GeoData):
    def __init__(self, geometry, crs="WGS84"):
        super().__init__(geometry, crs)

    def convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)

        converted_coords = [
            transformer.transform(x, y) for x, y in list(self.geometry.coords)
        ]
        self.geometry = shapely.LineString(converted_coords)


class Point(GeoData):
    def __init__(self, geometry, crs="WGS84"):
        super().__init__(geometry, crs)

    def convert_to_crs(self, crs):
        transformer = pyproj.Transformer.from_crs(self.crs, crs, always_xy=True)
        x, y = transformer.transform(self.geometry.x, self.geometry.y)
        self.geometry = shapely.Point(x, y)


class Polygon(GeoData):
    def __init__(self, geometry, crs="WGS84"):
        super().__init__(geometry, crs)

    def convert_to_crs(self, crs):
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
    cgal_polygon = bindings.Polygon_2(cgal_points)
    return cgal_polygon


def get_sweep_offset(overlap=0.1, height=30, field_of_view=90):
    if overlap < 0 or overlap > 1:
        raise ValueError("Overlap percentage has to be a float between 0 and 1!")

    return abs(
        2 * height * math.tan((field_of_view * math.pi / 180.0) / 2) * (1 - overlap)
    )


def generate_sweep_pattern(
    polygon: shapely.Polygon,
    sweep_offset,
    direction=None,
    clockwise=True,
    disconnected_sweeps=True,
):
    # Make sure that the orientation of the polygon is counterclockwise and the interior is clockwise
    cgal_poly = shapely_polygon_to_cgal(orient(polygon=polygon))
    segments = bindings.generate_sweeps(cgal_poly, sweep_offset)

    if disconnected_sweeps:
        lines = [
            shapely.LineString(
                [[seg.source.x, seg.source.y], [seg.target.x, seg.target.y]]
            )
            for seg in segments
        ]
        result = lines
    else:
        # Combine all segments into a single LineString
        combined_line = []
        for seg in segments:
            combined_line.append([seg.source.x, seg.source.y])
            combined_line.append([seg.target.x, seg.target.y])
        result = [combined_line]
    return shapely.MultiLineString(result)


def decompose_polygon(boundary: shapely.Polygon, obstacles: shapely.MultiPolygon):
    pwh = bindings.PolygonWithHoles(shapely_polygon_to_cgal(boundary))

    for poly in obstacles:
        pwh.add_hole(shapely_polygon_to_cgal(poly))
    return bindings.decompose(pwh)
