import logging

import shapely
import shapely.plotting as shplt
from shapely.affinity import rotate

from matplotlib import pyplot as plt
from trajgenpy import Geometries, Logging
import geojson
import numpy as np

log = Logging.get_logger()


def custom_example():
    logging.warning("This is a warning message")
    poly = shapely.Polygon(
        [
            [10.470650724002098, 55.38520996875883],
            [10.469291130631092, 55.38477524425895],
            [10.468346742421403, 55.3855938172706],
            [10.466164880005067, 55.384530131529345],
            [10.46587993528641, 55.38429426690277],
            [10.46378762692467, 55.381052141937516],
            [10.46446335297145, 55.38091801121436],
            [10.465041383686241, 55.380941137233464],
            [10.465733392288826, 55.3808902599732],
            [10.466262575337083, 55.3807653791485],
            [10.467736146595996, 55.38083475743292],
            [10.469681911959213, 55.38064974840299],
            [10.470015704344348, 55.381177021857354],
            [10.4717823615984, 55.38068212504581],
            [10.472189425482355, 55.380964263239036],
            [10.473744409517451, 55.38081625656886],
            [10.470553028670054, 55.381861541821166],
            [10.469535368961601, 55.38115852115337],
            [10.467695440207848, 55.381491532500746],
            [10.469193435299047, 55.38253217489091],
            [10.467003431605235, 55.3834109182535],
            [10.468371166253775, 55.384391387802054],
            [10.47115548321807, 55.38393815489968],
            [10.471619536045012, 55.384285017281144],
            [10.470650724002098, 55.38520996875883],
        ]
    )


    #----------------------------------测试数据----------------------------
    def transform_points(ox, oy, polygon):
        # Calculate the min and max of the original polygon's x (lon) and y (lat) coords
        poly_x = [point[0] for point in polygon.exterior.coords]
        poly_y = [point[1] for point in polygon.exterior.coords]
        poly_x_min, poly_x_max = min(poly_x), max(poly_x)
        poly_y_min, poly_y_max = min(poly_y), max(poly_y)

        # Calculate the min and max of the ox and oy coords
        ox_min, ox_max = min(ox), max(ox)
        oy_min, oy_max = min(oy), max(oy)

        # Transform the ox and oy values to the new coordinate system
        transformed_ox = [np.interp(x, (ox_min, ox_max), (poly_x_min, poly_x_max)) for x in ox]
        transformed_oy = [np.interp(y, (oy_min, oy_max), (poly_y_min, poly_y_max)) for y in oy]

        # Create a list of tuples from the transformed coordinates
        transformed_points = list(zip(transformed_ox, transformed_oy))
    
        # Create a new Polygon object from the transformed coordinates
        transformed_polygon = shapely.Polygon(transformed_points)

        return transformed_polygon

    # Use the defined function to transform the ox and oy lists and create a new polygon
    
    #case1
    #ox = [0.0, 20.0, 50.0, 100.0, 130.0, 40.0, 0.0]
    #oy = [0.0, -20.0, 0.0, 30.0, 60.0, 80.0, 0.0]

    #case2
    #ox = [0.0, 20.0, 50.0, 100.0, 130.0, 40.0, 0.0]
    #oy = [0.0, -20.0, 0.0, 0.0, 60.0, 80.0, 0.0]
   
    #case3
    #ox = [0.0, 50.0, 50.0, 0.0, 0.0]
    #oy = [0.0, 0.0, 30.0, 30.0, 0.0]

    #case4
    #ox = [0.0, 20.0, 50.0, 200.0, 130.0, 40.0, 0.0]
    #oy = [0.0, -80.0, 0.0, 30.0, 60.0, 80.0, 0.0]

    #case5
    #ox = [0.0, 200.0, 200.0, 100.0, 0.0, 0.0]
    #oy = [0.0, 0.0, 200.0, 20.0, 200.0, 0.0]

    #case6
    #ox = [50, 100, 150, 160, 100, 50, 30, 30, 50]
    #oy = [50, 30, 50, 100, 150, 160, 150, 50, 50]

    #case7
    #ox = [50, 100, 150, 160, 100, 50, 30, 30, 50]
    #oy = [50, 30, 50, 100, 150, 160, 150, 50, 50]
    #ox = ox[::-1] # 反转列表
    #oy = oy[::-1]

    #transformed_polygon = transform_points(ox, oy, poly)
    #poly = transformed_polygon
    #-----------------------------------------------------------------------
    # 设置旋转角
    alpha = 45

    #-----------------------------------------------------------------------

    geo_poly = Geometries.GeoPolygon(poly)
    geo_poly.set_crs("EPSG:3857")
    
    geo_poly.geometry = rotate(geo_poly.get_geometry(), angle=alpha, origin=(0, 0))

    polygon_list = Geometries.decompose_polygon(geo_poly.get_geometry(), obstacles=None)
    # Define the coordinates for the inner polygon
    hole = shapely.Polygon(
        [
            [10.465643258937206, 55.38227460993366],
            [10.464960949768766, 55.38131078946813],
            [10.466980216090207, 55.380975542061805],
            [10.46771784762305, 55.38200746639646],
            [10.465643258937206, 55.38227460993366],
        ]
    )

    hole = Geometries.GeoPolygon(hole)
    hole.set_crs("EPSG:3857")

    polygon_list = Geometries.decompose_polygon(
        geo_poly.get_geometry()#, obstacles=hole.get_geometry()
    )

    offset = Geometries.get_sweep_offset(0.1, 30, 100)
    print(offset)
    result = []
    for decomposed_poly in polygon_list:
        sweeps_connected = Geometries.generate_sweep_pattern(
            decomposed_poly, offset, clockwise=True, connect_sweeps=True
        )
        assert len(sweeps_connected) == 1

        sweeps_disconnected = Geometries.generate_sweep_pattern(
            decomposed_poly, offset, clockwise=False, connect_sweeps=True
        )

        # Plotting
        result.extend(sweeps_disconnected)
    multi_traj = Geometries.GeoMultiTrajectory(result, "EPSG:3857")

    # # Export the polygon and the buildings as a GeoJSON file
    # geojson_collection = geojson.FeatureCollection(
    #     [
    #         geo_poly.to_geojson(id="boundary"),
    #         hole.to_geojson(id="obstacles"),
    #         multi_traj.to_geojson(id="tasks"),
    #     ]
    # )
    # with open("environment.geojson", "w") as f:
    #     geojson.dump(geojson_collection, f)

    # Convert the geometry back to WGS84(geodesic)
    # multi_traj.set_crs("WGS84")
    # geo_poly.set_crs("WGS84")
    # hole.set_crs("WGS84")
    
    #for decomposed_poly in polygon_list:
    #    d_poly = Geometries.GeoPolygon(decomposed_poly)
    #    d_poly.plot(facecolor="grey", edgecolor="green",linewidth=2, alpha=0.4)
    
    #hole.plot(facecolor="red", edgecolor="red", linewidth=2, alpha=0.2)
    multi_traj.concatenate_trajectories(geo_poly.get_geometry())

    geo_poly.geometry = rotate(geo_poly.get_geometry(), angle=-alpha, origin=(0, 0))
    geo_poly.plot(facecolor="grey", linewidth=2, alpha=0.4)
    multi_traj.geometry = rotate(multi_traj.get_geometry(), angle=-alpha, origin=(0, 0))
    multi_traj.plot(color="black", linewidth=1)
    plt.axis("equal")
    plt.axis("off")
    plt.savefig("output.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    custom_example()
