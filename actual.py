import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from time import perf_counter
# from scipy.spatial import Voronoi, _adjust_bounds, voronoi_plot_2d
from sklearn.cluster import DBSCAN
# from scipy.spatial import ConvexHull
# from scipy.spatial import cKDTree
from shapely.geometry import Polygon
import shapely.geometry as geometry
from shapely.geometry import Point, LineString, GeometryCollection
from shapely.ops import split, unary_union, linemerge
from geometry import *
from exceptions import *
from scipy.interpolate import splprep, splev
# import matplotlib.animation as animation
import alphashape

def interpolate_points(start, end, num_points):
    # Linearly interpolate between two points in 2D space
    x_vals = np.linspace(start[0], end[0], num_points)
    y_vals = np.linspace(start[1], end[1], num_points)
    return np.column_stack((x_vals, y_vals))

def calculate_density(line):
    # Calculate distances between consecutive points in the line
    distances = np.sqrt(np.sum(np.diff(line, axis=0)**2, axis=1))
    # Return the average distance
    return np.mean(distances)

def cloest_point_to_corner(x, y, x1, y1, lines):
    topleft = (x, y1)
    topright = (x1, y1)
    bottomleft = (x, y)
    bottomright = (x1, y)
    closest_points = {}
    for corner in [topleft, topright, bottomleft, bottomright]:
        closest_points[corner] = None
        min_distance = float('inf')
        for label, points in lines.items():
            for point in points:
                distance = np.linalg.norm(np.array(corner) - np.array(point))
                if distance < min_distance:
                    min_distance = distance
                    closest_points[corner] = (label, point)

    leftmost_label = closest_points[topleft][0]
    rightmost_label = closest_points[topright][0]

    if closest_points[topleft][0] != closest_points[bottomleft][0]:
        lines[leftmost_label] = np.vstack([lines[leftmost_label], lines[closest_points[bottomleft][0]]])
        del lines[closest_points[bottomleft][0]]
    if closest_points[topright][0] != closest_points[bottomright][0]:
        lines[rightmost_label] = np.vstack([lines[rightmost_label], lines[closest_points[bottomright][0]]])
        del lines[closest_points[bottomright][0]]

    return leftmost_label, rightmost_label

def clustering_DBSCAN(lane):
    clustering = DBSCAN(eps=2, min_samples=2).fit(lane)
    labels = clustering.labels_

    lines = {}
    for label, point in zip(labels, lane):
        if label in lines:
            lines[label].append(point)
        else:
            lines[label] = [point]
    return lines

def bounding_box(points):
    x_coordinates, y_coordinates = zip(*points)

    return min(x_coordinates), min(y_coordinates), max(x_coordinates), max(y_coordinates)

def smooth_line(line):
    # Separate the x and y coordinates
    x, y = line.T

    num_points = int((max(y) - min(y)) * 1.5)

    
    # Perform spline interpolation
    tck, u = splprep([x, y], s=1)
    unew = np.linspace(0, 1, num_points)
    out = splev(unew, tck)
    
    # Combine the new x and y coordinates
    smoothed_line = np.vstack(out).T
    return smoothed_line

def custom_csv_reader(file_path):
    data = defaultdict(list)
    
    with open(file_path, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            if len(values) < 3:
                continue  # Skip lines with insufficient data
            
            try:
                timestamp = values[0]
                lane = values[1]
                coordinates = [(float(values[i]), float(values[i + 1])) for i in range(2, len(values) - 1, 2)]
                data[timestamp].append(coordinates)
            except ValueError:
                continue  # Skip lines with conversion errors
    
    return data

# Read and parse the CSV file
file_path = 'lanes.csv'
parsed_data = custom_csv_reader(file_path)
count = 10
for timestamp, smoothed_lines in parsed_data.items():
    if count != 10:
        count += 1
        continue
    # print(len(smoothed_lines))
    # for line in smoothed_lines:
    #     print(len(line[1]))

    # break
    # Create a new figure
    fig, ax = plt.subplots()
    # polygons = [Polygon(line) for line in smoothed_lines]
    polygons = [LineString(line) for line in smoothed_lines]

    # # sort lines[0] by y coordinate and store in lane0
    # leftmostlane = np.array(lines[leftmost_label])
    # leftmostlane = leftmostlane[np.argsort(leftmostlane[:,1])]
    # rightmostlane = np.array(lines[rightmost_label])
    # rightmostlane = rightmostlane[np.argsort(-rightmostlane[:,1])]

    # # combine leftmost and rightmost
    # poly = np.vstack([leftmostlane, rightmostlane])
    # plt.plot(poly[:,0], poly[:,1])
    # bdry = Polygon(poly)

    # Calculate the alpha shape boundary
    alpha = 0.2  # Adjust the alpha value as needed
    alphalines = []
    for line in smoothed_lines:
        alphalines += line
    bdry = alphashape.alphashape(alphalines, alpha)
    
    for line in smoothed_lines:
        x, y = np.array(line).T
        ax.plot(x, y)
    
    ax.plot(*bdry.exterior.xy, color='red')
    plt.show()

    fig, ax = plt.subplots()
    # mtpolygon = geometry.MultiPolygon(polygons)
    mtpolygon = geometry.MultiLineString(polygons)

    # et = perf_counter()
    # print("Elapsed Time for Pre Processing: {:.10f}".format((et-st)*10**3), " ms")

    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

    # st = perf_counter()
    centerline = Centerline(input_geometry=mtpolygon, bdry=bdry)
    # et = perf_counter()
    # print("Elapsed Time for Calculate & Post Processing Voronoi: {:.10f}".format((et-st)*10**3), " ms")

    ctl = centerline.centerlines

    print(len(ctl))

    coords = [cor.xy for cor in ctl]

    for i, cor in enumerate(coords):
        # ax.plot(cor[0], cor[1], color=colors[i % len(colors)], marker='.', linewidth=1)
        ax.plot(cor[0], cor[1], color='b', marker='.', linewidth=1)

    for i, line in enumerate(smoothed_lines):
        # print("# L - ", len(line))
        x, y = np.array(line).T
        ax.plot(x, y)
        
    ## identify each line by branch nodes
    # st = perf_counter()
    multiline = geometry.MultiLineString(ctl)
    # print(len(multiline))
    merged = linemerge(multiline)
    # print(merged)
    # et = perf_counter()
    # print("Elapsed Time for Identify Lanes: {:.10f}".format((et-st)*10**3), " ms")

    for i, line in enumerate(merged.geoms):
        x, y = line.xy
        ax.plot(x, y, color=colors[i % len(colors)], marker='.', linewidth=1)

    ax.axis('equal')
    plt.show()
    break




