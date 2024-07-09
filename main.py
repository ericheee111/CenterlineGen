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
import networkx as nx
from scipy.interpolate import splprep, splev
# import matplotlib.animation as animation
import geopandas as gpd

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

# def interpolate_points_with_density(lines, closest_points, xMostLabel, xBottomLabel):
#     # Get the last point of the top left line and the first point of the bottom left line
#     top_left_last_point = lines[xMostLabel][-1]
#     bottom_left_first_point = lines[closest_points[xBottomLabel][0]][0]

#     # Calculate the density of each line to determine interpolation points
#     density_top = calculate_density(lines[xMostLabel])
#     density_bottom = calculate_density(lines[closest_points[xBottomLabel][0]])
#     average_density = (density_top + density_bottom) / 2

#     # Determine number of interpolated points based on average density and distance between the points
#     distance_between_lines = np.linalg.norm(top_left_last_point - bottom_left_first_point)
#     num_interpolated_points = max(int(distance_between_lines / average_density), 1)  # At least one point
    
#     # Create interpolated points
#     interpolated_points = interpolate_points(top_left_last_point, bottom_left_first_point, num_interpolated_points)
#     return interpolated_points

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

## Main Starts Here
fig = plt.figure()
ax = fig.add_subplot(111)

lane = np.genfromtxt("data3.csv", delimiter=",")
lane = np.vstack([lane])

st = perf_counter()

x, y, x1, y1 = bounding_box(lane)
print(x, y, x1, y1)

lines = clustering_DBSCAN(lane)

# Find the closest points to the corners
leftmost_label, rightmost_label = cloest_point_to_corner(x, y, x1, y1, lines)

all_lines = [np.array(lines[label]) for label in lines.keys()]
sorted_lines = [line[np.argsort(line[:, 1])] for line in all_lines]

# Apply the smoothing function to each line
smoothed_lines = [smooth_line(line) for line in sorted_lines]

# polygons = [Polygon(line) for line in smoothed_lines]
polygons = [LineString(line) for line in smoothed_lines]

# sort lines[0] by y coordinate and store in lane0
leftmostlane = np.array(lines[leftmost_label])
leftmostlane = leftmostlane[np.argsort(leftmostlane[:,1])]
rightmostlane = np.array(lines[rightmost_label])
rightmostlane = rightmostlane[np.argsort(-rightmostlane[:,1])]

# combine leftmost and rightmost
poly = np.vstack([leftmostlane, rightmostlane])
# plt.plot(poly[:,0], poly[:,1])
bdry = Polygon(poly)

# mtpolygon = geometry.MultiPolygon(polygons)
mtpolygon = geometry.MultiLineString(polygons)

et = perf_counter()
print("Elapsed Time for Pre Processing: {:.10f}".format((et-st)*10**3), " ms")

colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

st = perf_counter()
centerline = Centerline(input_geometry=mtpolygon, bdry=bdry)
et = perf_counter()
print("Elapsed Time for Calculate & Post Processing Voronoi: {:.10f}".format((et-st)*10**3), " ms")

ctl = centerline.centerlines

coords = [cor.xy for cor in ctl]

for i, cor in enumerate(coords):
    # ax.plot(cor[0], cor[1], color=colors[i % len(colors)], marker='.', linewidth=1)
    ax.plot(cor[0], cor[1], color='b', marker='.', linewidth=1)

for i, line in enumerate(smoothed_lines):
    # print("# L - ", len(line))
    x, y = line.T
    ax.plot(x, y)

# def init():
#     ax.clear()
#     return[]

# def animate(num):
#     ax.clear()
#     ax.axis('equal')
#     for i in range(num):
#         cor = coords[i]
#         plt.plot(cor[0], cor[1], color=colors[i % len(colors)], marker='.', linewidth=1)
#     return[]

# ani = animation.FuncAnimation(fig, animate, frames=len(coords), init_func=init, blit=True, repeat=True)


## find branches
# st = perf_counter()

# Create a GeoDataFrame from LineStrings
# gdf = gpd.GeoDataFrame(geometry=ctl)
# sindex = gdf.sindex

# G = nx.Graph()

# # Add edges and nodes to the graph with intersection check optimization
# for idx, line in gdf.iterrows():
#     start, end = map(tuple, [line.geometry.coords[0], line.geometry.coords[-1]])
#     G.add_edge(start, end, line=line)

#     # Use spatial index to find potential intersections
#     possible_matches_index = list(sindex.intersection(line.geometry.bounds))
#     for other_idx in possible_matches_index:
#         if other_idx == idx:
#             continue
#         other_line = gdf.iloc[other_idx].geometry
#         if line.geometry.intersects(other_line):
#             intersection = line.geometry.intersection(other_line)
#             # Ensure that the intersection is added as a node
#             if 'Point' == intersection.geom_type:
#                 G.add_node(tuple(intersection.coords[0]))
#             elif 'MultiPoint' == intersection.type:
#                 for pt in intersection.geoms:
#                     G.add_node(tuple(pt.coords[0]))


# # Find branches
# branch_nodes = [node for node, degree in G.degree() if degree > 2]

# et = perf_counter()
# print("Elapsed Time for Find Branch: {:.10f}".format((et-st)*10**3), " ms")

# Optionally, extract subgraphs that contain these branch points
# print("Num of Branches: ", len(branch_nodes))
# for node in branch_nodes:
#     # Get all edges connected to this node
#     connected_edges = G.edges(node, data=True)
#     print(f"Branch at node {node}:")
    # for edge in connected_edges:
    #     print(f" ------------ Connects to {edge[1]} via {edge[2]['line']}")

# Highlight branch nodes in red
# for point in branch_nodes:
#     ax.scatter(point[0], point[1], facecolors='none', edgecolors='red', s=100)

## identify each line by branch nodes
st = perf_counter()
multiline = geometry.MultiLineString(ctl)
# print(multiline)
merged = linemerge(multiline)
# print(merged)
et = perf_counter()
print("Elapsed Time for Identify Lanes: {:.10f}".format((et-st)*10**3), " ms")

for i, line in enumerate(merged.geoms):
    x, y = line.xy
    ax.plot(x, y, color=colors[i % len(colors)], marker='.', linewidth=1)

def read_lines_from_file(filename):
    lines = []
    with open(filename, 'r') as file:
        for line in file:
            points = list(map(float, line.replace(',', ' ').split()))
            lines.append(points)
    return lines

def parse_points(lines):
    parsed_lines = []
    for points in lines:
        x = points[0::2]
        y = points[1::2]
        parsed_lines.append((x, y))
    return parsed_lines

filename = 'data.txt'
lines = read_lines_from_file(filename)
parsed_lines = parse_points(lines)
for line in parsed_lines:
    x, y = line
    ax.plot(x, y, color='black', linewidth=1)



ax.axis('equal')
plt.show()





