import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import shapely.geometry as geometry
from shapely.geometry import Point, LineString, GeometryCollection
from shapely.ops import split, unary_union, linemerge
import numpy as np

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

def plot_lines(lines):
    for line in lines:
        x, y = line
        plt.plot(x, y)

    plt.axis('equal')
    plt.show()

def convert_to_linestring(parsed_lines):
    linestrings = []
    for line in parsed_lines:
        points = list(zip(line[0], line[1]))
        linestrings.append(LineString(points))
    return linestrings

if __name__ == "__main__":
    filename = 'data.txt'
    lines = read_lines_from_file(filename)
    parsed_lines = parse_points(lines)
    print(len(parsed_lines))
    plot_lines(parsed_lines)

    linestrings = convert_to_linestring(parsed_lines)        
        # Do something with the line string

    multiline = geometry.MultiLineString(linestrings)
    # print(len(multiline))
    merged = linemerge(multiline)
    # print(merged)
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    for i, line in enumerate(merged.geoms):
        x, y = line.xy
        plt.plot(x, y, color=colors[i % len(colors)], marker='.', linewidth=1)

    plt.axis('equal')
    plt.show()
