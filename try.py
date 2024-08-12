from numpy import array
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
from time import perf_counter
import numpy as np

borders = [np.array([1.5484909, -0.21024238]),
                np.array([1.5413352, -0.21024236]),
                np.array([1.5378445, -0.21024236]),
                np.array([1.5343539, -0.21024235]),
                np.array([1.5308806, -0.19333011])]
voronoi_diagram = Voronoi(borders)
vertices = voronoi_diagram.vertices
ridges = voronoi_diagram.ridge_vertices

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
voronoi_plot_2d(voronoi_diagram, ax=ax)

# Add the borders as a LineString
border_line = LineString(borders)
ax.plot(*border_line.xy, color='red', linewidth=2)

# Set the aspect ratio and limits of the plot
ax.set_aspect('equal')
ax.set_xlim(border_line.bounds[0] - 0.1, border_line.bounds[2] + 0.1)
ax.set_ylim(border_line.bounds[1] - 0.1, border_line.bounds[3] + 0.1)

# Show the plot
plt.show()
