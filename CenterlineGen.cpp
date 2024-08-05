#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <unordered_set>
#include <chrono>

#define JC_VORONOI_IMPLEMENTATION
#define JC_VORONOI_CLIP_IMPLEMENTATION
//#define OUTPIC

#include "DBSCAN.h"
#include "jc_voronoi.h"
#include "jc_voronoi_clip.h"
#include "stb_image_write.h"
#include "strtree.hpp"

#include "PolynomialRegression.h"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/index/predicates.hpp>
#include <boost/geometry/index/adaptors/query.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/index/rtree.hpp>
//#include <boost/timer/timer.hpp>
#include <boost/geometry/algorithms/buffer.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Simple_cartesian.h>

#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <cassert>



typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> bst_point;
typedef boost::geometry::model::segment<bst_point> bst_segment;
typedef boost::geometry::model::linestring<bst_point> bst_linestring;
typedef boost::geometry::model::polygon<bst_point> bst_polygon;
typedef std::pair<bst_segment, uint32_t> bst_value;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::FT                                                FT;
typedef K::Point_2                                           CGAL_Point;
typedef K::Segment_2                                         CGAL_Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>          Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;

typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;

typedef struct VoronoiContext_ {
    uint32_t num_points;
    jcv_point* points;
    jcv_rect rect;
    jcv_diagram diagram;
    jcv_clipper clipper;
    jcv_clipper* clipper_ptr;
} vorcon;

std::map<long long, std::vector<std::vector<PointDB>>> readActualCSV(const std::string& filename) {
    std::cout << "read file" << std::endl;
    std::ifstream file(filename);
    std::map<long long, std::vector<std::vector<PointDB>>> data;

    if (!file.is_open()) {
        std::cerr << "Unable to open file" << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;

        long long timestamp;
        std::vector<std::vector<PointDB>> lanes;

        try {
            // Read timestamp
            if (!std::getline(ss, cell, ',')) throw std::invalid_argument("Missing timestamp");
            timestamp = std::stoll(cell);
            if (!std::getline(ss, cell, ',')) throw std::invalid_argument("Missing lane number");
            int laneNumber = std::stoi(cell);

            // Read lane coordinates
            while (std::getline(ss, cell, ',')) {
                std::vector<PointDB> lane;

                do {
                    PointDB coordinate;

                    // Read x coordinate
                    coordinate.x = std::stod(cell);
                    coordinate.clusterID = laneNumber;

                    // Check if next cell is available for y coordinate
                    if (!std::getline(ss, cell, ',')) break;
                    coordinate.y = std::stod(cell);

                    lane.push_back(coordinate);
                } while (std::getline(ss, cell, ','));

                lanes.push_back(lane);
            }
            data[timestamp].insert(data[timestamp].end(), lanes.begin(), lanes.end());
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing line: " << line << std::endl;
            std::cerr << "Exception: " << e.what() << std::endl;
            // Optionally, continue processing other lines or stop here
            continue;
        }
    }

    file.close();
    return data;
}

/* ----------------------------------------------------------------------- */
static void plot(int x, int y, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    if (x < 0 || y < 0 || x >(width - 1) || y >(height - 1))
        return;
    int index = y * width * nchannels + x * nchannels;
    for (int i = 0; i < nchannels; ++i)
    {
        image[index + i] = color[i];
    }
}

// http://members.chello.at/~easyfilter/bresenham.html
static void draw_line(int x0, int y0, int x1, int y1, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; // error value e_xy

    for (;;)
    {  // loop
        plot(x0, y0, image, width, height, nchannels, color);
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; } // e_xy+e_x > 0
        if (e2 <= dx) { err += dx; y0 += sy; } // e_xy+e_y < 0
    }
}

// http://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
static inline int orient2d(const jcv_point* a, const jcv_point* b, const jcv_point* c)
{
    return ((int)b->x - (int)a->x) * ((int)c->y - (int)a->y) - ((int)b->y - (int)a->y) * ((int)c->x - (int)a->x);
}

static inline int min2(int a, int b)
{
    return (a < b) ? a : b;
}

static inline int max2(int a, int b)
{
    return (a > b) ? a : b;
}

static inline int min3(int a, int b, int c)
{
    return min2(a, min2(b, c));
}
static inline int max3(int a, int b, int c)
{
    return max2(a, max2(b, c));
}

// Remaps the point from the input space to image space
static inline jcv_point remap(const jcv_point* pt, const jcv_point* min, const jcv_point* max, const jcv_point* scale)
{
    jcv_point p;
    p.x = (pt->x - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static inline jcv_point remap(const PointDB* pt, const jcv_point* min, const jcv_point* max, const jcv_point* scale)
{
    jcv_point p;
    p.x = (pt->x - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static inline jcv_point remap(const CGAL_Point* pt, const jcv_point* min, const jcv_point* max, const jcv_point* scale)
{
    jcv_point p;
    p.x = (pt->x() - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y() - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static void draw_triangle(const jcv_point* v0, const jcv_point* v1, const jcv_point* v2, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    int area = orient2d(v0, v1, v2);
    if (area == 0)
        return;

    // Compute triangle bounding box
    int minX = min3((int)v0->x, (int)v1->x, (int)v2->x);
    int minY = min3((int)v0->y, (int)v1->y, (int)v2->y);
    int maxX = max3((int)v0->x, (int)v1->x, (int)v2->x);
    int maxY = max3((int)v0->y, (int)v1->y, (int)v2->y);

    // Clip against screen bounds
    minX = max2(minX, 0);
    minY = max2(minY, 0);
    maxX = min2(maxX, width - 1);
    maxY = min2(maxY, height - 1);

    // Rasterize
    jcv_point p;
    for (p.y = (jcv_real)minY; p.y <= (jcv_real)maxY; p.y++) {
        for (p.x = (jcv_real)minX; p.x <= (jcv_real)maxX; p.x++) {
            // Determine barycentric coordinates
            int w0 = orient2d(v1, v2, &p);
            int w1 = orient2d(v2, v0, &p);
            int w2 = orient2d(v0, v1, &p);

            // If p is on or inside all edges, render pixel.
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                plot((int)p.x, (int)p.y, image, width, height, nchannels, color);
            }
        }
    }
}

/* ----------------------------------------------------------------------- */
std::vector<PointDB> readCSV(const std::string& filename) {
    std::vector<PointDB> points;
    std::ifstream file(filename);
    std::string line, value;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return points;
    }

    while (getline(file, line)) {
        std::stringstream ss(line);
        double x, y;

        // Read the x value
        getline(ss, value, ',');
        x = std::stof(value);

        // Read the y value
        getline(ss, value, ',');
        y = std::stof(value);

        // Create a Point and add to the vector
        points.emplace_back(x, y);
    }

    file.close();
    return points;
}

std::vector<std::pair<double, double>> findBoundingBox(const std::vector<PointDB>& points) {
    // Initialize min and max values
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::min();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::min();

    // Iterate through each point
    for (const auto& point : points) {
        double x = point.x;
        double y = point.y;

        // Update min and max values
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    // Return the bounding box as a vector of pairs
    return { {minX, minY}, {maxX, maxY} };
}

std::vector<std::pair<double, double>> findBoundingBox(const std::vector<std::vector<PointDB>>& lanes) {
    // Initialize min and max values
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::min();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::min();

    // Iterate through each point
    for (const auto& points : lanes) {
        for (const auto& point : points) {
            double x = point.x;
            double y = point.y;

            // Update min and max values
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
    }
    

    // Return the bounding box as a vector of pairs
    return { {minX, minY}, {maxX, maxY} };
}

std::vector<CGAL_Point> convert_to_cgal_points(const std::vector<std::vector<PointDB>>& lane_lines, vorcon& vc) {
	size_t num_points = 0;
	for (const auto& line : lane_lines) {
		num_points += line.size();
	}

    std::vector<CGAL_Point> cgalpoints;

	jcv_point* points = new jcv_point[num_points];
	size_t idx = 0;
	for (const auto& line : lane_lines) {
		for (const auto& point : line) {
			points[idx] = { point.x, point.y };
            cgalpoints.push_back(CGAL_Point(point.x, point.y));
            idx++;
		}
	}
    vc.num_points = num_points;
    vc.points = points;
    return cgalpoints;
}

template <class OutputIterator>
void alpha_edges(const Alpha_shape_2& A, OutputIterator out) {
    for (Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(); it != A.alpha_shape_edges_end(); ++it) {
        *out++ = A.segment(*it);
    }
}

std::vector<jcv_point> segments_to_path(const std::vector<CGAL_Segment>& segments) {
    std::map<jcv_point, std::vector<jcv_point>> adjacency_list;
    std::set<std::pair<jcv_point, jcv_point>> visited_segments;

    // Construct adjacency list
    for (const auto& segment : segments) {
        jcv_point source = { segment.source().x(), segment.source().y() };
        jcv_point target = { segment.target().x(), segment.target().y() };
        adjacency_list[source].push_back(target);
        adjacency_list[target].push_back(source);
    }

    // Start from the first segment's source
    std::vector<jcv_point> path;
    jcv_point start = { segments.front().source().x(), segments.front().source().y() };
    jcv_point current = start;
    path.push_back(current);

    // To handle branches and backtracking
    std::vector<std::pair<jcv_point, size_t>> stack;

    while (true) {
        bool found_next = false;

        for (size_t i = 0; i < adjacency_list[current].size(); ++i) {
            jcv_point next = adjacency_list[current][i];
            if (visited_segments.find({ current, next }) == visited_segments.end()) {
                stack.push_back({ current, i });
                visited_segments.insert({ current, next });
                visited_segments.insert({ next, current });
                current = next;
                path.push_back(current);
                found_next = true;
                break;
            }
        }

        if (!found_next) {
            if (current == start) {
                break; // We have returned to the start, so the polygon is complete
            }
            if (stack.empty()) {
                path.clear(); // No more points to backtrack to and not closed, path is invalid
                break;
            }
            auto prev_point = stack.back().first;
            auto index = stack.back().second;
            stack.pop_back();
            current = prev_point;
            path.pop_back();
        }
    }

    if (!path.empty()) {
        path.push_back(start); // Close the path only if a valid path is found
    }

    return path;
}

//void testing() {
//    std::vector<std::pair<jcv_point, jcv_point>> centerlines11;
//    std::vector<jcv_point> boundary11 = {
//        {50, 50}, {250, 50}, {250, 250}, {50, 250}, {50, 50}
//    };
//
//    std::vector<std::vector<jcv_point>> lanelines = {
//        {{80, 50}, {80, 80}, {80, 120}, {80, 150}, {80, 250}},
//        {{120, 50}, {120, 80}, {120, 120}, {120, 150}, {120, 250}},
//        {{180, 50}, {180, 80}, {180, 120}, {180, 150}, {180, 250}}
//    };
//
//    for (int i = 0; i < 50; ++i) {
//        // Generate random points within the boundary
//        jcv_point start = { rand() % (149 - 51 + 1) + 51, rand() % (149 - 51 + 1) + 51 };
//        jcv_point end = { rand() % (149 - 51 + 1) + 51, rand() % (149 - 51 + 1) + 51 };
//
//        // Ensure the start point is not the same as the end point
//        while (start.x == end.x && start.y == end.y) {
//            end = { (double)(rand() % (149 - 51 + 1)) + 51, (double)(rand() % (149 - 51 + 1) + 51) };
//        }
//
//        centerlines11.push_back({ {start}, {end} });
//    }
//    for (int i = 0; i < 50; ++i) {
//        // Generate random points within the boundary
//        jcv_point start = { rand() % (75 - 30 + 1) + 30, rand() % (65 - 30 + 1) + 33 };
//        jcv_point end = { rand() % (90 - 39 + 1) + 40, rand() % (110 - 87 + 1) + 87 };
//        centerlines11.push_back({ {start}, {end} });
//    }
//
//    for (int i = 0; i < 50; ++i) {
//        // Generate random points within the boundary
//        jcv_point start = { rand() % (250 - 200 + 1) + 200, rand() % (250 - 170 + 1) + 180 };
//        jcv_point end = { rand() % (250 - 200 + 1) + 200, rand() % (250 - 170 + 1) + 180 };
//        centerlines11.push_back({ {start}, {end} });
//    }
//
//    for (int i = 0; i < 50; ++i) {
//        // Generate random points within the boundary
//        jcv_point start = { rand() % (200 - 130 + 1) + 130, rand() % (190 - 120 + 1) + 120 };
//        jcv_point end = { rand() % (200 - 130 + 1) + 130, rand() % (190 - 120 + 1) + 120 };
//        centerlines11.push_back({ {start}, {end} });
//    }
//
//
//    unsigned char color_red[] = { 255, 0, 0 };
//    unsigned char color_blue[] = { 0, 0, 255 };
//    unsigned char color_green[] = { 0, 255, 0 };
//    int width = 256;
//    int height = 256;
//    size_t imagesize = (size_t)(width * height * 3);
//    unsigned char* image = (unsigned char*)malloc(imagesize);
//    memset(image, 0, imagesize);
//    jcv_point dimensions;
//    dimensions.x = (jcv_real)width;
//    dimensions.y = (jcv_real)height;
//
//    for (int i = 1; i < boundary11.size(); i++) {
//        jcv_point min = { 0,0 };
//        jcv_point max = { 200, 200 };
//        /*auto pt1 = remap(&boundary11[i], &min, &max, &dimensions);
//        auto pt2 = remap(&boundary11[i - 1], &min, &max, &dimensions);*/
//        auto pt1 = boundary11[i];
//        auto pt2 = boundary11[i - 1];
//        draw_line(pt1.x, pt1.y, pt2.x, pt2.y, image, width, height, 3, color_blue);
//    }
//
//    for (auto& lane : lanelines) {
//        for (int i = 1; i < lane.size(); i++) {
//            auto pt1 = lane[i];
//            auto pt2 = lane[i - 1];
//            draw_line(pt1.x, pt1.y, pt2.x, pt2.y, image, width, height, 3, color_blue);
//        }
//    }
//
//
//    for (int i = 0; i < centerlines11.size(); i++) {
//        jcv_point min = { 0,0 };
//        jcv_point max = { 200, 200 };
//        /*auto pt1 = remap(&centerlines11[i].first, &min, &max, &dimensions);
//        auto pt2 = remap(&centerlines11[i].second, &min, &max, &dimensions);*/
//        auto pt1 = centerlines11[i].first;
//        auto pt2 = centerlines11[i].second;
//        draw_line(pt1.x, pt1.y, pt2.x, pt2.y, image, width, height, 3, color_red);
//    }
//
//        
//    auto filteredCenterlines = filterCenterlines(centerlines11, boundary11, lanelines);
//        
//    for (const auto& line : filteredCenterlines) {
//        //std::cout << "Filtered line: ((" << line.first.x << ", " << line.first.y << "), (" << line.second.x << ", " << line.second.y << "))\n";
//        jcv_point min = { 0,0 };
//        jcv_point max = { 200, 200 };
//        /*auto pt1 = remap(&line.first, &min, &max, &dimensions);
//        auto pt2 = remap(&line.second, &min, &max, &dimensions);*/
//        auto pt1 = line.first;
//        auto pt2 = line.second;
//
//        draw_line(pt1.x, pt1.y, pt2.x, pt2.y, image, width, height, 3, color_green);
//    }
//        
//
//    // flip image
//    int stride = width * 3;
//    uint8_t* row = (uint8_t*)malloc((size_t)stride);
//    for (int y = 0; y < height / 2; ++y)
//    {
//        memcpy(row, &image[y * stride], (size_t)stride);
//        memcpy(&image[y * stride], &image[(height - 1 - y) * stride], (size_t)stride);
//        memcpy(&image[(height - 1 - y) * stride], row, (size_t)stride);
//    }
//
//    char path[512];
//    sprintf_s(path, "geosstrtreeexample.png");
//
//    stbi_write_png(path, width, height, 3, image, stride);
//    std::cout << "done " << path << std::endl;
//
//    free(image);
//}

//void testing() {
//    // Initialize the GEOS library
//    geos::geom::GeometryFactory::Ptr factory = geos::geom::GeometryFactory::create();
//
//    // Example vector of vector of jcv_point for demonstration purposes
//    std::vector<std::vector<jcv_point>> lines = {
//        {{74.151,-13.1434}, {74.4412,-13.2563}},
//        {{74.4412,-13.2563}, {75.11,-13.3341}},
//        {{75.11,-13.3341}, {75.4286,-13.4584}},
//        {{75.4286,-13.4584}, {76.067,-13.5327}, {76.4145,-13.6692}},
//        {{86.8501,-15.9683}, {87.244,-16.132}},
//        {{87.244,-16.132}, {87.8278,-16.1999}},
//        {{87.8278,-16.1999}, {88.2279,-16.3668}, {88.8051,-16.4339}, {89.2115,-16.6037}, {89.7823,-16.6701}},
//        {{89.7823,-16.6701}, {90.1951,-16.843}}
//        
//    };
//
//    // Create a MultiLineString from the vector of lines
//    std::unique_ptr<geos::geom::MultiLineString> multiLine = makeMultiLineString(lines, factory.get());
//
//    // Use LineMerger to merge the MultiLineString
//    geos::operation::linemerge::LineMerger lineMerger;
//    lineMerger.add(multiLine.get());
//
//    // Get the merged result
//    auto merged(lineMerger.getMergedLineStrings());
//
//    // Print the result (for demonstration purposes)
//    for (const auto& line : merged) {
//		std::cout << line->toString() << std::endl;
//	}
//}

// Function to convert PointDB to jcv_point
inline jcv_point convertPointDBToJcvPoint(const PointDB& point) {
    return { point.x, point.y };
}

// Function to convert vector of vector of PointDB to vector of vector of jcv_point
std::vector<std::vector<jcv_point>> pointDBtoJCV(const std::vector<std::vector<PointDB>>& points) {
    std::vector<std::vector<jcv_point>> result;
    result.reserve(points.size()); // Reserve space to improve performance

    for (const auto& row : points) {
        std::vector<jcv_point> newRow;
        newRow.reserve(row.size()); // Reserve space to improve performance
        for (const auto& point : row) {
            newRow.push_back(convertPointDBToJcvPoint(point));
        }
        result.push_back(std::move(newRow));
    }

    return result;
}

void runactualgeos()
{
    std::string filename = "lanes.csv";
    std::map<long long, std::vector<std::vector<PointDB>>> lanedata = readActualCSV(filename);
    std::cout << "Number of timestamps: " << lanedata.size() << std::endl;
    int count = 0;
    uint32_t totaltime = 0;
    for (const auto& entry : lanedata) {
        /*if (count != 10) {
            count++;
            continue;
        }*/
        //std::cout << "Timestamp: " << entry.first << std::endl;
        std::vector<std::vector<PointDB>> inlanes = entry.second; // all lanes in one timestemp

        std::vector<std::vector<PointDB>> lanes;
        for (const auto& inlane : inlanes) {
            std::vector<PointDB> line;
            for (uint32_t i = 0; i < inlane.size(); i += 1) {
                line.push_back(inlane[i]);
            }
            lanes.push_back(line);
        }

        //std::cout << "Number of lanes: " << lanes.size() << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::pair<double, double>> bounding_box = findBoundingBox(lanes);

        auto lanelines = pointDBtoJCV(lanes);

        vorcon vc;
        int width = 512;
        int height = 512;
        size_t imagesize = (size_t)(width * height * 3);
        unsigned char* image = (unsigned char*)malloc(imagesize);
        memset(image, 0, imagesize);

        jcv_point dimensions;
        dimensions.x = (jcv_real)width;
        dimensions.y = (jcv_real)height;

        std::vector<CGAL_Point> points = convert_to_cgal_points(lanes, vc);
        vc.rect = { bounding_box[0].first, bounding_box[0].second, bounding_box[1].first, bounding_box[1].second };
        memset(&vc.diagram, 0, sizeof(jcv_diagram));

        /*auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Elapsed time -- preprocess: " << elapsed.count() << " ms" << std::endl;*/

        //start = std::chrono::high_resolution_clock::now();

        Alpha_shape_2 A(points.begin(), points.end(),
            FT(60), // 11 - 50+
            Alpha_shape_2::GENERAL);

        std::vector<CGAL_Segment> segs;
        alpha_edges(A, std::back_inserter(segs));
        auto alpha_bdry = segments_to_path(segs);
      
        /*end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Elapsed time -- alphashape: " << elapsed.count() << " ms" << std::endl;*/

        //start = std::chrono::high_resolution_clock::now();

        jcv_diagram_generate(vc.num_points, vc.points, &vc.rect, 0, &vc.diagram);

        /*end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Elapsed time -- voronoi gen: " << elapsed.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();*/

        //const jcv_site* sites = jcv_diagram_get_sites(&vc.diagram);
        //for (int i = 0; i < vc.diagram.numsites; ++i)
        //{
        //    const jcv_site* site = &sites[i];

        //    srand((unsigned int)site->index); // for generating colors for the triangles

        //    unsigned char color_tri[3];
        //    unsigned char basecolor = 120;
        //    color_tri[0] = basecolor + (unsigned char)(rand() % (235 - basecolor));
        //    color_tri[1] = basecolor + (unsigned char)(rand() % (235 - basecolor));
        //    color_tri[2] = basecolor + (unsigned char)(rand() % (235 - basecolor));

        //    jcv_point s = remap(&site->p, &vc.diagram.min, &vc.diagram.max, &dimensions);

        //    const jcv_graphedge* e = site->edges;
        //    while (e)
        //    {
        //        jcv_point p0 = remap(&e->pos[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
        //        jcv_point p1 = remap(&e->pos[1], &vc.diagram.min, &vc.diagram.max, &dimensions);

        //        draw_triangle(&s, &p0, &p1, image, width, height, 3, color_tri);
        //        e = e->next;
        //    }
        //}

        // If all you need are the edges
        const jcv_edge* edge = jcv_diagram_get_edges(&vc.diagram);
        int edgecount = 1;
        std::vector<std::pair<jcv_point, jcv_point>> edge_lines;
        unsigned char color_red[] = { 255, 0, 0 };
        unsigned char color_blue[] = { 75, 75, 230 };
        unsigned char color_green[] = { 0, 255, 0 };
        unsigned char color_white[] = { 255, 255, 255 };
        unsigned char color_orange[] = { 255, 153, 51 };
        unsigned char color_yellow[] = { 250, 255, 0 };
        unsigned char color_cyan[] = { 0, 255, 255 };
        unsigned char color_purple[] = { 255, 0, 255 };
        unsigned char color_lightpurple[] = { 178, 102, 255 };
        std::vector<unsigned char*> colors = { color_blue, color_orange, color_white, color_yellow, color_cyan, color_purple, color_lightpurple };

        while (edge)
        {
            edge_lines.push_back(make_pair(edge->pos[0], edge->pos[1]));
            edge = jcv_diagram_get_next_edge(edge);
            
            edgecount++;
        }

        /*end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Elapsed time -- get edge: " << elapsed.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();*/

        auto filteredCenterlines = filterCenterlines(edge_lines, alpha_bdry, lanelines);
        //std::cout << "Number of centerlines: " << filteredCenterlines.size() << std::endl;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << " -------- Elapsed time: " << elapsed.count() << " ms" << std::endl;
        totaltime += elapsed.count();
        std::cout << "avg time: " << totaltime / (count + 1) << " ms" << std::endl;

        int cor = 0;
        for (const auto& seg : filteredCenterlines) {
            //std::cout << "ctl size: " << seg.size() << std::endl;
            auto color = colors[cor % 6];
            for (int i = 1; i < seg.size(); i++) {
				jcv_point p0 = remap(&seg[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
				jcv_point p1 = remap(&seg[i], &vc.diagram.min, &vc.diagram.max, &dimensions);

                draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);
			}
			cor++;
        }

        /*jcv_delauney_iter delauney;
        jcv_delauney_begin(&vc.diagram, &delauney);
        jcv_delauney_edge delauney_edge;
        unsigned char color_delauney[] = { 64, 64, 255 };
        while (jcv_delauney_next(&delauney, &delauney_edge))
        {
            jcv_point p0 = remap(&delauney_edge.pos[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p1 = remap(&delauney_edge.pos[1], &vc.diagram.min, &vc.diagram.max, &dimensions);
            draw_line((int)p0.x, (int)p0.y, (int)p1.x, (int)p1.y, image, width, height, 3, color_delauney);
        }*/
        


        /* ----------- draw lane ------------- */

        // draw this frame lane lines
        for (const auto& lane : lanes) {
            jcv_point p0 = remap(&lane[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
            for (const auto& point : lane) {
                jcv_point p = remap(&point, &vc.diagram.min, &vc.diagram.max, &dimensions);
                draw_line((double)p0.x, (double)p0.y, (double)p.x, (double)p.y, image, width, height, 3, color_red);
                p0 = p;
            }
        }

        for (uint32_t i = 1; i < alpha_bdry.size(); i++)
        {
            auto curr = remap(&alpha_bdry[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto prev = remap(&alpha_bdry[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
            draw_line(curr.x, curr.y, prev.x, prev.y, image, width, height, 3, color_green);
        }

        // flip image
        int stride = width * 3;
        uint8_t* row = (uint8_t*)malloc((size_t)stride);
        for (int y = 0; y < height / 2; ++y)
        {
            memcpy(row, &image[y * stride], (size_t)stride);
            memcpy(&image[y * stride], &image[(height - 1 - y) * stride], (size_t)stride);
            memcpy(&image[(height - 1 - y) * stride], row, (size_t)stride);
        }

        char path[512];
        sprintf_s(path, "img/out%d.png", count);
        //sprintf_s(path, "out%d.png", count);

        stbi_write_png(path, width, height, 3, image, stride);
        std::cout << "done " << path << std::endl;



        free(image);
        count++;
        //break;
    }

}

double dotProduct(const jcv_point& p1, const jcv_point& p2) {
    return p1.x * p2.x + p1.y * p2.y;
}

// Function to calculate the magnitude squared of a point
double magnitudeSquared(const jcv_point& p) {
    return p.x * p.x + p.y * p.y;
}

// Function to calculate the perpendicular distance from a point to a line segment
double perpendicularDistance(const jcv_point& start, const jcv_point& end, const jcv_point& point) {
    jcv_point lineVec = { end.x - start.x, end.y - start.y };
    jcv_point pointVec = { point.x - start.x, point.y - start.y };

    double lineLenSq = magnitudeSquared(lineVec);
    double t = dotProduct(pointVec, lineVec) / lineLenSq;

    // Clamp t to the range [0, 1]
    t = std::max(0.0, std::min(1.0, t));

    jcv_point projection = { start.x + t * lineVec.x, start.y + t * lineVec.y };
    jcv_point distVec = { point.x - projection.x, point.y - projection.y };

    return std::sqrt(magnitudeSquared(distVec));
}

void runactualvordist()
{
    std::string filename = "lanes.csv";
    std::map<long long, std::vector<std::vector<PointDB>>> lanedata = readActualCSV(filename);
    std::cout << "Number of timestamps: " << lanedata.size() << std::endl;
    int count = 0;
    double totaltime = 0;
    for (const auto& entry : lanedata) {
        std::vector<std::vector<PointDB>> inlanes = entry.second; // all lanes in one timestemp

        std::vector<std::vector<PointDB>> lanes;
        for (const auto& inlane : inlanes) {
            std::vector<PointDB> line;
            for (uint32_t i = 0; i < inlane.size(); i += 1) {
                line.push_back(inlane[i]);
            }
            lanes.push_back(line);
        }

        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::pair<double, double>> bounding_box = findBoundingBox(lanes);

        auto lanelines = pointDBtoJCV(lanes);

        vorcon vc;
        int width = 512;
        int height = 512;
        size_t imagesize = (size_t)(width * height * 3);
        unsigned char* image = (unsigned char*)malloc(imagesize);
        memset(image, 0, imagesize);

        jcv_point dimensions;
        dimensions.x = (jcv_real)width;
        dimensions.y = (jcv_real)height;

        std::vector<CGAL_Point> points = convert_to_cgal_points(lanes, vc);
        vc.rect = { bounding_box[0].first, bounding_box[0].second, bounding_box[1].first, bounding_box[1].second };
        memset(&vc.diagram, 0, sizeof(jcv_diagram));

        Alpha_shape_2 A(points.begin(), points.end(),
            FT(60), // 11 - 50+
            Alpha_shape_2::GENERAL);

        std::vector<CGAL_Segment> segs;
        alpha_edges(A, std::back_inserter(segs));
        auto alpha_bdry = segments_to_path(segs);

        /*jcv_clipping_polygon clpoly;
        jcv_clipper* clipper = 0;
        clpoly.points = alpha_bdry.data();
        clpoly.num_points = alpha_bdry.size();

        jcv_clipper polyclipper;
        polyclipper.test_fn = jcv_clip_polygon_test_point;
        polyclipper.clip_fn = jcv_clip_polygon_clip_edge;
        polyclipper.fill_fn = jcv_clip_polygon_fill_gaps;
        polyclipper.ctx = &clpoly;

        clipper = &polyclipper;*/

        jcv_diagram_generate(vc.num_points, vc.points, &vc.rect, 0, &vc.diagram);        

        // If all you need are the edges
        const jcv_edge* edge = jcv_diagram_get_edges(&vc.diagram);
        int edgecount = 1;
        std::vector<std::pair<jcv_point, jcv_point>> edge_lines;
        unsigned char color_red[] = { 255, 0, 0 };
        unsigned char color_blue[] = { 75, 75, 230 };
        unsigned char color_green[] = { 0, 255, 0 };
        unsigned char color_white[] = { 255, 255, 255 };
        unsigned char color_orange[] = { 255, 153, 51 };
        unsigned char color_yellow[] = { 250, 255, 0 };
        unsigned char color_cyan[] = { 0, 255, 255 };
        unsigned char color_purple[] = { 255, 0, 255 };
        unsigned char color_lightpurple[] = { 178, 102, 255 };
        std::vector<unsigned char*> colors = { color_blue, color_orange, color_white, color_yellow, color_cyan, color_purple, color_lightpurple };

        while (edge)
        {
            auto site = edge->sites;
            auto s0 = site[0];
            auto s1 = site[1];
            if (s0 == NULL || s1 == NULL) {
                edge = jcv_diagram_get_next_edge(edge);
                continue;
            }
            auto p0 = s0->p;
            auto p1 = s1->p;
            auto dist0 = perpendicularDistance(edge->pos[0], edge->pos[1], p0);
            auto dist1 = perpendicularDistance(edge->pos[0], edge->pos[1], p1);
            if (dist0 > 1.25 && dist1 > 1.25) {
                edge_lines.push_back(make_pair(edge->pos[0], edge->pos[1]));
            }
            
            edge = jcv_diagram_get_next_edge(edge);

            edgecount++;
        }

        auto filteredCenterlines = filterCenterlinesByBoundary(edge_lines, alpha_bdry);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << " -------- Elapsed time: " << elapsed.count() << " ms" << std::endl;
        totaltime += elapsed.count();
        std::cout << "avg time: " << totaltime / (count + 1) << " ms" << std::endl;

        int cor = 0;
        /*for (const auto& seg : filteredCenterlines) {
            auto color = colors[cor % 6];
            for (int i = 1; i < seg.size(); i++) {
                jcv_point p0 = remap(&seg[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
                jcv_point p1 = remap(&seg[i], &vc.diagram.min, &vc.diagram.max, &dimensions);

                draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);
            }
            cor++;
        }*/
        for (const auto& seg : edge_lines) {
            auto color = colors[cor % 6];
            jcv_point p0 = remap(&seg.first, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p1 = remap(&seg.second, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);
            cor++;
        }


        /* ----------- draw lane ------------- */

        // draw this frame lane lines
        for (const auto& lane : lanes) {
            jcv_point p0 = remap(&lane[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
            for (const auto& point : lane) {
                jcv_point p = remap(&point, &vc.diagram.min, &vc.diagram.max, &dimensions);
                draw_line((double)p0.x, (double)p0.y, (double)p.x, (double)p.y, image, width, height, 3, color_red);
                p0 = p;
            }
        }

        for (uint32_t i = 1; i < alpha_bdry.size(); i++)
        {
            auto curr = remap(&alpha_bdry[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto prev = remap(&alpha_bdry[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
            draw_line(curr.x, curr.y, prev.x, prev.y, image, width, height, 3, color_green);
        }

        // flip image
        int stride = width * 3;
        uint8_t* row = (uint8_t*)malloc((size_t)stride);
        for (int y = 0; y < height / 2; ++y)
        {
            memcpy(row, &image[y * stride], (size_t)stride);
            memcpy(&image[y * stride], &image[(height - 1 - y) * stride], (size_t)stride);
            memcpy(&image[(height - 1 - y) * stride], row, (size_t)stride);
        }

        char path[512];
        sprintf_s(path, "img/out%d.png", count);
        //sprintf_s(path, "out%d.png", count);

        stbi_write_png(path, width, height, 3, image, stride);
        std::cout << "done " << path << std::endl;



        free(image);
        count++;
        break;
    }

}

void testing() {
    geos::io::WKTReader reader;
    geos::io::WKTWriter writer;
    GeometryFactory::Ptr factory = GeometryFactory::create();
    std::unique_ptr<Geometry> boundary(reader.read("POLYGON ((30 10, 40 40, 20 40, 10 20, 30 10))")); // Simple square boundary

    int gridSize = 10; // Define the size of the grid
    std::vector<std::unique_ptr<Geometry>> gridCells = createGridCells(boundary.get(), gridSize, factory.get());

    // Define lanelines
    std::vector<std::unique_ptr<Geometry>> lanelines;
    lanelines.push_back(std::unique_ptr<Geometry>(reader.read("LINESTRING (35 15, 45 45)")));
    lanelines.push_back(std::unique_ptr<Geometry>(reader.read("LINESTRING (25 25, 35 35)")));

    // Define centerlines
    std::vector<std::unique_ptr<Geometry>> centerlines;
    centerlines.push_back(std::unique_ptr<Geometry>(reader.read("LINESTRING (10 10, 20 20)")));
    centerlines.push_back(std::unique_ptr<Geometry>(reader.read("LINESTRING (20 20, 30 30)")));

    TemplateSTRtree<const Geometry*> gridStrTree;
    TemplateSTRtree<const Geometry*> laneStrTree;
    for (const auto& cell : gridCells) {
        gridStrTree.insert(cell->getEnvelopeInternal(), cell.get());
    }

    // Split and index lanelines
    splitAndIndexLanelines(gridCells, lanelines, laneStrTree);

    // Output the grid cells
    std::cout << "Grid Cells:" << std::endl;
    for (const auto& cell : gridCells) {
        std::cout << writer.write(cell.get()) << std::endl;
    }

    // Example point query
    std::unique_ptr<Point> point(factory->createPoint(Coordinate(15.0, 15.0))); // Example point coordinates
    if (isPointInBoundary(point.get(), gridStrTree)) {
        std::cout << "Point is within the boundary." << std::endl;
    }
    else {
        std::cout << "Point is outside the boundary." << std::endl;
    }

    // Check if centerlines are too close to lanelines
    double minDistance = 5.0; // Minimum distance threshold
    std::cout << "Centerline Proximity Check:" << std::endl;
    for (const auto& centerline : centerlines) {
        if (isCenterlineTooCloseToLanelines(centerline.get(), laneStrTree, minDistance)) {
            std::cout << "A centerline is too close to a laneline." << std::endl;
        }
        else {
            std::cout << "A centerline is at a safe distance from lanelines." << std::endl;
        }
    }
}

int main() {
    //runsampledata();
    //runactualdata();
    //testing();
    //runactualgeos();
    runactualvordist();
	return 0;
}               

