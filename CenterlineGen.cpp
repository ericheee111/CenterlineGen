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
// #define OUTPIC

#include "DBSCAN.h"
#include "jc_voronoi.h"
#include "jc_voronoi_clip.h"
#include "stb_image_write.h"
#include "strtree.hpp"
#include "vor_manager.h"

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
// #include <boost/timer/timer.hpp>
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2 CGAL_Point;
typedef K::Segment_2 CGAL_Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;

typedef struct VoronoiContext_
{
    uint32_t num_points;
    jcv_point *points;
    jcv_rect rect;
    jcv_diagram diagram;
    jcv_clipper clipper;
    jcv_clipper *clipper_ptr;
} vorcon;

std::map<long long, std::vector<std::vector<jcv_point>>> readActualCSV(const std::string &filename)
{
    std::cout << "read file" << std::endl;
    std::ifstream file(filename);
    std::map<long long, std::vector<std::vector<jcv_point>>> data;

    if (!file.is_open())
    {
        std::cerr << "Unable to open file" << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string cell;

        long long timestamp;
        std::vector<std::vector<jcv_point>> lanes;

        try
        {
            // Read timestamp
            if (!std::getline(ss, cell, ','))
                throw std::invalid_argument("Missing timestamp");
            timestamp = std::stoll(cell);
            if (!std::getline(ss, cell, ','))
                throw std::invalid_argument("Missing lane number");
            int laneNumber = std::stoi(cell);

            // Read lane coordinates
            while (std::getline(ss, cell, ','))
            {
                std::vector<jcv_point> lane;

                do
                {
                    jcv_point coordinate;

                    // Read x coordinate
                    coordinate.x = std::stod(cell);

                    // Check if next cell is available for y coordinate
                    if (!std::getline(ss, cell, ','))
                        break;
                    coordinate.y = std::stod(cell);

                    lane.push_back(coordinate);
                } while (std::getline(ss, cell, ','));

                lanes.push_back(lane);
            }
            data[timestamp].insert(data[timestamp].end(), lanes.begin(), lanes.end());
        }
        catch (const std::exception &e)
        {
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
static void plot(int x, int y, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
{
    if (x < 0 || y < 0 || x > (width - 1) || y > (height - 1))
        return;
    int index = y * width * nchannels + x * nchannels;
    for (int i = 0; i < nchannels; ++i)
    {
        image[index + i] = color[i];
    }
}

// http://members.chello.at/~easyfilter/bresenham.html
static void draw_line(int x0, int y0, int x1, int y1, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
{
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2; // error value e_xy

    for (;;)
    { // loop
        plot(x0, y0, image, width, height, nchannels, color);
        if (x0 == x1 && y0 == y1)
            break;
        e2 = 2 * err;
        if (e2 >= dy)
        {
            err += dy;
            x0 += sx;
        } // e_xy+e_x > 0
        if (e2 <= dx)
        {
            err += dx;
            y0 += sy;
        } // e_xy+e_y < 0
    }
}

// http://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
static inline int orient2d(const jcv_point *a, const jcv_point *b, const jcv_point *c)
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
static inline jcv_point remap(const jcv_point *pt, const jcv_point *min, const jcv_point *max, const jcv_point *scale)
{
    jcv_point p;
    p.x = (pt->x - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static inline jcv_point remap(const PointDB *pt, const jcv_point *min, const jcv_point *max, const jcv_point *scale)
{
    jcv_point p;
    p.x = (pt->x - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static inline jcv_point remap(const CGAL_Point *pt, const jcv_point *min, const jcv_point *max, const jcv_point *scale)
{
    jcv_point p;
    p.x = (pt->x() - min->x) / (max->x - min->x) * scale->x;
    p.y = (pt->y() - min->y) / (max->y - min->y) * scale->y;
    return p;
}

static void draw_triangle(const jcv_point *v0, const jcv_point *v1, const jcv_point *v2, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
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
    for (p.y = (jcv_real)minY; p.y <= (jcv_real)maxY; p.y++)
    {
        for (p.x = (jcv_real)minX; p.x <= (jcv_real)maxX; p.x++)
        {
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
std::vector<PointDB> readCSV(const std::string &filename)
{
    std::vector<PointDB> points;
    std::ifstream file(filename);
    std::string line, value;

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return points;
    }

    while (getline(file, line))
    {
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

double dotProduct(const jcv_point &p1, const jcv_point &p2)
{
    return p1.x * p2.x + p1.y * p2.y;
}

double crossProduct(const jcv_point &p1, const jcv_point &p2)
{
    return p1.x * p2.y - p1.y * p2.x;
}

// Function to calculate the magnitude squared of a point
double magnitudeSquared(const jcv_point &p)
{
    return p.x * p.x + p.y * p.y;
}

// Function to calculate the perpendicular distance from a point to a line segment
double perpendicularDistance(const jcv_point &start, const jcv_point &end, const jcv_point &point)
{
    jcv_point lineVec = {end.x - start.x, end.y - start.y};
    jcv_point pointVec = {point.x - start.x, point.y - start.y};

    double lineLenSq = magnitudeSquared(lineVec);
    double t = dotProduct(pointVec, lineVec) / lineLenSq;

    // Clamp t to the range [0, 1]
    t = std::max(0.0, std::min(1.0, t));

    jcv_point projection = {start.x + t * lineVec.x, start.y + t * lineVec.y};
    jcv_point distVec = {point.x - projection.x, point.y - projection.y};

    return std::sqrt(magnitudeSquared(distVec));
}

int calcShadowDist(const jcv_point &pt, const jcv_point &line_start_pt, const jcv_point &line_end_pt,
                   jcv_point *project_pt, double *coef, double *dist)
{
    if (project_pt == nullptr)
    {
        return 1;
    }

    jcv_point direction;
    direction.x = line_end_pt.x - line_start_pt.x;
    direction.y = line_end_pt.y - line_start_pt.y;
    auto line_len = std::sqrt(magnitudeSquared(direction));
    if (line_len == 0)
    {
        *project_pt = line_start_pt;
        *coef = 0;
        std::cout << "line_len == 0, start == end" << std::endl;
        return 1;
    }

    auto unit_direction = direction;
    unit_direction.x /= line_len;
    unit_direction.y /= line_len;
    jcv_point v1;
    v1.x = pt.x - line_start_pt.x;
    v1.y = pt.y - line_start_pt.y;
    auto project_len = dotProduct(v1, unit_direction);
    auto _coef = project_len / line_len;
    project_pt->x = (line_end_pt.x - line_start_pt.x) * _coef + line_start_pt.x;
    project_pt->y = (line_end_pt.y - line_start_pt.y) * _coef + line_start_pt.y;
    if (coef != nullptr && dist != nullptr)
    {
        *coef = _coef;
        *dist = crossProduct(v1, direction) / line_len;
    }

    return 0;
}

double calcDist(const jcv_point &pt, const jcv_point &line_start_pt, const jcv_point &line_end_pt)
{
    jcv_point project_pt;
    double coef, dist;
    calcShadowDist(pt, line_start_pt, line_end_pt, &project_pt, &coef, &dist);
    if (coef < 0)
    {
        // std::cout << "coef < 0" << std::endl;
        dist = std::hypot(pt.x - line_start_pt.x, pt.y - line_start_pt.y);
    }
    else if (coef > 1)
    {
        // std::cout << "coef > 1" << std::endl;
        dist = std::hypot(pt.x - line_end_pt.x, pt.y - line_end_pt.y);
    }

    return dist;
}

std::vector<std::pair<double, double>> findBoundingBox(const std::vector<PointDB> &points)
{
    // Initialize min and max values
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::min();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::min();

    // Iterate through each point
    for (const auto &point : points)
    {
        double x = point.x;
        double y = point.y;

        // Update min and max values
        if (x < minX)
            minX = x;
        if (x > maxX)
            maxX = x;
        if (y < minY)
            minY = y;
        if (y > maxY)
            maxY = y;
    }

    // Return the bounding box as a vector of pairs
    return {{minX, minY}, {maxX, maxY}};
}

std::vector<std::pair<double, double>> findBoundingBox(const std::vector<std::vector<jcv_point>> &lanes)
{
    // Initialize min and max values
    double minX = DBL_MAX;
    double maxX = -DBL_MAX;
    double minY = DBL_MAX;
    double maxY = -DBL_MAX;

    // Iterate through each point
    for (const auto &points : lanes)
    {
        for (const auto &point : points)
        {
            double x = point.x;
            double y = point.y;

            // Update min and max values
            if (x < minX)
                minX = x;
            if (x > maxX)
                maxX = x;
            if (y < minY)
                minY = y;
            if (y > maxY)
                maxY = y;
        }
    }

    // Return the bounding box as a vector of pairs
    return {{minX, minY}, {maxX, maxY}};
}

std::vector<CGAL_Point> convert_to_cgal_points(const std::vector<std::vector<jcv_point>> &lane_lines, vorcon &vc)
{
    size_t num_points = 0;
    for (const auto &line : lane_lines)
    {
        num_points += line.size();
    }

    std::vector<CGAL_Point> cgalpoints;

    jcv_point *points = new jcv_point[num_points];
    size_t idx = 0;
    for (const auto &line : lane_lines)
    {
        for (const auto &point : line)
        {
            points[idx] = {point.x, point.y};
            cgalpoints.push_back(CGAL_Point(point.x, point.y));
            idx++;
        }
    }
    vc.num_points = num_points;
    vc.points = points;
    return cgalpoints;
}

template <class OutputIterator>
void alpha_edges(const Alpha_shape_2 &A, OutputIterator out)
{
    for (Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(); it != A.alpha_shape_edges_end(); ++it)
    {
        *out++ = A.segment(*it);
    }
}

std::vector<jcv_point> segments_to_path(const std::vector<CGAL_Segment> &segments)
{
    std::map<jcv_point, std::vector<jcv_point>> adjacency_list;
    std::set<std::pair<jcv_point, jcv_point>> visited_segments;

    // Construct adjacency list
    for (const auto &segment : segments)
    {
        jcv_point source = {segment.source().x(), segment.source().y()};
        jcv_point target = {segment.target().x(), segment.target().y()};
        adjacency_list[source].push_back(target);
        adjacency_list[target].push_back(source);
    }

    // Start from the first segment's source
    std::vector<jcv_point> path;
    jcv_point start = {segments.front().source().x(), segments.front().source().y()};
    jcv_point current = start;
    path.push_back(current);

    // To handle branches and backtracking
    std::vector<std::pair<jcv_point, size_t>> stack;

    while (true)
    {
        bool found_next = false;

        for (size_t i = 0; i < adjacency_list[current].size(); ++i)
        {
            jcv_point next = adjacency_list[current][i];
            if (visited_segments.find({current, next}) == visited_segments.end())
            {
                stack.push_back({current, i});
                visited_segments.insert({current, next});
                visited_segments.insert({next, current});
                current = next;
                path.push_back(current);
                found_next = true;
                break;
            }
        }

        if (!found_next)
        {
            if (current == start)
            {
                break; // We have returned to the start, so the polygon is complete
            }
            if (stack.empty())
            {
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

    if (!path.empty())
    {
        path.push_back(start); // Close the path only if a valid path is found
    }

    return path;
}

jcv_point calc_midpoint(const jcv_point &p1, const jcv_point &p2)
{
    return {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
}

struct ComparePoints
{
    bool operator()(const jcv_point &p1, const jcv_point &p2)
    {
        return p1.x > p2.x; // Smaller x first
    }
};

std::vector<jcv_point> findBoundary(const std::vector<std::vector<jcv_point>> lanes)
{
    std::priority_queue<jcv_point, std::vector<jcv_point>, ComparePoints> pq;
    for (const auto &lane : lanes)
    {
        for (const auto &point : lane)
        {
            pq.push(point);
        }
    }

    std::vector<jcv_point> left;
    std::vector<jcv_point> right;
    while (!pq.empty())
    {
        auto pt = pq.top();
        double start = pt.x;
        double end = start + 2;
        jcv_point miny = pt;
        jcv_point maxy = pt;
        pq.pop();
        while (!pq.empty() && pq.top().x < end)
        {
            auto next = pq.top();
            pq.pop();
            if (next.y < miny.y)
            {
                miny = next;
            }
            else if (next.y > maxy.y)
            {
                maxy = next;
            }
        }
        left.push_back(miny);
        right.push_back(maxy);
    }
    std::reverse(right.begin(), right.end());
    left.insert(left.end(), right.begin(), right.end());
    left.push_back(left.front());
    return left;
}

void test() {
    std::cout << "run" << std::endl;
    std::vector<std::vector<jcv_point>> lanelines = { {{1.5484909e5, -0.21024238e5},
                                                          {1.5413352e5, -0.21024236e5},
                                                          {1.5378445e5, -0.21024236e5},
                                                          {1.5343539e5, -0.21024235e5},
                                                          {1.5308806e5, -0.19333011e5}} }; // all lanes in one timestemp
    vorcon vc;
    int width = 256;
    int height = 256;
    size_t imagesize = (size_t)(width * height * 3);
    unsigned char* image = (unsigned char*)malloc(imagesize);
    memset(image, 0, imagesize);

    jcv_point dimensions;
    dimensions.x = (jcv_real)width;
    dimensions.y = (jcv_real)height;

    std::vector<CGAL_Point> points = convert_to_cgal_points(lanelines, vc); // convert + store points in vc
    memset(&vc.diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(vc.num_points, (const jcv_point*)vc.points, 0, 0, &vc.diagram);
    const jcv_edge* edge = jcv_diagram_get_edges(&vc.diagram);
    std::vector<std::pair<jcv_point, jcv_point>> edge_lines;
    while (edge)
    {
        
        edge_lines.push_back(make_pair(edge->pos[0], edge->pos[1]));
        edge = jcv_diagram_get_next_edge(edge);
    }

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
    int cor = 0;
    std::cout << "edge line size: " << edge_lines.size() << std::endl;
    for (auto& e : edge_lines) {
         /*auto color = colors[cor % 6];*/
        unsigned char color[] = { rand()%255, rand() % 255, rand() % 255 };

         jcv_point p0 = remap(&e.first, &vc.diagram.min, &vc.diagram.max, &dimensions);
         jcv_point p1 = remap(&e.second, &vc.diagram.min, &vc.diagram.max, &dimensions);
         std::cout << "edge: " << p0.x << " " << p0.y << "; " << p1.x << " " << p1.y << std::endl;

         draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);


         cor++;
    }

    for (auto& pt : lanelines[0]) {
        jcv_point p = remap(&pt, &vc.diagram.min, &vc.diagram.max, &dimensions);
        plot(p.x+1, p.y+1, image, width, height, 3, color_white);
        plot(p.x-1, p.y-1, image, width, height, 3, color_white);
        plot(p.x+1, p.y-1, image, width, height, 3, color_white);
        plot(p.x-1, p.y+1, image, width, height, 3, color_white);
        plot(p.x+1, p.y, image, width, height, 3, color_white);
        plot(p.x-1, p.y, image, width, height, 3, color_white);
        plot(p.x, p.y+1, image, width, height, 3, color_white);
        plot(p.x, p.y-1, image, width, height, 3, color_white);
        plot(p.x, p.y, image, width, height, 3, color_white);
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
    sprintf_s(path, "testimg.png");
    // sprintf_s(path, "out%d.png", count);

    stbi_write_png(path, width, height, 3, image, stride);
    std::cout << "done " << path << std::endl;

    jcv_diagram_free(&vc.diagram);
    free(image);

}

void runactualvordist()
{
    std::string filename = "lanes.csv";
    std::map<long long, std::vector<std::vector<jcv_point>>> lanedata = readActualCSV(filename);
    std::cout << "Number of timestamps: " << lanedata.size() << std::endl;
    int count = 0;
    double totaltime = 0;
    for (const auto &entry : lanedata)
    {
        /*if (count != 95) {
            count++;
            continue;
        }*/
          std::vector<std::vector<jcv_point>> lanelines = entry.second; // all lanes in one timestemp
        /*std::vector<std::vector<jcv_point>> lanelines = {{{1.5484909e5, -0.21024238e5},
                                                          {1.5413352e5, -0.21024236e5},
                                                          {1.5378445e5, -0.21024236e5},
                                                          {1.5343539e5, -0.21024235e5},
                                                          {1.5308806e5, -0.19333011e5}}};*/ // all lanes in one timestemp

        /*std::vector<std::vector<jcv_point>> lanelines;
        for (const auto& inlane : inlanes) {
            std::vector<jcv_point> line;
            for (uint32_t i = 0; i < inlane.size(); i  += 1) {
                line.push_back(inlane[i]);
            }
            lanelines.push_back(line);
        }*/

        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::pair<double, double>> bounding_box = findBoundingBox(lanelines);
        // std::cout << "bounding box: " << bounding_box[0].first << " " << bounding_box[0].second << " " << bounding_box[1].first << " " << bounding_box[1].second << std::endl;

        vorcon vc;
        int width = 512;
        int height = 512;
        size_t imagesize = (size_t)(width * height * 3);
        unsigned char *image = (unsigned char *)malloc(imagesize);
        memset(image, 0, imagesize);

        jcv_point dimensions;
        dimensions.x = (jcv_real)width;
        dimensions.y = (jcv_real)height;

        std::vector<CGAL_Point> points = convert_to_cgal_points(lanelines, vc); // convert + store points in vc
        // vc.rect = { bounding_box[0].first, bounding_box[0].second, bounding_box[1].first, bounding_box[1].second };
        memset(&vc.diagram, 0, sizeof(jcv_diagram));

        // Alpha_shape_2 A(points.begin(), points.end(),
        //     FT(60), // 11 - 50+
        //     Alpha_shape_2::GENERAL);

        // std::vector<CGAL_Segment> segs;
        // alpha_edges(A, std::back_inserter(segs));
        // auto alpha_bdry = segments_to_path(segs);

        // auto alpha_bdry = findBoundary(lanelines);
        std::vector<jcv_point> alpha_bdry;
        //std::cout << "before voronoi" << std::endl;
        jcv_diagram_generate(vc.num_points, (const jcv_point *)vc.points, 0, 0, &vc.diagram);
        //std::cout << "after voronoi" << std::endl;
        //  If all you need are the edges
        const jcv_edge *edge = jcv_diagram_get_edges(&vc.diagram);
        //int edgecount = 1;
        std::vector<std::pair<jcv_point, jcv_point>> edge_lines;

        vor_manager vm;

        while (edge)
        {
            auto site = edge->sites;
            auto s0 = site[0];
            auto s1 = site[1];
            if (s0 == NULL || s1 == NULL)
            {
                edge = jcv_diagram_get_next_edge(edge);
                continue;
            }
            auto p0 = s0->p;
            auto p1 = s1->p;
            auto dist0 = fabs(calcDist(p0, edge->pos[0], edge->pos[1]));
            auto sitedist = std::hypot(p0.x - p1.x, p0.y - p1.y);
            if (dist0 > 1.3 && dist0 < 6 && sitedist >= 2.5)
            {
                edge_lines.push_back(make_pair(edge->pos[0], edge->pos[1]));
                vm.addEdge(edge, s0, s1);
            }

            edge = jcv_diagram_get_next_edge(edge);

            //edgecount++;
        }
        //std::cout << "edge cnt: " << edge_lines.size() << std::endl;

        /*for (int i = 0; i < edge_lines.size(); i++) {
            auto spair = vm.getEdgeSites(edge_lines[i].first, edge_lines[i].second);
            auto s0 = spair.first;
            auto s1 = spair.second;
            std::cout << "edge " << edge_lines[i].first.x << " " << edge_lines[i].first.y
                << ", " << edge_lines[i].second.x << " " << edge_lines[i].second.y << "; s0: "
                << s0->p.x << " " << s0->p.y << ", s1: " << s1->p.x << " " << s1->p.y << std::endl;
        }*/

        auto filteredCenterlines = filterCenterlinesByBoundary(edge_lines, alpha_bdry);

        // std::cout << "centerlines size: " << filteredCenterlines.size() << std::endl;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << " -------- Elapsed time: " << elapsed.count() << " ms" << std::endl;
        totaltime += elapsed.count();
        std::cout << "avg time: " << totaltime / (count + 1) << " ms" << std::endl;

        unsigned char color_red[] = {255, 0, 0};
        unsigned char color_blue[] = {75, 75, 230};
        unsigned char color_green[] = {0, 255, 0};
        unsigned char color_white[] = {255, 255, 255};
        unsigned char color_orange[] = {255, 153, 51};
        unsigned char color_yellow[] = {250, 255, 0};
        unsigned char color_cyan[] = {0, 255, 255};
        unsigned char color_purple[] = {255, 0, 255};
        unsigned char color_lightpurple[] = {178, 102, 255};
        std::vector<unsigned char *> colors = {color_blue, color_orange, color_white, color_yellow, color_cyan, color_purple, color_lightpurple};

        /*for (const auto& seg : edge_lines) {
            auto color = colors[cor % 6];
            jcv_point p0 = remap(&seg.first, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p1 = remap(&seg.second, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);
            cor++;
        }*/

        /*for (const auto& line : filteredCenterlines) {
            //auto color = colors[cor % 6];
            //std::cout << "line size: " << line.size() << std::endl;
            if (line.size() <= 2) {

                continue;
            }
            //std::cout << "-------------------------------------------------------------------" << std::endl;
            auto end00 = line[0];
            auto end01 = line[1];
            auto end10 = line[line.size() - 2];
            auto end11 = line[line.size() - 1];
            auto spair = vm.getEdgeSites(end00, end01);
            auto s00 = spair.first;
            auto s01 = spair.second;
            std::cout << "edge1: " << end00.x << " " << end00.y
                << ", " << end01.x << " " << end01.y << "; s0: "
                << s00->p.x << " " << s00->p.y << ", s1: " << s01->p.x << " " << s01->p.y << std::endl;

            jcv_point p0 = remap(&s00->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p1 = remap(&s01->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto mid0 = calc_midpoint(end00, end01);
            jcv_point pm0 = remap(&mid0, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(pm0.x, pm0.y, p1.x, p1.y, image, width, height, 3, color_white);
            draw_line(pm0.x, pm0.y, p0.x, p0.y, image, width, height, 3, color_white);

            auto spair1 = vm.getEdgeSites(end10, end11);
            auto s10 = spair1.first;
            auto s11 = spair1.second;
            std::cout << "edge2: " << end10.x << " " << end10.y
                << ", " << end11.x << " " << end11.y << "; s0: "
                << s10->p.x << " " << s10->p.y << ", s1: " << s11->p.x << " " << s11->p.y << std::endl;

            jcv_point p10 = remap(&s10->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p11 = remap(&s11->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto mid1 = calc_midpoint(end10, end11);
            jcv_point pm1 = remap(&mid1, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(pm1.x, pm1.y, p11.x, p11.y, image, width, height, 3, color_white);
            draw_line(pm1.x, pm1.y, p10.x, p10.y, image, width, height, 3, color_white);
        }
        std::cout << "ctl size: " << filteredCenterlines.size() << std::endl;*/
        for (auto it = filteredCenterlines.begin(); it != filteredCenterlines.end();)
        {
            // auto color = colors[cor % 6];
            // std::cout << "line size: " << line.size() << std::endl;
            auto line = filteredCenterlines.at(it - filteredCenterlines.begin());
            if (line.size() <= 2)
            {
                it = filteredCenterlines.erase(it);
                continue;
            }
            ++it;
            // std::cout << "-------------------------------------------------------------------" << std::endl;
            auto end00 = line[0];
            auto end01 = line[1];
            auto end10 = line[line.size() - 2];
            auto end11 = line[line.size() - 1];
            auto spair = vm.getEdgeSites(end00, end01);
            auto s00 = spair.first;
            auto s01 = spair.second;
            /*std::cout << "edge1: " << end00.x << " " << end00.y
                << ", " << end01.x << " " << end01.y << "; s0: "
                << s00->p.x << " " << s00->p.y << ", s1: " << s01->p.x << " " << s01->p.y << std::endl;*/

            jcv_point p0 = remap(&s00->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p1 = remap(&s01->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto mid0 = calc_midpoint(end00, end01);
            jcv_point pm0 = remap(&mid0, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(pm0.x, pm0.y, p1.x, p1.y, image, width, height, 3, color_white);
            draw_line(pm0.x, pm0.y, p0.x, p0.y, image, width, height, 3, color_white);

            auto spair1 = vm.getEdgeSites(end10, end11);
            auto s10 = spair1.first;
            auto s11 = spair1.second;
            /*std::cout << "edge2: " << end10.x << " " << end10.y
                << ", " << end11.x << " " << end11.y << "; s0: "
                << s10->p.x << " " << s10->p.y << ", s1: " << s11->p.x << " " << s11->p.y << std::endl;*/

            jcv_point p10 = remap(&s10->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            jcv_point p11 = remap(&s11->p, &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto mid1 = calc_midpoint(end10, end11);
            jcv_point pm1 = remap(&mid1, &vc.diagram.min, &vc.diagram.max, &dimensions);

            draw_line(pm1.x, pm1.y, p11.x, p11.y, image, width, height, 3, color_white);
            draw_line(pm1.x, pm1.y, p10.x, p10.y, image, width, height, 3, color_white);
        }

        int cor = 0;
        for (const auto &seg : filteredCenterlines)
        {
            auto color = colors[cor % 6];
            for (int i = 1; i < seg.size(); i++)
            {
                jcv_point p0 = remap(&seg[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
                jcv_point p1 = remap(&seg[i], &vc.diagram.min, &vc.diagram.max, &dimensions);

                draw_line(p0.x, p0.y, p1.x, p1.y, image, width, height, 3, color);
            }
            cor++;
        }
        // std::cout << "aft ctl size: " << filteredCenterlines.size() << std::endl;

        /* ----------- draw lane ------------- */

        // draw this frame lane lines
        for (const auto &lane : lanelines)
        {
            jcv_point p0 = remap(&lane[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
            for (const auto &point : lane)
            {
                jcv_point p = remap(&point, &vc.diagram.min, &vc.diagram.max, &dimensions);
                draw_line((double)p0.x, (double)p0.y, (double)p.x, (double)p.y, image, width, height, 3, color_red);
                p0 = p;
            }
        }

        /*for (uint32_t i = 1; i < alpha_bdry.size(); i++)
        {
            auto curr = remap(&alpha_bdry[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
            auto prev = remap(&alpha_bdry[i - 1], &vc.diagram.min, &vc.diagram.max, &dimensions);
            draw_line(curr.x, curr.y, prev.x, prev.y, image, width, height, 3, color_green);
        }*/

        // flip image
        int stride = width * 3;
        uint8_t *row = (uint8_t *)malloc((size_t)stride);
        for (int y = 0; y < height / 2; ++y)
        {
            memcpy(row, &image[y * stride], (size_t)stride);
            memcpy(&image[y * stride], &image[(height - 1 - y) * stride], (size_t)stride);
            memcpy(&image[(height - 1 - y) * stride], row, (size_t)stride);
        }

        char path[512];
        sprintf_s(path, "img/out%d.png", count);
        // sprintf_s(path, "out%d.png", count);

        stbi_write_png(path, width, height, 3, image, stride);
        std::cout << "done " << path << std::endl;

        jcv_diagram_free(&vc.diagram);
        free(image);
        count++;
         //break;
    }
}

void testing()
{
    jcv_point p0 = {0, 0};
    jcv_point p1 = {2, 3};
    jcv_point p2 = {3, 4};
    auto x = calcDist(p0, p1, p2);
    std::cout << "pt: " << p0.x << " " << p0.y << ", line st~ed: " << p1.x << " " << p1.y << ", " << p2.x << " " << p2.y << ", dist: " << x << std::endl;
}

int main()
{
    // testing();
    runactualvordist();
    //test();
    return 0;
}
