#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <map>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstdlib>
#include <unordered_set>

#define JC_VORONOI_IMPLEMENTATION
#define JC_VORONOI_CLIP_IMPLEMENTATION
#define OUTPIC

#include "DBSCAN.h"
#include "jc_voronoi.h"
#include "jc_voronoi_clip.h"
#include "stb_image_write.h"
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

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> bst_point;
typedef boost::geometry::model::segment<bst_point> bst_segment;
typedef boost::geometry::model::linestring<bst_point> bst_linestring;
typedef boost::geometry::model::polygon<bst_point> bst_polygon;
typedef std::pair<bst_segment, uint32_t> bst_value;

typedef struct VoronoiContext_ {
    uint32_t num_points;
    jcv_point* points;
    jcv_rect rect;
    jcv_diagram diagram;
    jcv_clipper clipper;
    jcv_clipper* clipper_ptr;
} vorcon;

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

void printResults(vector<PointDB>& points, int num_points)
{
    int i = 0;
    printf("Number of points: %u\n"
        " x     y     cluster_id\n"
        "-----------------------------\n"
        , num_points);
    while (i < num_points)
    {
        printf("%5.6lf %5.6lf : %d\n",
            points[i].x,
            points[i].y,
            points[i].clusterID);
        ++i;
    }
}

bool compareByY(const PointDB& a, const PointDB& b) {
    return a.y < b.y;
}

std::vector<PointDB> moving_average(std::vector<PointDB>& data, size_t window_size) {
    size_t n = data.size();
    window_size += 2;

    // Extend data with padding at both ends
    std::vector<PointDB> padded_data = data;
    padded_data.insert(padded_data.begin(), window_size, data.front());
    padded_data.insert(padded_data.end(), window_size, data.back());

    // Calculate moving average including padding
    std::vector<PointDB> smoothed(n);
    for (size_t i = 0; i < n; ++i) {
        double sum_x = 0.0, sum_y = 0.0;
        for (size_t j = 0; j < 2 * window_size + 1; ++j) {
            sum_x += padded_data[i + j].x;
            sum_y += padded_data[i + j].y;
        }
        smoothed[i].x = sum_x / (2 * window_size + 1);
        smoothed[i].y = sum_y / (2 * window_size + 1);
    }

    return smoothed;
}

std::vector<std::vector<PointDB>> apply_moving_average(std::vector<std::vector<PointDB>>& lines) {
    std::vector<std::vector<PointDB>> ret;
    for (auto& line : lines) {
        size_t window_size = std::min(2, (int)(line.size() / 5));
        ret.push_back(moving_average(line, window_size));
    }
    return ret;
}

std::vector<std::pair<int, PointDB>> find_closest_to_corner(const std::vector<std::pair<double, double>>& bounding_box, std::vector<std::vector<PointDB>>& lane_lines) {
    // Extract corners from bounding box
    double x = bounding_box[0].first;
    double y = bounding_box[0].second;
    double x1 = bounding_box[1].first;
    double y1 = bounding_box[1].second;

    PointDB topleft = { x, y1 };
    PointDB topright = { x1, y1 };
    PointDB bottomleft = { x, y };
    PointDB bottomright = { x1, y };

    std::vector<PointDB> corners = { topleft, topright, bottomleft, bottomright };
    std::vector<std::pair<int, PointDB>> closest_points(corners.size(), std::make_pair(-1, PointDB{ 0, 0 }));

    for (size_t i = 0; i < corners.size(); ++i) {
        double min_distance = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < lane_lines.size(); ++j) {
            for (const auto& point : lane_lines[j]) {
                double distance = std::sqrt(std::pow(corners[i].x - point.x, 2) + std::pow(corners[i].y - point.y, 2));
                if (distance < min_distance) {
                    min_distance = distance;
                    closest_points[i] = { static_cast<int>(j), point };
                }
            }
        }
    }

    return closest_points;
}

void merge_lines_on_side(const std::vector<std::pair<double, double>>& bounding_box, std::vector<std::vector<PointDB>>& lane_lines) {
    
    std::vector<std::pair<int, PointDB>> closest_points = find_closest_to_corner(bounding_box, lane_lines);
    
    int leftmost_label = closest_points[0].first;
    int rightmost_label = closest_points[1].first;

    // Merge lanes if necessary
    if (closest_points[0].first != closest_points[2].first) {
        lane_lines[leftmost_label].insert(lane_lines[leftmost_label].end(), lane_lines[closest_points[2].first].begin(), lane_lines[closest_points[2].first].end());
        lane_lines.erase(lane_lines.begin() + closest_points[2].first);
    }
    if (closest_points[1].first != closest_points[3].first) {
        lane_lines[rightmost_label].insert(lane_lines[rightmost_label].end(), lane_lines[closest_points[3].first].begin(), lane_lines[closest_points[3].first].end());
        lane_lines.erase(lane_lines.begin() + closest_points[3].first);
    }
}

std::vector<PointDB> smooth_line(const std::vector<PointDB>& line) {
    // Extract x and y coordinates from the line and swap them
    size_t n = line.size();
    double* x = new double[n];
    double* y = new double[n];
    for (size_t i = 0; i < n; i++) {
        x[i] = line[i].y;  // Swap x and y here for interpolation
        y[i] = line[i].x;
    }

    // Number of points for interpolation
    double min_x = *std::min_element(x, x + n); // Note swapped roles
    double max_x = *std::max_element(x, x + n);
    size_t num_points = static_cast<size_t>((max_x - min_x) * 2);

    // Setup for B-spline interpolation using GSL
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline, x, y, n);

    // Calculate interpolated values
    std::vector<PointDB> smoothed_line(num_points);
    for (size_t i = 0; i < num_points; i++) {
        double xi = min_x + i * (max_x - min_x) / (num_points - 1);
        double yi = gsl_spline_eval(spline, xi, acc);
        smoothed_line[i] = { yi, xi };  // Swap x and y back after interpolation
    }

    // Cleanup
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    delete[] x;
    delete[] y;

    return smoothed_line;
}

// Function to smooth lines using B-spline
std::vector<std::vector<PointDB>> smooth_lines(const std::vector<std::vector<PointDB>>& lines) {
    std::vector<std::vector<PointDB>> smoothed_lines;

    for (const auto& line : lines) {
        /*for (auto& pt : line)
            std::cout << "x: " << pt.x << "; y: " << pt.y << std::endl;*/
        smoothed_lines.push_back(smooth_line(line));
        //std::cout << " ------------------------------------------------------  " << std::endl;
    }

    return smoothed_lines;
}

void convert_to_jcv_points(const std::vector<std::vector<PointDB>>& lane_lines, vorcon& vc) {
	size_t num_points = 0;
	for (const auto& line : lane_lines) {
		num_points += line.size();
	}

	jcv_point* points = new jcv_point[num_points];
	size_t idx = 0;
	for (const auto& line : lane_lines) {
		for (const auto& point : line) {
			points[idx] = { point.x, point.y };
            idx++;
		}
	}
    vc.num_points = num_points;
    vc.points = points;
}

jcv_point* convert_to_jcv_points(const std::vector<PointDB>& lane) {
	jcv_point* points = new jcv_point[lane.size()];
	for (size_t i = 0; i < lane.size(); i++) {
		points[i] = { lane[i].x, lane[i].y };
	}
	return points;
}

void savetofile(std::vector<std::vector<PointDB>> lines) {
    std::ofstream file("data.txt");
    if (!file.is_open()) return;

    for (const auto& line : lines) {
        for (const auto& point : line) {
            file << point.x << "," << point.y << " ";
        }
        file << std::endl; // New line for each line series
    }

    file.close();
}

jcv_point calc_diff(std::vector<std::pair<double, double>> bounding) {
    double width = bounding[1].first - bounding[0].first;
    double height = bounding[1].second - bounding[0].second;
    return { width, height };
}

jcv_point bst_val_to_jcv(const bst_value val, bool zo) {
    auto line = val.first;
    if (zo == true) { // p0
        auto p0 = line.first;
        return { p0.get<0>(), p0.get<1>() };
    }
	else { // p1
		auto p1 = line.second;
		return { p1.get<0>(), p1.get<1>() };
	}
}

std::vector<uint32_t> del_intersect(std::vector<uint32_t> qry_idx, std::vector<uint32_t> to_delete_idx, uint32_t n) {
    const int BITS_PER_BLOCK = 64;
    int num_blocks = (n + BITS_PER_BLOCK - 1) / BITS_PER_BLOCK;
    std::vector<uint64_t> binqry(num_blocks, 0);
    std::vector<uint64_t> bindel(num_blocks, 0);

    for (int index : to_delete_idx) {
        int block = index / BITS_PER_BLOCK;
        int bit = index % BITS_PER_BLOCK;
        bindel[block] |= (1ULL << bit);
    }

    for (int index : qry_idx) {
		int block = index / BITS_PER_BLOCK;
		int bit = index % BITS_PER_BLOCK;
		binqry[block] |= (1ULL << bit);
	}

    std::vector<uint64_t> result(num_blocks, 0);
    for (int i = 0; i < num_blocks; ++i) {
        result[i] = binqry[i] ^ bindel[i];
    }

    std::vector<uint32_t> ret;
    for (int i = 0; i < num_blocks; ++i) {
		for (int j = 0; j < BITS_PER_BLOCK; ++j) {
			if (result[i] & (1ULL << j)) {
				ret.push_back(i * BITS_PER_BLOCK + j);
			}
		}
	}

    return ret;
}

int main()
{
    std::string filename = "data3.csv";
    std::vector<PointDB> lane = readCSV(filename);

    std::vector<std::pair<double, double>> bounding_box = findBoundingBox(lane);
    DBSCAN ds(MINIMUM_POINTS, EPSILON, lane);
    ds.run();

    //printResults(ds.m_points, ds.getTotalPointSize());
    lane = ds.m_points;
    //std::cout << ds.getClusterCount() << std::endl;

    // separate each line
    std::vector<std::vector<PointDB>> lane_lines(ds.getClusterCount());
    for (uint32_t i = 0; i < lane.size(); i++) {
        lane_lines[lane[i].clusterID - 1].push_back(lane[i]);
    }

    merge_lines_on_side(bounding_box, lane_lines);

    // sort by y
    for (uint32_t i = 0; i < lane_lines.size(); i++) {
        std::sort(lane_lines[i].begin(), lane_lines[i].end(), compareByY);
    }
    

    /*for (const auto& line : lane_lines) {
        for (auto& pt : line)
            std::cout << "x: " << pt.x << "; y: " << pt.y << std::endl;
        std::cout << std::endl;
    }*/

    lane_lines = apply_moving_average(lane_lines);

    lane_lines = smooth_lines(lane_lines);

    std::vector<std::pair<int, PointDB>> closest_points = find_closest_to_corner(bounding_box, lane_lines);

    int leftmost_label = closest_points[0].first;
    int rightmost_label = closest_points[1].first;
    
    // get clipping polygon
    std::vector<PointDB> bdry = lane_lines[leftmost_label];
    std::vector<PointDB> rightmostlane = lane_lines[rightmost_label];
    std::reverse(rightmostlane.begin(), rightmostlane.end());
    bdry.insert(bdry.end(), rightmostlane.begin(), rightmostlane.end());
    bdry.push_back(bdry[0]);
    /*std::vector<std::vector<PointDB>> leftmostlane2;
    leftmostlane2.push_back(leftmostlane);
    savetofile(leftmostlane2);*/



    /*  -------------------done preprocessing start generate voronoi---------------------------*/
    vorcon vc;

    convert_to_jcv_points(lane_lines, vc);

    vc.rect = { bounding_box[0].first - 1, bounding_box[0].second - 1, bounding_box[1].first + 1, bounding_box[1].second + 1 };
    memset(&vc.diagram, 0, sizeof(jcv_diagram));

    jcv_diagram_generate(vc.num_points, vc.points, &vc.rect, 0, &vc.diagram);

    /* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
    auto diff = calc_diff(bounding_box);
#ifdef OUTPIC
    int width = 512;
    int height = 512;
#else
    int width = diff.x;
    int height = diff.y;
#endif // 

    const char* outputfile = "example.png";

    size_t imagesize = (size_t)(width * height * 3);
    unsigned char* image = (unsigned char*)malloc(imagesize);
    memset(image, 0, imagesize);

    jcv_point dimensions;
    dimensions.x = (jcv_real)width;
    dimensions.y = (jcv_real)height;

    unsigned char color_pt[] = { 255, 255, 255 };

    unsigned char color_edge[3];
    unsigned char basecolor = 120;

    // If all you need are the edges
    const jcv_edge* edge = jcv_diagram_get_edges(&vc.diagram);
    int edgecount = 1;
    std::vector<std::pair<jcv_point, jcv_point>> edge_lines;
    while (edge)
    {
        color_edge[0] = 50;
        color_edge[1] = 50;
        color_edge[2] = 255;


        jcv_point p0 = remap(&edge->pos[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
        jcv_point p1 = remap(&edge->pos[1], &vc.diagram.min, &vc.diagram.max, &dimensions);
        //draw_line((double)p0.x, (double)p0.y, (double)p1.x, (double)p1.y, image, width, height, 3, color_edge);

        /* --------- get edges --------- */
#ifndef OUTPIC
        std::cout << "p0: " << p0.x << ", " << p0.y << "; p1: " << p1.x << ", " << p1.y << std::endl;
#endif
        edge_lines.push_back(make_pair(p0, p1));
        /* ------- end get edge ---------*/

        edge = jcv_diagram_get_next_edge(edge);
        edgecount++;
    }
    std::cout << "edgecount: " << edgecount << std::endl;

    

    

    bst_polygon bdry_poly;
    std::vector<jcv_point> scaled_bdry;
    for (size_t i = 0; i < bdry.size(); i++) {
        // fit to scale
        auto pt = remap(&bdry[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
        scaled_bdry.push_back(pt);
		boost::geometry::append(bdry_poly.outer(), bst_point(pt.x, pt.y));
	}

    std::vector<bst_linestring> road_lines;
    for (size_t i = 0; i < lane_lines.size(); i++) {
        bst_linestring lstr;
        for (size_t j = 0; j < lane_lines[i].size(); j++) {
            auto pt = remap(&lane_lines[i][j], &vc.diagram.min, &vc.diagram.max, &dimensions);
            lstr.push_back(bst_point(pt.x, pt.y));
        }
        road_lines.push_back(lstr);
    }

    std::vector<bst_value> line_to_check;
    for (size_t i = 0; i < edge_lines.size(); i++) {
        auto pt0 = edge_lines[i].first;
        auto pt1 = edge_lines[i].second;
        bst_segment line = { bst_point(pt0.x, pt0.y), bst_point(pt1.x, pt1.y) };
        bst_value val = make_pair(line, i);
        line_to_check.push_back(val);
	}

    boost::geometry::index::rtree<bst_value, boost::geometry::index::rstar<16>> r_tree(line_to_check.begin(), line_to_check.end());

    std::vector<std::vector<bst_value>> to_delete;
    for (size_t i = 0; i < road_lines.size(); i++) {
        auto line = road_lines[i];
        std::vector<bst_value> to_delete_line;
        /*boost::geometry::model::multi_polygon<bst_polygon> buffer;
        boost::geometry::buffer(line, buffer, (double)1.0);*/

        boost::geometry::index::query(r_tree, boost::geometry::index::intersects(line), std::back_inserter(to_delete_line));
        to_delete.push_back(to_delete_line);
    }

    std::vector<bst_value> qry_result;

    boost::geometry::index::query(r_tree, boost::geometry::index::intersects(bdry_poly), std::back_inserter(qry_result));

    //std::cout << "qry result size: " << qry_result.size() << std::endl;

    std::vector<uint32_t> qry_idx;
    for (uint32_t i = 0; i < qry_result.size(); i++) {
        qry_idx.push_back(qry_result[i].second);
	}

    std::vector<uint32_t> to_delete_idx;
    for (size_t i = 0; i < to_delete.size(); i++) {
		for (size_t j = 0; j < to_delete[i].size(); j++) {
			to_delete_idx.push_back(to_delete[i][j].second);
		}
	}

    /*std::sort(qry_idx.begin(), qry_idx.end());
    std::sort(to_delete_idx.begin(), to_delete_idx.end());
    std::cout << " ------------qry idx---------------" << std::endl;
    for (auto& i : qry_idx) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << " ------------to delete idx---------------" << std::endl;
    for (auto& i : to_delete_idx) {
		std::cout << i << " ";
	}
    std::cout << std::endl;*/
    /*std::unordered_set<uint32_t> to_delete_set(to_delete_idx.begin(), to_delete_idx.end());
    qry_idx.erase(
        std::remove_if(qry_idx.begin(), qry_idx.end(), [&to_delete_set](uint32_t i) { 
            return to_delete_set.find(i) != to_delete_set.end(); 
        }), 
        qry_idx.end());*/


    std::vector<uint32_t> to_keep_idx = del_intersect(qry_idx, to_delete_idx, line_to_check.size());




    
    /*std::cout << "qry idx size: " << qry_idx.size() << std::endl;
    std::cout << " ------------new qry idx---------------" << std::endl;
    for (auto& i : qry_idx) {
        std::cout << i << " ";
    }
    std::cout << std::endl;*/
    
    for (uint32_t i = 0; i < to_keep_idx.size(); i++) {
        jcv_point p0 = bst_val_to_jcv(line_to_check[to_keep_idx[i]], true);
        jcv_point p1 = bst_val_to_jcv(line_to_check[to_keep_idx[i]], false);
        color_edge[0] = 255;
        color_edge[1] = 50;
        color_edge[2] = 50;
        //std::cout << "p0: " << p0.x << ", " << p0.y << "; p1: " << p1.x << ", " << p1.y << std::endl;
        draw_line((double)p0.x, (double)p0.y, (double)p1.x, (double)p1.y, image, width, height, 3, color_edge);
    }



    /* -------------- image ------------------ */
    // Plot the sites ---- lane points
    for (uint32_t i = 0; i < vc.num_points; ++i)
    {
        jcv_point p = remap(&vc.points[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
        plot((double)p.x, (double)p.y, image, width, height, 3, color_pt);
    }

    for (uint32_t i = 1; i < scaled_bdry.size(); i++)
    {
        auto curr = scaled_bdry[i];
        auto prev = scaled_bdry[i - 1];
        draw_line(curr.x, curr.y, prev.x, prev.y, image, width, height, 3, color_pt);
    }

    /* ------------------------------------------------------------ */

    

    /* ------------------------------------------------------------- */

    //free(clippoints);
    /*free(vc.points);
    free(&vc.rect);*/

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
    sprintf_s(path, "%s", outputfile);

    stbi_write_png(path, width, height, 3, image, stride);
    printf("Done.\n");

    free(image);
    jcv_diagram_free(&vc.diagram);

    return 0;

}



