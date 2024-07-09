#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <map>
#include "DBSCAN.h"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <numeric>
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"
#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "jc_voronoi_clip.h"
#include <cstdlib>
#include "stb_image_write.h"

extern "C" int wrap_stbi_write_png(char const* filename, int w, int h, int comp, const void* data, int stride_in_bytes);

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

static void relax_points(const jcv_diagram* diagram, jcv_point* points)
{
    const jcv_site* sites = jcv_diagram_get_sites(diagram);
    for (int i = 0; i < diagram->numsites; ++i)
    {
        const jcv_site* site = &sites[i];
        jcv_point sum = site->p;
        int count = 1;

        const jcv_graphedge* edge = site->edges;

        while (edge)
        {
            sum.x += edge->pos[0].x;
            sum.y += edge->pos[0].y;
            ++count;
            edge = edge->next;
        }

        points[site->index].x = sum.x / (jcv_real)count;
        points[site->index].y = sum.y / (jcv_real)count;
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

/* ----------------------------------------------------------------------- */



std::vector<Point> readCSV(const std::string& filename) {
    std::vector<Point> points;
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

std::vector<std::pair<double, double>> findBoundingBox(const std::vector<Point>& points) {
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

void printResults(vector<Point>& points, int num_points)
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

bool compareByY(const Point& a, const Point& b) {
    return a.y < b.y;
}

std::vector<Point> moving_average(std::vector<Point>& data, size_t window_size) {
    size_t n = data.size();
    window_size += 2;

    // Extend data with padding at both ends
    std::vector<Point> padded_data = data;
    padded_data.insert(padded_data.begin(), window_size, data.front());
    padded_data.insert(padded_data.end(), window_size, data.back());

    // Calculate moving average including padding
    std::vector<Point> smoothed(n);
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

std::vector<std::vector<Point>> apply_moving_average(std::vector<std::vector<Point>>& lines) {
    std::vector<std::vector<Point>> ret;
    for (auto& line : lines) {
        size_t window_size = std::min(2, (int)(line.size() / 5));
        ret.push_back(moving_average(line, window_size));
    }
    return ret;
}

std::vector<std::pair<int, Point>> find_closest_to_corner(const std::vector<std::pair<double, double>>& bounding_box, std::vector<std::vector<Point>>& lane_lines) {
    // Extract corners from bounding box
    double x = bounding_box[0].first;
    double y = bounding_box[0].second;
    double x1 = bounding_box[1].first;
    double y1 = bounding_box[1].second;

    Point topleft = { x, y1 };
    Point topright = { x1, y1 };
    Point bottomleft = { x, y };
    Point bottomright = { x1, y };

    std::vector<Point> corners = { topleft, topright, bottomleft, bottomright };
    std::vector<std::pair<int, Point>> closest_points(corners.size(), std::make_pair(-1, Point{ 0, 0 }));

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

void merge_lines_on_side(const std::vector<std::pair<double, double>>& bounding_box, std::vector<std::vector<Point>>& lane_lines) {
    
    std::vector<std::pair<int, Point>> closest_points = find_closest_to_corner(bounding_box, lane_lines);
    
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

std::vector<Point> smooth_line(const std::vector<Point>& line) {
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
    std::vector<Point> smoothed_line(num_points);
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
std::vector<std::vector<Point>> smooth_lines(const std::vector<std::vector<Point>>& lines) {
    std::vector<std::vector<Point>> smoothed_lines;

    for (const auto& line : lines) {
        /*for (auto& pt : line)
            std::cout << "x: " << pt.x << "; y: " << pt.y << std::endl;*/
        smoothed_lines.push_back(smooth_line(line));
        //std::cout << " ------------------------------------------------------  " << std::endl;
    }

    return smoothed_lines;
}

void convert_to_jcv_points(const std::vector<std::vector<Point>>& lane_lines, vorcon& vc) {
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

jcv_point* convert_to_jcv_points(const std::vector<Point>& lane) {
	jcv_point* points = new jcv_point[lane.size()];
	for (size_t i = 0; i < lane.size(); i++) {
		points[i] = { lane[i].x, lane[i].y };
	}
	return points;
}

void savetofile(std::vector<std::vector<Point>> lines) {
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

int main()
{
    std::string filename = "data3.csv";
    std::vector<Point> lane = readCSV(filename);

    std::vector<std::pair<double, double>> bounding_box = findBoundingBox(lane);
    DBSCAN ds(MINIMUM_POINTS, EPSILON, lane);
    ds.run();

    //printResults(ds.m_points, ds.getTotalPointSize());
    lane = ds.m_points;
    //std::cout << ds.getClusterCount() << std::endl;

    // separate each line
    std::vector<std::vector<Point>> lane_lines(ds.getClusterCount());
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

    std::vector<std::pair<int, Point>> closest_points = find_closest_to_corner(bounding_box, lane_lines);

    int leftmost_label = closest_points[0].first;
    int rightmost_label = closest_points[1].first;
    
    // get clipping polygon
    std::vector<Point> bdry = lane_lines[leftmost_label];
    std::vector<Point> rightmostlane = lane_lines[rightmost_label];
    std::reverse(rightmostlane.begin(), rightmostlane.end());
    bdry.insert(bdry.end(), rightmostlane.begin(), rightmostlane.end());
    /*std::vector<std::vector<Point>> leftmostlane2;
    leftmostlane2.push_back(leftmostlane);
    savetofile(leftmostlane2);*/



    /*  -------------------done preprocessing start generate voronoi---------------------------*/
    vorcon vc;

    jcv_clipping_polygon polygon;
    jcv_clipper* clipper = 0;

    polygon.num_points = bdry.size();
    polygon.points = convert_to_jcv_points(bdry);

    jcv_clipper polygonclipper;
    polygonclipper.test_fn = jcv_clip_polygon_test_point;
    polygonclipper.clip_fn = jcv_clip_polygon_clip_edge;
    polygonclipper.fill_fn = jcv_clip_polygon_fill_gaps;
    polygonclipper.ctx = &polygon;
    clipper = &polygonclipper;

    convert_to_jcv_points(lane_lines, vc);

    vc.rect = { bounding_box[0].first - 1, bounding_box[0].second - 1, bounding_box[1].first + 1, bounding_box[1].second + 1 };
    memset(&vc.diagram, 0, sizeof(jcv_diagram));

    jcv_diagram_generate(vc.num_points, vc.points, &vc.rect, 0, &vc.diagram);

    /* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */
    int width = 512;
    int height = 512;
    int count = 200;
    int numrelaxations = 0;
    const char* inputfile = 0;
    const char* clipfile = 0; // a file with clipping points
    const char* outputfile = "example.png";

    size_t imagesize = (size_t)(width * height * 3);
    unsigned char* image = (unsigned char*)malloc(imagesize);
    memset(image, 0, imagesize);

    jcv_point dimensions;
    dimensions.x = (jcv_real)width;
    dimensions.y = (jcv_real)height;

    unsigned char color_pt[] = { 255, 64, 64 };
    unsigned char color_line[] = { 220, 220, 220 };
    unsigned char color_delauney[] = { 64, 64, 255 };

    unsigned char color_edge[3];
    unsigned char basecolor = 120;
    
    // If all you need are the edges
    const jcv_edge* edge = jcv_diagram_get_edges(&vc.diagram);
    int edgecount = 1;
    while (edge)
    {
        color_edge[0] = rand() % 255;
        color_edge[1] = rand() % 255;
        color_edge[2] = rand() % 255;


        jcv_point p0 = remap(&edge->pos[0], &vc.diagram.min, &vc.diagram.max, &dimensions);
        jcv_point p1 = remap(&edge->pos[1], &vc.diagram.min, &vc.diagram.max, &dimensions);
        draw_line((int)p0.x, (int)p0.y, (int)p1.x, (int)p1.y, image, width, height, 3, color_edge);
        edge = jcv_diagram_get_next_edge(edge);
        edgecount++;
    }
    std::cout << "edgecount: " << edgecount << std::endl;

    printf("Done.\n"); // rendering

    jcv_diagram_free(&vc.diagram);

    // Plot the sites ---- lane points
    for (int i = 0; i < vc.num_points; ++i)
    {
        jcv_point p = remap(&vc.points[i], &vc.diagram.min, &vc.diagram.max, &dimensions);
        plot((int)p.x, (int)p.y, image, width, height, 3, color_pt);
    }

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
    printf("Writing %s\n", path);

    stbi_write_png(path, width, height, 3, image, stride);
    printf("Done.\n");

    free(image);


    return 0;

}