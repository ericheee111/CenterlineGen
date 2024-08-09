#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
//#include <map>
#include <array>
#include "jc_voronoi.h"

// Hash function for Point with double coordinates
struct PointHash {
    std::size_t operator()(const jcv_point& p) const {
        auto hash1 = std::hash<double>{}(p.x);
        auto hash2 = std::hash<double>{}(p.y);
        return hash1 ^ (hash2 << 1);
    }
};

// Line segment structure
struct LineSeg {
    jcv_point p1, p2;

    // Constructor to ensure consistent ordering of points
    LineSeg(jcv_point point1, jcv_point point2) {
        if (point1.x < point2.x || (point1.x == point2.x && point1.y < point2.y)) {
            p1 = point1;
            p2 = point2;
        }
        else {
            p1 = point2;
            p2 = point1;
        }
    }

    // Equality operator for LineSegment considering unordered points
    bool operator==(const LineSeg& other) const {
        return p1 == other.p1 && p2 == other.p2;
    }
};

struct LineSegmentHash {
    std::size_t operator()(const LineSeg& ls) const {
        PointHash pointHash;
        size_t h1 = pointHash(ls.p1);
        size_t h2 = pointHash(ls.p2);
        return h1 ^ (h2 << 1);
    }
};


class edge_info {
public:
    edge_info() : edge(nullptr), site1(nullptr), site2(nullptr) {}
    edge_info(const jcv_edge* edge, jcv_site* site1, jcv_site* site2) : edge(edge), site1(site1), site2(site2) {}
	const jcv_edge* edge;
	jcv_site* site1;
	jcv_site* site2;
};

class vor_manager
{
public:
    vor_manager() {};
    ~vor_manager() {};

    void addEdge(const jcv_edge* edge, jcv_site* site1, jcv_site* site2);
    std::pair<jcv_site*, jcv_site*> getEdgeSites(jcv_point p1, jcv_point p2);


private:
    std::unordered_map<LineSeg, edge_info, LineSegmentHash> lineMap;
    //std::map<std::array<double,4>, edge_info> lineMap;
};

