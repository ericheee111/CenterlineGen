#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

#define MINIMUM_POINTS 2     // minimum number of cluster
#define EPSILON (2*2)  // distance for clustering, metre^2

using namespace std;

typedef struct Points_
{
    double x, y;  // X, Y position
    int32_t clusterID;  // clustered ID
    Points_() : x(0), y(0), clusterID(-1) {}
    Points_(double x, double y) : x(x), y(y), clusterID(-1) {}
}Point;

class DBSCAN {
public:
    DBSCAN(int32_t minPts, double eps, vector<Point> points) {
        m_minPoints = minPts;
        m_epsilon = eps;
        m_points = points;
        m_pointSize = (int32_t)points.size();
        m_clusterCount = 0;
    }
    ~DBSCAN() {}

    int32_t run();
    vector<int32_t> calculateCluster(Point point);
    int32_t expandCluster(Point point, int32_t clusterID);
    inline double calculateDistance(const Point& pointCore, const Point& pointTarget);

    int32_t getTotalPointSize() { return m_pointSize; }
    int32_t getMinimumClusterSize() { return m_minPoints; }
    int32_t getEpsilonSize() { return m_epsilon; }
    int32_t getClusterCount() { return m_clusterCount - 1; }

public:
    vector<Point> m_points;

private:
    int32_t m_pointSize;
    int32_t m_minPoints;
    int32_t m_epsilon;
    int32_t m_clusterCount;
};

#endif // DBSCAN_H