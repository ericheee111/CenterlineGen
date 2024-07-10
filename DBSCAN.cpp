#include "DBSCAN.h"

int32_t DBSCAN::run()
{
    int32_t clusterID = 1;
    vector<PointDB>::iterator iter;
    for (iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if (iter->clusterID == UNCLASSIFIED)
        {
            if (expandCluster(*iter, clusterID) != FAILURE)
            {
                clusterID += 1;
            }
        }
    }
    m_clusterCount = clusterID;

    return 0;
}

int32_t DBSCAN::expandCluster(PointDB point, int32_t clusterID)
{
    vector<int32_t> clusterSeeds = calculateCluster(point);

    if (clusterSeeds.size() < m_minPoints)
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else
    {
        int32_t index = 0, indexCorePoint = 0;
        vector<int32_t>::iterator iterSeeds;
        for (iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            m_points.at(*iterSeeds).clusterID = clusterID;
            if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y)
            {
                indexCorePoint = index;
            }
            ++index;
        }
        clusterSeeds.erase(clusterSeeds.begin() + indexCorePoint);

        for (vector<int32_t>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i)
        {
            vector<int32_t> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));

            if (clusterNeighors.size() >= m_minPoints)
            {
                vector<int32_t>::iterator iterNeighors;
                for (iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors)
                {
                    if (m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE)
                    {
                        if (m_points.at(*iterNeighors).clusterID == UNCLASSIFIED)
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return SUCCESS;
    }
}

vector<int32_t> DBSCAN::calculateCluster(PointDB point)
{
    int32_t index = 0;
    vector<PointDB>::iterator iter;
    vector<int32_t> clusterIndex;
    for (iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if (calculateDistance(point, *iter) <= m_epsilon)
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    return clusterIndex;
}

inline double DBSCAN::calculateDistance(const PointDB& pointCore, const PointDB& pointTarget)
{
    return pow(pointCore.x - pointTarget.x, 2) + pow(pointCore.y - pointTarget.y, 2);
}

