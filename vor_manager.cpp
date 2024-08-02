#include "vor_manager.h"

void vor_manager::addEdge(jcv_edge* edge, jcv_site* site1, jcv_site* site2)
{
	LineSegment ls = LineSegment(edge->pos[0], edge->pos[1]);
	edge_info ei = edge_info(edge, site1, site2);
	lineMap[ls] = ei;
}

std::pair<jcv_site*, jcv_site*> vor_manager::getEdgeSites(jcv_point p1, jcv_point p2)
{
	LineSegment ls = LineSegment(p1, p2);
	auto edge = lineMap[ls];
	return std::make_pair(edge.site1, edge.site2);
}
