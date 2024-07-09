#include "voronoi.h"

voronoi::voronoi(int num_points, jcv_point* points, jcv_rect* rect, jcv_clipper* clipper, jcv_diagram* diagram)
{
	m_NumPoints = num_points;
	m_ShapePoints = points;
	m_BoundingRect = rect;
	m_Clipper = clipper;
	m_Diagram = *diagram;
}

void voronoi::run()
{
	jcv_diagram_generate(m_NumPoints, m_ShapePoints, m_BoundingRect, m_Clipper, &m_Diagram);
}
