#pragma once

#include "jc_voronoi.h"
#include "jc_voronoi_clip.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

class voronoi
{
public:
	voronoi(int num_points, jcv_point* points, jcv_rect* rect, jcv_clipper* clipper, jcv_diagram* diagram);
	~voronoi() {};
	void run();

private:
	jcv_diagram m_Diagram;
	jcv_point* m_ShapePoints;
	int m_NumPoints;
	jcv_rect* m_BoundingRect;
	jcv_clipper* m_Clipper;

};

