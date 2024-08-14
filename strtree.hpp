#include <iostream>
#include <vector>
//#include <geos/geom/Coordinate.h>
//#include <geos/geom/Envelope.h>
//#include <geos/geom/GeometryFactory.h>
//#include <geos/geom/Geometry.h>
//#include <geos/geom/Point.h>
//#include <geos/geom/Polygon.h>
//#include <geos/index/strtree/TemplateSTRtree.h>
#include "jc_voronoi.h"

//#include <geos/geom/CoordinateSequence.h>
//#include <geos/geom/CoordinateArraySequence.h>
//#include <geos/geom/GeometryFactory.h>
//#include <geos/geom/LineString.h>
//#include <geos/operation/linemerge/LineMerger.h>
//#include <geos/operation/overlay/OverlayOp.h>
//#include <geos/operation/distance/DistanceOp.h>
//#include <geos/geom/prep/PreparedGeometryFactory.h>
//#include <geos/geom/prep/PreparedPolygon.h>
#include <iostream>
//#include <geos/index/strtree/STRtree.h>
#include <chrono>
#include <geos_c.h>
#include <cassert>
#include "vor_manager.h"



//using namespace geos::geom;
//using namespace geos::index::strtree;
//using namespace geos::operation::linemerge;
//using namespace geos::operation::overlay;
//using namespace geos::operation::distance;
//using namespace geos::geom::prep;

//// Function to convert jcv_point to GEOS Coordinate
//Coordinate toCoordinate(const jcv_point& point) {
//    return Coordinate(point.x, point.y);
//}
//
//// Function to create a GEOS LineString from a pair of jcv_point
//std::unique_ptr<LineString> makeLineString(const std::pair<jcv_point, jcv_point>& line, const GeometryFactory* factory) {
//    auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
//    coordseq.add(toCoordinate(line.first)); // Add the first point
//    coordseq.add(toCoordinate(line.second)); // Add the second point
//    return std::unique_ptr<LineString>(factory->createLineString(coordseq));
//}
//
//// Function to create a GEOS Polygon from a vector of jcv_point (boundary)
//geos::geom::Polygon* makePolygon(const std::vector<jcv_point>& boundary, const GeometryFactory* factory) {
//    auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
//    for (const auto& point : boundary) {
//        coordseq.add(toCoordinate(point)); // Add each point to the sequence
//    }
//    // Close the polygon by adding the first point at the end if not already closed
//    if (!(coordseq.front() == coordseq.back())) {
//        coordseq.add(coordseq.front());
//    }
//    auto polygon = factory->createLinearRing(coordseq);
//    return factory->createPolygon(polygon, nullptr);
//}
//
//// Function to create a GEOS MultiLineString from a vector of vector of jcv_point
//std::unique_ptr<MultiLineString> makeMultiLineString(const std::vector<std::vector<jcv_point>>& lines, const GeometryFactory* factory) {
//    std::vector<std::unique_ptr<LineString>> lineStrings;
//    for (const auto& line : lines) {
//        auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
//        for (const auto& point : line) {
//            coordseq.add(toCoordinate(point)); // Add each point to the sequence
//        }
//        lineStrings.push_back(std::unique_ptr<LineString>(factory->createLineString(coordseq)));
//    }
//    return std::unique_ptr<MultiLineString>(factory->createMultiLineString(std::move(lineStrings)));
//}
//
//std::unique_ptr<MultiLineString> makeMultiLineString(const std::vector<std::pair<jcv_point, jcv_point>>& lines, const GeometryFactory* factory) {
//    std::vector<std::unique_ptr<LineString>> lineStrings;
//    for (const auto& line : lines) {
//        auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
//        auto pt0 = toCoordinate({ (float)line.first.x, (float)line.first.y });
//        auto pt1 = toCoordinate({ (float)line.second.x, (float)line.second.y });
//        
//        coordseq.add(pt0); // Add each point to the sequence
//        coordseq.add(pt1); // Add each point to the sequence
//        
//        auto lstr = factory->createLineString(coordseq);
//        //auto lstrr = lstr.get();
//        for (size_t i = 0; i < lstr->getNumPoints(); i++) {
//            auto pt = lstr->getPointN(i);
//            int x = 123;
//            pt.get()->setUserData((void*)&x);
//        }
//        // xxxxxxxxxxxxxxx
//        //auto tstlstr = lstr.get();
//        for (size_t i = 0; i < lstr->getNumPoints(); i++) {
//            auto pt = lstr->getPointN(i);
//            auto dt = pt.get()->getUserData();
//            std::cout << "userdata: " << *(int*)dt << std::endl;
//        }
//        //x xxxxxxxxxxxxxxxxxxxxxxxxxxxx
//        lineStrings.push_back(std::unique_ptr<LineString>(lstr));
//    }
//    return std::unique_ptr<MultiLineString>(factory->createMultiLineString(std::move(lineStrings)));
//}

GEOSGeometry* makeMultiLineString(const std::vector<std::pair<jcv_point, jcv_point>>& lines, GEOSContextHandle_t ctx) {
    std::vector<GEOSGeometry*> lineStrings;

    for (const auto& line : lines) {
        // Create a coordinate sequence for each line
        GEOSCoordSequence* coordseq = GEOSCoordSeq_create_r(ctx, 2, 2);  // 2 points with 2 dimensions (x, y)

        // Set coordinates for the start and end points
        GEOSCoordSeq_setX_r(ctx, coordseq, 0, (float)line.first.x);
        GEOSCoordSeq_setY_r(ctx, coordseq, 0, (float)line.first.y);
        GEOSCoordSeq_setX_r(ctx, coordseq, 1, (float)line.second.x);
        GEOSCoordSeq_setY_r(ctx, coordseq, 1, (float)line.second.y);

        assert(coordseq != NULL);

        // Create the LineString geometry from the coordinate sequence
        GEOSGeometry* lineString = GEOSGeom_createLineString_r(ctx, coordseq);

        assert(lineString != NULL);

        lineStrings.push_back(lineString);
    }

    //std::cout << "Number of lines: " << lineStrings.size() << std::endl;

    // Create the MultiLineString from the array of LineStrings
    GEOSGeometry* multiLineString = GEOSGeom_createCollection_r(ctx, GEOS_MULTILINESTRING, lineStrings.data(), lineStrings.size());

    assert(multiLineString != NULL);
    //auto xx = GEOSGetNumGeometries(multiLineString);
    //std::cout << "Number of geometries: " << xx << std::endl;

    return multiLineString;  // Caller is responsible for freeing this
}


void savetofile(std::vector<std::vector<jcv_point>> lines) {
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

//std::unique_ptr<GeometryCollection> makeGeometryCollection(const std::vector<std::vector<jcv_point>>& lanes, const GeometryFactory* factory) {
//    std::vector<std::unique_ptr<Geometry>> geometries;
//    for (const auto& lane : lanes) {
//        for (const auto& point : lane) {
//            geometries.push_back(std::unique_ptr<Geometry>(factory->createPoint(toCoordinate(point))));
//        }
//    }
//    return std::unique_ptr<GeometryCollection>(factory->createGeometryCollection(std::move(geometries)));
//}

std::vector<std::vector<jcv_point>> MergeCenterlines(const std::vector<std::pair<jcv_point, jcv_point>>& centerlines) {
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, NULL);
    GEOSContext_setErrorHandler_r(ctx, NULL);

    GEOSGeometry* multiLineString = makeMultiLineString(centerlines, ctx);
    assert(multiLineString != NULL);
    GEOSGeometry* merged = GEOSLineMerge_r(ctx, multiLineString);
    assert(merged != NULL);

    std::vector<std::vector<jcv_point>> ctls;

    int num_geometries = GEOSGetNumGeometries_r(ctx, merged);
    //std::cout << "merged num geometries: " << num_geometries << std::endl;
    for (int i = 0; i < num_geometries; i++) {
		auto line = GEOSGeom_getCoordSeq_r(ctx, GEOSGetGeometryN_r(ctx, merged, i));
		std::vector<jcv_point> points;
        unsigned int size, size1;
		GEOSCoordSeq_getSize_r(ctx, line, &size);
        //size1 = GEOSGetNumGeometries_r(ctx, line1);
        //std::cout << "line size: " << size << std::endl;
		for (int j = 0; j < size; j++) {
			double x, y;
			GEOSCoordSeq_getX_r(ctx, line, j, &x);
			GEOSCoordSeq_getY_r(ctx, line, j, &y);

			points.push_back({ (float)x, (float)y });
		}

		ctls.push_back(points);
        //std::cout << "line[" << i << "] size: " << points.size() << std::endl;
	}

    GEOSGeom_destroy_r(ctx, merged);
    GEOSGeom_destroy_r(ctx, multiLineString);

    // Finish GEOS
    GEOS_finish_r(ctx);
    return ctls;
}

// Main function to filter centerlines based on boundary
//std::vector<std::vector<jcv_point>> filterCenterlines(
//    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
//    const std::vector<jcv_point>& boundary, 
//    const std::vector<std::vector<jcv_point>>& lanes) {
//
//    // Create GeometryFactory
//    //GeometryFactory::Ptr factory = GeometryFactory::create();
//
//    //// Create STRtree index
//    //TemplateSTRtree<Geometry*> index;
//
//    //// Insert centerlines into the index
//    //std::vector<std::unique_ptr<Geometry>> geometries;
//    //for (const auto& line : centerlines) {
//    //    auto lineString = makeLineString(line, factory.get());
//    //    index.insert(lineString->getEnvelopeInternal(), lineString.get());
//    //    geometries.push_back(std::move(lineString));
//    //}
//
//    //// Create boundary polygon
//    //auto boundaryPolygon = makePolygon(boundary, factory.get());
//
//    //auto preparedBoundaryPolygon = PreparedGeometryFactory::prepare(boundaryPolygon);
//
//    //auto lanelines = makeMultiLineString(lanes, factory.get());
//    ////auto lanelines = makeGeometryCollection(lanes, factory.get());
//    //auto preparedLanes = PreparedGeometryFactory::prepare(lanelines.get());
//
//    //// https://gis.stackexchange.com/questions/212007/using-spatial-index-to-intersect-points-with-polygon-when-points-and-polygon-ha
//    //// Query the index with the boundary polygon
//    //std::vector<Geometry*> queryResult;
//    //const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();
//    //
//    //auto start = std::chrono::high_resolution_clock().now();
//    //auto visitor = [&queryResult, &preparedBoundaryPolygon, &preparedLanes](Geometry* geom) {
//    //    /*if (geom->getEnvelope()->intersects(boundaryPolygon)) {
//    //        if (lanelines.get()->distance(geom) >= 1.25) {
//    //            queryResult.push_back(geom);
//    //        }
//
//    //    }*/
//    //    if (preparedBoundaryPolygon.get()->containsProperly(geom->getEnvelope().get())) {
//    //        if (preparedLanes.get()->distance(geom) >= 1.25) {
//    //            queryResult.push_back(geom);
//    //        }
//
//    //    }
//    //};
//    //
//    //index.query(*queryEnv, visitor);
//
//    //auto end = std::chrono::high_resolution_clock().now();
//    //std::chrono::duration<double, milli> elapsed = end - start;
//    //std::cout << "------- Elapsed time -- query: " << elapsed.count() << " ms" << std::endl;
//
//
//    //std::vector<std::vector<jcv_point>> centerline;
//
//    //for (const auto& geom : queryResult) {
//    //    auto line = dynamic_cast<geos::geom::LineString*>(geom);
//    //    std::vector<jcv_point> points;
//    //    for (int i = 0; i < line->getNumPoints(); i++) {
//    //        auto coord = line->getCoordinateN(i);
//    //        points.push_back({ (float)coord.x, (float)coord.y });
//    //    }
//    //    centerline.push_back(points);
//    //}
//    //std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerline, factory.get());
//    std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerlines, factory.get());
//    geos::operation::linemerge::LineMerger lineMerger;
//    lineMerger.add(multilinestring.get());
//    auto merged(lineMerger.getMergedLineStrings());
//    
//    std::vector<std::vector<jcv_point>> ctls;
//
//    for (const auto& line : merged) {
//        std::vector<jcv_point> points;
//        for (int i = 0; i < line->getNumPoints(); i++) {
//            auto coord = line->getCoordinateN(i);
//			points.push_back({ (float)coord.x, (float)coord.y });
//		}
//        ctls.push_back(points);
//	}
//    //savetofile(ctls);
//
//    return ctls;
//}
//
//std::vector<std::vector<jcv_point>> filterCenterlinesByBoundary(
//    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
//    const std::vector<jcv_point>& boundary) {
//
//    // Create GeometryFactory
//    GeometryFactory::Ptr factory = GeometryFactory::create();
//
//    std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerlines, factory.get());
//    geos::operation::linemerge::LineMerger lineMerger;
//    lineMerger.add(multilinestring.get());
//    auto merged(lineMerger.getMergedLineStrings());
//    
//    std::vector<std::vector<jcv_point>> ctls;
//
//    for (const auto& line : merged) {
//        std::vector<jcv_point> points;
//        for (int i = 0; i < line->getNumPoints(); i++) {
//            auto coord = line->getCoordinateN(i);
//			points.push_back({ (float)coord.x, (float)coord.y });
//		}
//        ctls.push_back(points);
//	}
//    //savetofile(ctls);
//
//    return ctls;
//}

