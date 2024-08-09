#include <iostream>
#include <vector>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/Polygon.h>
#include <geos/index/strtree/TemplateSTRtree.h>
#include "jc_voronoi.h"

#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/operation/linemerge/LineMerger.h>
#include <geos/operation/overlay/OverlayOp.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include <geos/geom/prep/PreparedPolygon.h>
#include <iostream>
#include <geos/index/strtree/STRtree.h>
#include <chrono>



using namespace geos::geom;
using namespace geos::index::strtree;
using namespace geos::operation::linemerge;
using namespace geos::operation::overlay;
using namespace geos::operation::distance;
using namespace geos::geom::prep;

// Function to convert jcv_point to GEOS Coordinate
Coordinate toCoordinate(const jcv_point& point) {
    return Coordinate(point.x, point.y);
}

// Function to create a GEOS LineString from a pair of jcv_point
std::unique_ptr<LineString> makeLineString(const std::pair<jcv_point, jcv_point>& line, const GeometryFactory* factory) {
    auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
    coordseq.add(toCoordinate(line.first)); // Add the first point
    coordseq.add(toCoordinate(line.second)); // Add the second point
    return std::unique_ptr<LineString>(factory->createLineString(coordseq));
}

// Function to create a GEOS Polygon from a vector of jcv_point (boundary)
geos::geom::Polygon* makePolygon(const std::vector<jcv_point>& boundary, const GeometryFactory* factory) {
    auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
    for (const auto& point : boundary) {
        coordseq.add(toCoordinate(point)); // Add each point to the sequence
    }
    // Close the polygon by adding the first point at the end if not already closed
    if (!(coordseq.front() == coordseq.back())) {
        coordseq.add(coordseq.front());
    }
    auto polygon = factory->createLinearRing(coordseq);
    return factory->createPolygon(polygon, nullptr);
}

// Function to create a GEOS MultiLineString from a vector of vector of jcv_point
std::unique_ptr<MultiLineString> makeMultiLineString(const std::vector<std::vector<jcv_point>>& lines, const GeometryFactory* factory) {
    std::vector<std::unique_ptr<LineString>> lineStrings;
    for (const auto& line : lines) {
        auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
        for (const auto& point : line) {
            coordseq.add(toCoordinate(point)); // Add each point to the sequence
        }
        lineStrings.push_back(std::unique_ptr<LineString>(factory->createLineString(coordseq)));
    }
    return std::unique_ptr<MultiLineString>(factory->createMultiLineString(std::move(lineStrings)));
}

std::unique_ptr<MultiLineString> makeMultiLineString(const std::vector<std::pair<jcv_point, jcv_point>>& lines, const GeometryFactory* factory) {
    std::vector<std::unique_ptr<LineString>> lineStrings;
    for (const auto& line : lines) {
        auto coordseq = CoordinateArraySequence(); // Create a sequence of coordinates
        coordseq.add(toCoordinate({ (float)line.first.x, (float)line.first.y })); // Add each point to the sequence
        coordseq.add(toCoordinate({ (float)line.second.x, (float)line.second.y })); // Add each point to the sequence
        lineStrings.push_back(std::unique_ptr<LineString>(factory->createLineString(coordseq)));
    }
    return std::unique_ptr<MultiLineString>(factory->createMultiLineString(std::move(lineStrings)));
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

std::vector<std::unique_ptr<Geometry>> createGridCells(const Geometry* boundary, int gridSize, const GeometryFactory* factory) {
    std::vector<std::unique_ptr<Geometry>> gridCells;
    const Envelope* env = boundary->getEnvelopeInternal();
    double minX = env->getMinX();
    double minY = env->getMinY();
    double maxX = env->getMaxX();
    double maxY = env->getMaxY();
    double cellWidth = (maxX - minX) / gridSize;
    double cellHeight = (maxY - minY) / gridSize;

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            double x1 = minX + i * cellWidth;
            double y1 = minY + j * cellHeight;
            double x2 = x1 + cellWidth;
            double y2 = y1 + cellHeight;
            Envelope cellEnv(x1, x2, y1, y2);
            std::unique_ptr<Geometry> cellPolygon(factory->toGeometry(&cellEnv));
            std::unique_ptr<Geometry> clippedPolygon(OverlayOp::overlayOp(cellPolygon.get(), boundary, OverlayOp::OpCode::opINTERSECTION));
            if (!clippedPolygon->isEmpty()) {
                gridCells.push_back(std::move(clippedPolygon));
            }
        }
    }

    return gridCells;
}

void splitAndIndexLanelines(const std::vector<std::unique_ptr<Geometry>>& gridCells, const std::vector<std::unique_ptr<Geometry>>& lanelines, TemplateSTRtree<const Geometry*>& laneStrTree) {
    GeometryFactory::Ptr factory = GeometryFactory::create();
    for (const auto& cell : gridCells) {
        for (const auto& laneline : lanelines) {
            std::unique_ptr<Geometry> clippedLine(OverlayOp::overlayOp(laneline.get(), cell.get(), OverlayOp::OpCode::opINTERSECTION));
            if (!clippedLine->isEmpty()) {
                laneStrTree.insert(clippedLine->getEnvelopeInternal(), clippedLine.release());
            }
        }
    }
}

bool isPointInBoundary(const Point* point, TemplateSTRtree<const Geometry*>& strtree) {
    std::vector<void*> results;
    strtree.query(point->getEnvelopeInternal(), results);

    for (const auto& cell : results) {
        if (static_cast<Geometry*>(cell)->contains(point)) {
            return true;
        }
    }
    return false;
}

bool isCenterlineTooCloseToLanelines(const Geometry* centerline, TemplateSTRtree<const Geometry*>& laneStrTree, double minDistance) {
    std::vector<void*> results;
    laneStrTree.query(centerline->getEnvelopeInternal(), results);

    for (const auto& laneline : results) {
        if (DistanceOp::distance(centerline, static_cast<Geometry*>(laneline)) < minDistance) {
            return true;
        }
    }
    return false;
}

std::unique_ptr<GeometryCollection> makeGeometryCollection(const std::vector<std::vector<jcv_point>>& lanes, const GeometryFactory* factory) {
    std::vector<std::unique_ptr<Geometry>> geometries;
    for (const auto& lane : lanes) {
        for (const auto& point : lane) {
            geometries.push_back(std::unique_ptr<Geometry>(factory->createPoint(toCoordinate(point))));
        }
    }
    return std::unique_ptr<GeometryCollection>(factory->createGeometryCollection(std::move(geometries)));
}


// Main function to filter centerlines based on boundary
std::vector<std::vector<jcv_point>> filterCenterlines(
    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
    const std::vector<jcv_point>& boundary, 
    const std::vector<std::vector<jcv_point>>& lanes) {

    // Create GeometryFactory
    GeometryFactory::Ptr factory = GeometryFactory::create();

    //// Create STRtree index
    //TemplateSTRtree<Geometry*> index;

    //// Insert centerlines into the index
    //std::vector<std::unique_ptr<Geometry>> geometries;
    //for (const auto& line : centerlines) {
    //    auto lineString = makeLineString(line, factory.get());
    //    index.insert(lineString->getEnvelopeInternal(), lineString.get());
    //    geometries.push_back(std::move(lineString));
    //}

    //// Create boundary polygon
    //auto boundaryPolygon = makePolygon(boundary, factory.get());

    //auto preparedBoundaryPolygon = PreparedGeometryFactory::prepare(boundaryPolygon);

    //auto lanelines = makeMultiLineString(lanes, factory.get());
    ////auto lanelines = makeGeometryCollection(lanes, factory.get());
    //auto preparedLanes = PreparedGeometryFactory::prepare(lanelines.get());

    //// https://gis.stackexchange.com/questions/212007/using-spatial-index-to-intersect-points-with-polygon-when-points-and-polygon-ha
    //// Query the index with the boundary polygon
    //std::vector<Geometry*> queryResult;
    //const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();
    //
    //auto start = std::chrono::high_resolution_clock().now();
    //auto visitor = [&queryResult, &preparedBoundaryPolygon, &preparedLanes](Geometry* geom) {
    //    /*if (geom->getEnvelope()->intersects(boundaryPolygon)) {
    //        if (lanelines.get()->distance(geom) >= 1.25) {
    //            queryResult.push_back(geom);
    //        }

    //    }*/
    //    if (preparedBoundaryPolygon.get()->containsProperly(geom->getEnvelope().get())) {
    //        if (preparedLanes.get()->distance(geom) >= 1.25) {
    //            queryResult.push_back(geom);
    //        }

    //    }
    //};
    //
    //index.query(*queryEnv, visitor);

    //auto end = std::chrono::high_resolution_clock().now();
    //std::chrono::duration<double, milli> elapsed = end - start;
    //std::cout << "------- Elapsed time -- query: " << elapsed.count() << " ms" << std::endl;


    //std::vector<std::vector<jcv_point>> centerline;

    //for (const auto& geom : queryResult) {
    //    auto line = dynamic_cast<geos::geom::LineString*>(geom);
    //    std::vector<jcv_point> points;
    //    for (int i = 0; i < line->getNumPoints(); i++) {
    //        auto coord = line->getCoordinateN(i);
    //        points.push_back({ (float)coord.x, (float)coord.y });
    //    }
    //    centerline.push_back(points);
    //}
    //std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerline, factory.get());
    std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerlines, factory.get());
    geos::operation::linemerge::LineMerger lineMerger;
    lineMerger.add(multilinestring.get());
    auto merged(lineMerger.getMergedLineStrings());
    
    std::vector<std::vector<jcv_point>> ctls;

    for (const auto& line : merged) {
        std::vector<jcv_point> points;
        for (int i = 0; i < line->getNumPoints(); i++) {
            auto coord = line->getCoordinateN(i);
			points.push_back({ (float)coord.x, (float)coord.y });
		}
        ctls.push_back(points);
	}
    //savetofile(ctls);

    return ctls;
}

std::vector<std::vector<jcv_point>> filterCenterlinesByBoundary(
    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
    const std::vector<jcv_point>& boundary) {

    // Create GeometryFactory
    GeometryFactory::Ptr factory = GeometryFactory::create();

    //// Create STRtree index
    //TemplateSTRtree<Geometry*> index;

    //// Insert centerlines into the index
    //std::vector<std::unique_ptr<Geometry>> geometries;
    //for (const auto& line : centerlines) {
    //    auto lineString = makeLineString(line, factory.get());
    //    index.insert(lineString->getEnvelopeInternal(), lineString.get());
    //    geometries.push_back(std::move(lineString));
    //}

    //// Create boundary polygon
    //auto boundaryPolygon = makePolygon(boundary, factory.get());

    //auto preparedBoundaryPolygon = PreparedGeometryFactory::prepare(boundaryPolygon);

    //// https://gis.stackexchange.com/questions/212007/using-spatial-index-to-intersect-points-with-polygon-when-points-and-polygon-ha
    //// Query the index with the boundary polygon
    //std::vector<Geometry*> queryResult;
    //const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();
    //
    //auto start = std::chrono::high_resolution_clock().now();
    //auto visitor = [&queryResult, &preparedBoundaryPolygon](Geometry* geom) {
    //    if (preparedBoundaryPolygon.get()->containsProperly(geom->getEnvelope().get())) {
    //        queryResult.push_back(geom);
    //    }
    //};
    //
    //index.query(*queryEnv, visitor);

    //auto end = std::chrono::high_resolution_clock().now();
    //std::chrono::duration<double, milli> elapsed = end - start;
    //std::cout << "------- Elapsed time -- query: " << elapsed.count() << " ms" << std::endl;


    //std::vector<std::vector<jcv_point>> centerline;

    //for (const auto& geom : queryResult) {
    //    auto line = dynamic_cast<geos::geom::LineString*>(geom);
    //    std::vector<jcv_point> points;
    //    for (int i = 0; i < line->getNumPoints(); i++) {
    //        auto coord = line->getCoordinateN(i);
    //        points.push_back({ (float)coord.x, (float)coord.y });
    //    }
    //    centerline.push_back(points);
    //}
    //std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerline, factory.get());
    std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerlines, factory.get());
    geos::operation::linemerge::LineMerger lineMerger;
    lineMerger.add(multilinestring.get());
    auto merged(lineMerger.getMergedLineStrings());
    
    std::vector<std::vector<jcv_point>> ctls;

    for (const auto& line : merged) {
        std::vector<jcv_point> points;
        for (int i = 0; i < line->getNumPoints(); i++) {
            auto coord = line->getCoordinateN(i);
			points.push_back({ (float)coord.x, (float)coord.y });
		}
        ctls.push_back(points);
	}
    //savetofile(ctls);

    return ctls;
}

