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
#include <iostream>
#include <geos/index/strtree/STRtree.h>
#include <chrono>



using namespace geos::geom;
using namespace geos::index::strtree;
using namespace geos::operation::linemerge;

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

// Main function to filter centerlines based on boundary
std::vector<std::vector<jcv_point>> filterCenterlines(
    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
    const std::vector<jcv_point>& boundary, 
    const std::vector<std::vector<jcv_point>>& lanes) {

    // Create GeometryFactory
    GeometryFactory::Ptr factory = GeometryFactory::create();

    // Create STRtree index
    TemplateSTRtree<Geometry*> index;

    // Insert centerlines into the index
    std::vector<std::unique_ptr<Geometry>> geometries;
    for (const auto& line : centerlines) {
        auto lineString = makeLineString(line, factory.get());
        index.insert(lineString->getEnvelopeInternal(), lineString.get());
        geometries.push_back(std::move(lineString));
    }

    // Create boundary polygon
    auto boundaryPolygon = makePolygon(boundary, factory.get());

    auto lanelines = makeMultiLineString(lanes, factory.get());

    // https://gis.stackexchange.com/questions/212007/using-spatial-index-to-intersect-points-with-polygon-when-points-and-polygon-ha
    // Query the index with the boundary polygon
    std::vector<Geometry*> queryResult;
    const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();
    
    auto start = std::chrono::high_resolution_clock().now();
    auto visitor = [&queryResult, &boundaryPolygon, &lanelines](Geometry* geom) {
        if (geom->getEnvelope()->intersects(boundaryPolygon)) {
            if (lanelines.get()->distance(geom) >= 1.25) {
                queryResult.push_back(geom);
            }

        }
    };
    
    index.query(*queryEnv, visitor);

    auto end = std::chrono::high_resolution_clock().now();
    std::chrono::duration<double, milli> elapsed = end - start;
    std::cout << "------- Elapsed time -- query: " << elapsed.count() << " ms" << std::endl;


    std::vector<std::vector<jcv_point>> centerline;

    for (const auto& geom : queryResult) {
        auto line = dynamic_cast<geos::geom::LineString*>(geom);
        std::vector<jcv_point> points;
        for (int i = 0; i < line->getNumPoints(); i++) {
            auto coord = line->getCoordinateN(i);
            points.push_back({ (float)coord.x, (float)coord.y });
        }
        centerline.push_back(points);
    }
    std::unique_ptr<geos::geom::MultiLineString> multilinestring = makeMultiLineString(centerline, factory.get());
    geos::operation::linemerge::LineMerger lineMerger;
    lineMerger.add(multilinestring.get());
    auto merged(lineMerger.getMergedLineStrings());
    
    std::vector<std::vector<jcv_point>> ctls;

    for (const auto& line : merged) {
        std::vector<jcv_point> points;
        for (int i = 0; i < line->getNumPoints(); i++) {
            auto coord = line->getCoordinateN(i);
			points.push_back({ coord.x, coord.y });
		}
        ctls.push_back(points);
	}
    //savetofile(ctls);

    return ctls;
}

