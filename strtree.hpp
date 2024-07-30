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

    auto start = std::chrono::high_resolution_clock().now();
    // Insert centerlines into the index
    std::vector<std::unique_ptr<Geometry>> geometries;
    for (const auto& line : centerlines) {
        auto lineString = makeLineString(line, factory.get());
        index.insert(lineString->getEnvelopeInternal(), lineString.get());
        geometries.push_back(std::move(lineString));
    }

    auto end = std::chrono::high_resolution_clock().now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "------- Elapsed time -- put in strtree: " << elapsed.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock().now();

    // Create boundary polygon
    auto boundaryPolygon = makePolygon(boundary, factory.get());

    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    std::cout << "------- Elapsed time -- create polygon: " << elapsed.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock().now();

    auto lanelines = makeMultiLineString(lanes, factory.get());

    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    std::cout << "------- Elapsed time -- create multilinestring: " << elapsed.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock().now();

    // Query the index with the boundary polygon
    std::vector<Geometry*> queryResult;
    const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();

    auto visitor = [&queryResult, &boundaryPolygon, &lanelines](Geometry* geom) {
        //if (geom->within(boundaryPolygon)) {
        if (geom->within(boundaryPolygon) && geom->distance(lanelines.get()) >= 1.2) {
            queryResult.push_back(geom);
        }
        //}
    };
    
    index.query(*queryEnv, visitor);


    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    std::cout << "------- Elapsed time -- query: " << elapsed.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock().now();

    // Filter the original centerlines based on query result
    std::vector<std::pair<jcv_point, jcv_point>> filteredCenterlines;

    /*for (const auto& geom : queryResult) {
        const auto* lineString = dynamic_cast<const geos::geom::LineString*>(geom);
        jcv_point p0 = { (float)(lineString->getCoordinatesRO()->getAt(0).x), (float)(lineString->getCoordinatesRO()->getAt(0).y) };
        jcv_point p1 = { (float)(lineString->getCoordinatesRO()->getAt(1).x), (float)(lineString->getCoordinatesRO()->getAt(0).y) };
        filteredCenterlines.emplace_back(std::make_pair(p0, p1));
    }*/
    //std::cout << "------- Filtered centerlines: " << queryResult.size() << std::endl;
    std::vector<std::vector<jcv_point>> centerline;
    /*for (uint32_t i = 0; i < queryResult.size(); i++) {
        auto geom = dynamic_cast<geos::geom::LineString*>(queryResult[i]);
        jcv_point p0 = { (float)(geom->getCoordinatesRO()->getAt(0).x), (float)(geom->getCoordinatesRO()->getAt(0).y) };
        jcv_point p1 = { (float)(geom->getCoordinatesRO()->getAt(1).x), (float)(geom->getCoordinatesRO()->getAt(0).y) };
        centerline.push_back({ p0, p1 });
    }*/
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
    //std::cout << "------- Merged lines 2: " << merged.size() << std::endl;
    /*for (const auto& line : merged) {
        std::cout << line->toString() << std::endl;
    }*/
    
	end = std::chrono::high_resolution_clock().now();

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
 //   std::unique_ptr<geos::geom::MultiLineString> multiLine = makeMultiLineString(ctls, factory.get());

 //   // Use LineMerger to merge the MultiLineString
 //   geos::operation::linemerge::LineMerger lineMerger2;
 //   lineMerger2.add(multiLine.get());
 //   auto merged(lineMerger2.getMergedLineStrings());
 //   std::cout << "------- Merged lines 2: " << merged.size() << std::endl;
 //   for (const auto& line : merged) {
 //       std::cout << line->toString() << std::endl;
 //   }

    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    std::cout << "------- Elapsed time -- filter & merge: " << elapsed.count() << " ms" << std::endl;

    return ctls;
}

