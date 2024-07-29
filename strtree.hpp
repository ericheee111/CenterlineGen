#include <iostream>
#include <vector>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>
#include <geos/index/strtree/TemplateSTRtree.h>
#include <geos/geom/CoordinateSequence.h>
#include "jc_voronoi.h"

#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <iostream>
#include <geos/index/strtree/STRtree.h>
#include "geos_c.h"


using namespace geos::geom;
using namespace geos::index::strtree;

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


// Main function to filter centerlines based on boundary
std::vector<std::pair<jcv_point, jcv_point>> filterCenterlines(
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

    // Query the index with the boundary polygon
    std::vector<Geometry*> queryResult;
    auto visitor = [&queryResult, &boundaryPolygon, &lanelines](Geometry* geom) {
        /*if (geom->within(boundaryPolygon) && !geom->intersects(lanelines.get())) {
            queryResult.push_back(geom);
            
        }*/
        if (geom->within(boundaryPolygon) && geom->distance(lanelines.get()) >= 1.2) {
            queryResult.push_back(geom);

        }
    };
    const Envelope* queryEnv = boundaryPolygon->getEnvelopeInternal();
    index.query(*queryEnv, visitor);

    // Filter the original centerlines based on query result
    std::vector<std::pair<jcv_point, jcv_point>> filteredCenterlines;
    for (const auto& geom : queryResult) {
        for (const auto& line : centerlines) {
            auto lineString = makeLineString(line, factory.get());
            if (lineString->equals(geom)) {
                filteredCenterlines.push_back(line);
                break;
            }
        }
    }

    return filteredCenterlines;
}

