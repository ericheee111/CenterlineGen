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


using namespace geos::geom;
using namespace geos::index::strtree;

// Function to convert jcv_point to GEOS Coordinate
Coordinate toCoordinate(const jcv_point& point) {
    return Coordinate(point.x, point.y);
}

// Function to create a GEOS LineString from a pair of jcv_point
std::unique_ptr<LineString> makeLineString(const std::pair<jcv_point, jcv_point>& line, const GeometryFactory* factory) {
    auto coordseq = CoordinateSequence(); // Create a sequence of coordinates
    coordseq.add(toCoordinate(line.first)); // Add the first point
    coordseq.add(toCoordinate(line.second)); // Add the second point
    return std::unique_ptr<LineString>(factory->createLineString(coordseq));
}

// Function to create a GEOS Polygon from a vector of jcv_point (boundary)
std::unique_ptr<Polygon> makePolygon(const std::vector<jcv_point>& boundary, const GeometryFactory* factory) {
    auto coordseq = CoordinateSequence(); // Create a sequence of coordinates
    for (const auto& point : boundary) {
        coordseq.add(toCoordinate(point)); // Add each point to the sequence
    }
    // Close the polygon by adding the first point at the end if not already closed
    if (!(coordseq.front() == coordseq.back())) {
        coordseq.add(coordseq.front());
    }
    auto polygon = factory->createLinearRing(coordseq);
    return factory->createPolygon(std::move(polygon));
}

// Main function to filter centerlines based on boundary
std::vector<std::pair<jcv_point, jcv_point>> filterCenterlines(
    const std::vector<std::pair<jcv_point, jcv_point>>& centerlines,
    const std::vector<jcv_point>& boundary) {

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

    // Query the index with the boundary polygon
    std::vector<Geometry*> queryResult;
    auto visitor = [&queryResult, &boundaryPolygon](Geometry* geom) {
        if (geom->intersects(boundaryPolygon.get())) {
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


