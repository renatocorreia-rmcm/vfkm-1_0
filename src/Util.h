#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include "PolygonalPath.h"
#include "Vector.h"


/*
    CLUSTER CLASS
*/


typedef struct Cluster{
    std::string name;  // amount of curves
    std::vector<Cluster*> children;
    std::vector<int> indices;
    std::vector<float> curveErrors;
    std::pair<Vector*,Vector*> vectorField;  // { [x1, x2, ...], [y1, y2, ...] }  // each vector is linear, but contains 2 vectors of R^2 elements to represent the whole grid
    Cluster* parent;
    
    float error;
    float maxError;

    void clearChildren();
} Cluster;


/*
    UTIL CLASS
*/


class Util
{
    /*
        static class - no members
        just a container of methods
    */

public:

    Util();

    static void loadCurves(
        std::string filename, std::vector<PolygonalPath>&
    );
    static void loadCurves(
        std::string filename, std::vector<PolygonalPath>&,
        float &xmin, float &xmax, float &ymin, float &ymax, float &tmin, float &tmax  // dataset bounding box (to store, not read)
    );


    /* NOT USED */
    static void loadCurvesAndProject(
        std::string filename, std::vector<PolygonalPath>&,
        float &xmin, float &xmax, float &ymin, float &ymax, float &tmin, float &tmax
    );
    static void to_mercator(const float &lat, const float &lon, float &xMerc, float &yMerc);
};

#endif // UTIL_H
