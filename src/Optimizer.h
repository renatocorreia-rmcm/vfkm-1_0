#ifndef OPTIMIZER_H
#define OPTIMIZER_H

//#include <vector>
#include "Grid.h"



struct ProblemSettings
{
    // struct carrying problem data for solvers.
    // - curveIndices: respective to cluster
    // - totalCurveLength: summed length of all curves (used for normalization)
    
    Grid &grid;
    const std::vector<int> &curveIndices;
    const std::vector<CurveDescription> &curve_descriptions;
    float totalCurveLength;
    float smoothnessWeight;
    ProblemSettings(Grid &g,
                    const std::vector<int> &i,
                    const std::vector<CurveDescription> &cd,
                    float tcl, float sw):
        grid(g),
        curveIndices(i),
        curve_descriptions(cd),
        totalCurveLength(tcl),
        smoothnessWeight(sw) {}
};



class Optimizer
{
public:
    Optimizer(int size);
    ~Optimizer();

    // Multiply vector(s) by the system matrix A used in the conjugate gradient solver.
    // The system combines a Laplacian-based smoothness term and data terms coming from curve constraints.
    // This variant uses a ProblemSettings struct to access grid/curve parameters.
    
    // Inputs:
    //   x        - input vector (length = number of grid vertices)
    //   resultX  - output (A * x)
    //   diagM    - workspace / diagonal estimate (may be overwritten)
    //   prob     - problem settings 
    static void multiplyByA(const Vector& x, Vector &resultX, ProblemSettings &prob);

    void optimizeImplicitFastWithWeights(Grid &grid, int numberOfVectorFields,
                                         std::vector<PolygonalPath> curves,
                                         std::vector< std::pair<Vector*, Vector*> >& vectorFields,
                                         unsigned short *mapCurveToVectorField,
                                         float *mapCurveToError,
                                         float smoothnessWeight = 0.5);
};

#endif // OPTIMIZER_H
