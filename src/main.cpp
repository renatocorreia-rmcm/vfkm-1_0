#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <map>
#include <cassert>

#include "Optimizer.h"
#include "Util.h"

using namespace std;

void initExperiment(
    string filename, Cluster*& rootCluster, Grid*& g, vector<PolygonalPath>& curves, int gridResolution
)
{
    /*
    initialize root cluster, grid(using gridResolution), curves(using filename)
    */


    //initialize curves
    float xmin, xmax, ymin, ymax, tmin, tmax;  // bounding box of dataset
    Util::loadCurves(filename, curves, xmin, xmax, ymin, ymax, tmin, tmax);  // load `curves` array with each curveContent: vector<pair<2-uple,float>>  // load bounding box

    //initialize grid
    g = new Grid(xmin,ymin,xmax-xmin,ymax-ymin,gridResolution, gridResolution);  // always square grids

    // assign root cluster
    rootCluster = new Cluster;  // initial, top-level cluster that contains all of the curves (trajectories) before the clustering algorithm begins to partition them
    stringstream ss;
    ss << curves.size();
    rootCluster->name = ss.str();
    rootCluster->parent = NULL;

    // initialize rootcluster VF axis
    Vector* rootVFX = new Vector(gridResolution * gridResolution);
    rootVFX->setValues(0.0);
    Vector* rootVFY = new Vector(gridResolution * gridResolution);
    rootVFY->setValues(0.0);

    // vector field is represent by one array for each axis (an array of axis)
    rootCluster->vectorField = {rootVFX, rootVFY};

    // initialize indices and curveErrors
    for(size_t i = 0 ; i < curves.size() ; ++i){
        rootCluster->indices.push_back(i);
        rootCluster->curveErrors.push_back(0.0f);
    }
}

void saveExperiment(string directory, string currentFileLoaded, Cluster* root){
    stringstream ss;
    ss << directory << "/experiment.txt";

    ofstream experimentFile(ss.str().c_str());
    experimentFile << currentFileLoaded << endl;

    // write hierarchy
    queue<Cluster*> nodesToProcess;
    nodesToProcess.push(root);

    //hack
    map<Cluster*, string> mapClusterPath;
    mapClusterPath[root] = "r";

    experimentFile << "-1 r " << endl;

    while(!nodesToProcess.empty()){
        // take next cluster
        Cluster* c = nodesToProcess.front();
        nodesToProcess.pop();

        //write curve indices file
        string clusterName = mapClusterPath[c];

        stringstream curveFileName;
        curveFileName << directory << "/curves_" << clusterName << ".txt";

        ofstream curveIndicesFile(curveFileName.str().c_str());

        //cout << "NumCurves " << numberOfCurves << " errors " << c->curveErrors.size() << endl;
        int numberOfCurves = c->indices.size();
        assert(numberOfCurves == (int)c->curveErrors.size());

        for(int i = 0 ; i < numberOfCurves ; ++i){
            curveIndicesFile << c->indices.at(i) << " " << c->curveErrors.at(i) << endl;
        }
        curveIndicesFile.close();

        //write vector field file
        stringstream vectorFieldFileName;
        vectorFieldFileName << directory << "/vf_" << clusterName << ".txt";
        ofstream vectorFieldFile(vectorFieldFileName.str().c_str());
        Vector* xComponent = c->vectorField.first;
        Vector* yComponent = c->vectorField.second;
        int gridDimension = xComponent->getDimension();

        vectorFieldFile << gridDimension << endl;

        for(int i = 0 ; i < gridDimension ; ++i){
            vectorFieldFile << xComponent[0][i] << " " << yComponent[0][i] << endl;
        }
        vectorFieldFile.close();

        //process children
        int numberOfChildren = c->children.size();

        for(int i = 0 ; i < numberOfChildren ; ++i){
            Cluster* child = c->children.at(i);

            stringstream ss;
            ss << clusterName << "_" <<  i;
            mapClusterPath[child] = ss.str();

            experimentFile << clusterName << " " << ss.str() << endl;

            nodesToProcess.push(child);
        }
    }

}

int main(int argc, char *argv[]){
    /*
    arguments count
    arguments values
    */

    int rightNumberOfParameters = 6;
    if(argc != rightNumberOfParameters){
        //print correct usage
        cout << "./vfkm trajectoryFile gridResolution numberOfVectorFields smoothnessWeight outputDirectory" << endl;
        return 0;
    }

    // load arguments
    string filename(argv[1]);
    int    gridResolution = atoi(argv[2]);
    int    numberOfVectorFields = atoi(argv[3]);
    float  smoothnessWeight = atof(argv[4]);
    string outputDirectory(argv[5]);

    // initialize parameters
    vector<PolygonalPath> curves;
    Cluster* rootCluster = NULL;
    Grid* g = NULL;
    
    // load files + assign parameters
    cout << "Loading Files into parameters..." << endl;
    initExperiment(filename, rootCluster, g, curves, gridResolution);

    //optimize
    cout << "Optimizing..." << endl;
    Cluster* currentCluster = rootCluster;

    Optimizer op(g->getResolutionX() * g->getResolutionY());  // 
    int numberOfCurves = currentCluster->indices.size();

    // intialize maps
    unsigned short mapCurveToVF[numberOfCurves];
    float mapCurveToError[numberOfCurves];
    unsigned int mapCurveToIndexInCurveVector[numberOfCurves];

    vector<PolygonalPath> curvesInCurrentCluster;
    
    for(int i = 0 ; i < numberOfCurves ; ++i){
        /*
        here `current cluster` is the `root cluster`
        */
       mapCurveToError[i] = 0;
       mapCurveToVF[i] = -1;
       
       mapCurveToIndexInCurveVector[i] = currentCluster->indices.at(i);
       curvesInCurrentCluster.push_back(curves.at(currentCluster->indices.at(i)));
    }
    
    // initialize empty vector fields
    vector<pair<Vector*,Vector*> > vectorFields;
    int gridDimension  = g->getResolutionX() * g->getResolutionY();
    for(int i = 0 ; i < numberOfVectorFields ; ++i)
    {
        Vector* xComponent = new Vector(gridDimension);
        Vector* yComponent = new Vector(gridDimension);
        vectorFields.push_back(make_pair(xComponent, yComponent));
    }
    
    // optimize
    op.optimizeImplicitFastWithWeights(
        *g, numberOfVectorFields, curvesInCurrentCluster, vectorFields, &(mapCurveToVF[0]), mapCurveToError,smoothnessWeight
    );
    
    //count number of curves for each vf
    int numberOfChildren = vectorFields.size();
    int numberOfCurvesPerVF[numberOfChildren];
    for(int i = 0 ; i < numberOfChildren ; ++i){
        numberOfCurvesPerVF[i] = 0;
    }

    // initilialize maps
    vector<float> mapVectorFieldToError = vector<float>(vectorFields.size(),0);  // VF -> error
    vector<vector<int> > mapCurvesToClusters(numberOfChildren,vector<int>());  //  
    vector<vector<float> > mapCurveErrorsToClusters(numberOfChildren,vector<float>());
    
    // takes references and populate curves and errors of each cluster
    for(int i = 0 ; i < numberOfCurves ; ++i){
        vector<int>& curveCluster
            = mapCurvesToClusters.at(mapCurveToVF[i]);
        vector<float>& curveErrorsCluster
            = mapCurveErrorsToClusters.at(mapCurveToVF[i]);

        curveCluster.push_back(mapCurveToIndexInCurveVector[i]);
        curveErrorsCluster.push_back(mapCurveToError[i]);

        numberOfCurvesPerVF[mapCurveToVF[i]] += 1;

        float &error = mapVectorFieldToError.at(mapCurveToVF[i]);
        error += mapCurveToError[i];
    }

    // update cluster struct
    /*
    initialize a cluster object for each vector field created
    */
    currentCluster->clearChildren();  // root cluster
    for(int i = 0 ; i < numberOfChildren ; ++i){
        Cluster* c = new Cluster();
        currentCluster->children.push_back(c);

        c->children.clear();
        c->parent = currentCluster;

        stringstream ss;
        ss << currentCluster->name << ":" << i;
        c->name = ss.str();

        c->error = mapVectorFieldToError[i];
        c->indices = mapCurvesToClusters.at(i);
        c->curveErrors = mapCurveErrorsToClusters.at(i);
        c->vectorField = vectorFields.at(i);

        float maxE = -1000;
        for(int j = 0 ; j < (int)c->curveErrors.size() ; ++j)
        {
            float currentError = c->curveErrors.at(j);
            if(currentError > maxE)
                maxE = currentError;
        }
        c->maxError = maxE;
    }

    // save
    saveExperiment(outputDirectory, filename, rootCluster);  // root cluster is currentCluster; his childrens was assigneds

    // clean memory
    cout << "Cleaning Memory" << endl;
    for(int i = 0 ; i < numberOfVectorFields ; ++i){
        std::pair<Vector*, Vector*> &vfs = vectorFields.at(i);
        Vector* componentX  = vfs.first;
        if(componentX != NULL){
            delete componentX;
        }
            
        Vector* componentY  = vfs.second;
        if(componentY != NULL){
            delete componentY;
        }
    }

}
