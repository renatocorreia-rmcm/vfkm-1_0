#include "Util.h"
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <queue>
#include "Vector2D.h"

#define xDEBUG

using namespace std;



/*
    CLUSTER CLASS
*/



void Cluster::clearChildren(){

    // initialize toProcess (queue of clusters) with children (vector of clusters)
    queue<Cluster*> toProcess;
    int numberOfChildren = children.size();
    for(int i = 0 ; i < numberOfChildren ; ++i) toProcess.push(children.at(i));

    // delete each cluster of it
    while(!toProcess.empty()){
        // load cluster
        Cluster* cluster = toProcess.front();
        toProcess.pop();
        // delete its VF
        delete cluster->vectorField.first;
        delete cluster->vectorField.second;
        // add its children to toProcess queue
        int numberOfChildren = cluster->children.size();
        for(int i = 0 ; i < numberOfChildren ; ++i){
            toProcess.push(cluster->children.at(i));
        }
        // delete it
        delete cluster;
    }

    children.clear();
}



/*
    UTIL STATIC CLASS
*/

void Util::loadCurves(
    string filename, vector<PolygonalPath>& curves,
    float &xmin, float &xmax, float &ymin, float &ymax, float &tmin, float &tmax  // also loads bounding box of array
)
{
    /*
    read bounding box and store in given arguments.
    read raw points {space, time} and store they in curves array

    no tesselation yet,
    but removes consecutive time/space equal points 
    */

    cout << "loading curves" << endl;

    xmin = +FLT_MAX;
    ymin = +FLT_MAX;
    tmin = +FLT_MAX;
    xmax = -FLT_MAX;
    ymax = -FLT_MAX;
    tmax = -FLT_MAX;

    // int stride = 0;
    ifstream file (filename.c_str());
    int real_index = 0;
    
    if (file.is_open())
    {
      file >> xmin >> xmax >> ymin >> ymax >> tmin >> tmax;

        vector<pair<Vector2D,float> > curveContents;

        while (file.good()  && !file.eof()) 
        {
            /*
                pick lines (coordinates {space, time}) until end of file
            */

            float x,y,t;
            file >> x >> y >> t;

            if(x == 0 && y == 0 && t == 0)  // finished curve
            {  
                if(curveContents.size() >=2)  // if curve has minimum size
                {  

                    curves.push_back(PolygonalPath(curveContents));  // append new polygonal path to curves array
                }
                real_index++;
                curveContents.clear();  // clear curve container to store the next one
            } 
            else if (  // out of Bounding Box - end of curve
                x < xmin || x > xmax ||
                y < ymin || y > ymax ||
                t < tmin || t > tmax
            ) 
            {  
                if(curveContents.size() >=2){
                    curves.push_back(PolygonalPath(curveContents));
                }
                // WHY NOT real_index++ IN THIS CASE ??
                curveContents.clear();
            } 
            else  // append new point {space, time} to curve
            {
                pair<Vector2D, float> newPoint = {Vector2D(x,y), t};

                if (curveContents.size() == 0)  // is the first point
                    curveContents.push_back(newPoint);  // regular push

                else if (t == curveContents.back().second)  // repeated last timestamp
                    continue;  //  ignore this point
                else if (x == curveContents.back().first.X() && y == curveContents.back().first.Y())  // repeated last position
                    continue;  // ignore this point

                else  // regular point in trajectory
                    curveContents.push_back(newPoint);  // regular push
            }
        }

        file.close();
    }
    else
    {
        cerr << "Unable to open file " << filename << endl;
    }


    // output readed data

    int numberOfCurvesRead = curves.size();
    ofstream outfile("read_curves.txt");

    outfile << xmin << " " << xmax << " " << ymin << " " << ymax << " " << tmin << " " << tmax;  // bounding box

    for(int i = 0 ; i < numberOfCurvesRead ; ++i){
        //cout << "Curve " << i << " = " << curves.at(i).toString() << endl;
        PolygonalPath &curveContents = curves.at(i);
        int curveSize = curveContents.numberOfPoints();

        for(int j = 0 ; j < curveSize ; ++j){
            pair<Vector2D,float> pointTime = curveContents.getPoint(j);
            outfile << pointTime.first.X() << " " << pointTime.first.Y() << " " << pointTime.second << endl;
        }
        outfile << "0 0 0" << endl;
    }





    #ifdef DEBUG
        int numberOfCurvesRead = curves.size();

        cout << "numberOfCurvesRead = " << numberOfCurvesRead << endl;

        for(int i = 0 ; i < numberOfCurvesRead ; ++i){
            cout << "Curve " << i << " = " << curves.at(i).toString() << endl;
        }
    #endif   
}

void Util::loadCurves(std::string filename, std::vector<PolygonalPath>& curves){
    /*
        runs loadCurves without loading bounding box (See loadCurves)
    */
    float xmin, xmax, ymin, ymax, tmin, tmax;
    loadCurves(filename, curves, xmin,xmax,ymin,ymax,tmin,tmax);
}
