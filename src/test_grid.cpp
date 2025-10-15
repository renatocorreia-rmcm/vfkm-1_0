#include <iostream>
#include "Grid.h" 
#include "Vector.h"  
#include "Vector2D.h"
#include "PolygonalPath.h"

using namespace std;


/*
g++ -o test_grid test_grid.cpp Grid.cpp Vector.cpp Vector2D.cpp
*/



int main() {

    std::cout << "Starting Grid test\n\n";

    Grid g = Grid(0, 0, 1, 1, 3, 3);  // grid arguments when running `./vfkm ../data/trajectories.txt 3 2 0.5 ../resultclear`

    cout << "grid parameters: \n";
    cout << g.m_resolutionX << " " << g.m_resolutionY << " " << g.m_x << " " << g.m_y << " " << g.m_w << " " << g.m_h << " " << g.m_delta_x << " " << g.m_delta_y << '\n';

    std::vector<std::pair<Vector2D,float> > points = {
        {{0.05, 0.6}, 0},
	    {{0.45, 0.65}, 0.25},
	    {{0.75, 0.55}, 0.5},
	    {{0.95, 0.6}, 1}
    };

    PolygonalPath path = PolygonalPath(points);

    cout << "\n path before tesselation\n";


    for (int i = 0; i < path.points.size(); i++)
    {
        std::pair<Vector2D, float> point = path.getPoint(i);
        Vector2D pos = point.first;
        float time = point.second;

        cout << "( (" << pos.x << ", " << pos.y << "), " << time << " ), ";
    }

    g.clipLine(path);


    cout << "\n path after tesselation\n";

    for (int i = 0; i < path.points.size(); i++)
    {
        std::pair<Vector2D, float> point = path.getPoint(i);
        Vector2D pos = point.first;
        float time = point.second;

        cout << "( (" << pos.x << ", " << pos.y << "), " << time << " ), ";
    }


    return 0;
}