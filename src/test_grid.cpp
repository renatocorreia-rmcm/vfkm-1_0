#include <iostream>
#include "Grid.h" 
#include "Vector.h"  
#include "Vector2D.h"

using namespace std;


/*
g++ -o test_grid test_grid.cpp Grid.cpp Vector.cpp Vector2D.cpp
*/



int main() {

    std::cout << "Starting Grid test\n\n";

    Grid g = Grid(0, 0, 1, 1, 3, 3);  // grid arguments when running `./vfkm ../data/trajectories.txt 3 2 0.5 ../resultclear`

    cout << "grid parameters: \n";
    cout << g.m_resolutionX << " " << g.m_resolutionY << " " << g.m_x << " " << g.m_y << " " << g.m_w << " " << g.m_h << " " << g.m_delta_x << " " << g.m_delta_y << '\n';



    return 0;
}