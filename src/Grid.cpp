#include "Grid.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#define xDEBUG

using namespace std;


/*
    GRID CLASS
*/


Grid::Grid(float x, float y, float w, float h, int resolutionX, int resolutionY):
    // m_* flags member variable (attribute)
    m_resolutionX(resolutionX),
    m_resolutionY(resolutionY),
    m_x(x),
    m_y(y),
    m_w(w),
    m_h(h)
{
    // escalas (cordenada grid) / (cordenada mundo real)
    m_delta_x = (float)m_w / (float)(m_resolutionX - 1);
    m_delta_y = (float)m_h / (float)(m_resolutionY - 1);
}


TriangularFace Grid::getFaceWherePointLies(const Vector2D &v) const
{
    /*
        REQUIRES POINT IN GRID COORDINATE SYSTEM
    */
    if (  // coordinates out of range
        v.x < 0 || v.x > (m_resolutionX - 1.0) ||
        v.y < 0 || v.y > (m_resolutionY - 1.0)
    )
    {
        cerr << "BAD POINT!" << endl;
        exit(1);
    }
    TriangularFace face;
    
    int square_x = int(v.x), square_y = int(v.y);
    float fx = v.x - square_x, fy = v.y - square_y;
    if (fx > fy) {
        // bottom triangle
        face.indexV1 = vertexIndex(square_x,   square_y  );
        face.indexV2 = vertexIndex(square_x+1, square_y  );
        face.indexV3 = vertexIndex(square_x+1, square_y+1);
    } else {
        // top triangle
        face.indexV1 = vertexIndex(square_x,   square_y  );
        face.indexV2 = vertexIndex(square_x+1, square_y+1);
        face.indexV3 = vertexIndex(square_x,   square_y+1);
    }
    return face;
}

TriangularFace Grid::getFace(int index){
    int line = index / (m_resolutionX - 1);
    int column = index % (m_resolutionX - 1);

    int lowerleftVertex = (line/2) * (m_resolutionX) + column;
    TriangularFace tf;

    if(line % 2 == 0){
        //bottom line
        tf.indexV1 = lowerleftVertex;
        tf.indexV2 = lowerleftVertex + 1;
        tf.indexV3 = lowerleftVertex + m_resolutionX + 1;
    }
    else{
        //top line
        tf.indexV1 = lowerleftVertex;
        tf.indexV2 = lowerleftVertex + m_resolutionX;
        tf.indexV3 = lowerleftVertex + m_resolutionX + 1;
    }

    return tf;
}

TriangularFace Grid::getFace(int index, int resolutionX, int ){
    int line = index / (resolutionX - 1);
    int column = index % (resolutionX - 1);

    int lowerleftVertex = (line/2) * (resolutionX) + column;
    TriangularFace tf;

    if(line % 2 == 0){
        //bottom line
        tf.indexV1 = lowerleftVertex;
        tf.indexV2 = lowerleftVertex + 1;
        tf.indexV3 = lowerleftVertex + resolutionX + 1;
    }
    else{
        //top line
        tf.indexV1 = lowerleftVertex;
        tf.indexV2 = lowerleftVertex + resolutionX;
        tf.indexV3 = lowerleftVertex + resolutionX + 1;
    }

    return tf;
}


///////////////////////////////

void Grid::locate_point(PointLocation &l, const Vector2D& point) const
{
    Vector2D vertex1 = getGridVertex(l.face.indexV1);
    Vector2D vertex2 = getGridVertex(l.face.indexV2);
    Vector2D vertex3 = getGridVertex(l.face.indexV3);

    float det =
        (vertex2.X() - vertex1.X()) * (vertex3.Y() - vertex1.Y()) - 
        (vertex2.Y() - vertex1.Y()) * (vertex3.X() - vertex1.X());

    if(det == 0){
        cout << "det == 0!!!!" << endl;
        exit(1);
    }

    float beta  = ((vertex1.X() - vertex3.X()) * (point.Y() - vertex3.Y()) - (vertex1.Y() - vertex3.Y()) * (point.X() - vertex3.X()))/det;
    float gamma = ((vertex2.X() - vertex1.X()) * (point.Y() - vertex1.Y()) - (vertex2.Y() - vertex1.Y()) * (point.X() - vertex1.X()))/det;
    float alpha = 1.0 - gamma - beta;
    l.barycentric_coords[0] = alpha;
    l.barycentric_coords[1] = beta;
    l.barycentric_coords[2] = gamma;
}

void Grid::multiplyByLaplacian(Vector &firstComponent, Vector &secondComponent) const
{
    if(m_resolutionX * m_resolutionY != firstComponent.getDimension() ||
            m_resolutionX * m_resolutionY != secondComponent.getDimension()){
        cout << "Error while multiplying grid by vector. Incompatible dimensions." << endl;
        exit(1);
    }

    int numberOfVectors = m_resolutionX * m_resolutionY;
    VECTOR_TYPE newFirstComponent[numberOfVectors];
    VECTOR_TYPE newSecondComponent[numberOfVectors];

    //float deltaX = m_w / (m_resolutionX - 1.0);
    //float deltaY = m_h / (m_resolutionY - 1.0);

    float horizontalCotangentWeight = m_delta_x / m_delta_y;
    float verticalCotangentWeight = m_delta_y / m_delta_x;

    //numberOfVectors == numberOfVertices in the grid == resolution * resolution
    for(int i = 0 ; i < numberOfVectors ; ++i){
        int row = i / m_resolutionX;
        int col = i % m_resolutionX;

        bool canMoveLeft = col > 0;
        bool canMoveDown = row > 0;
        bool canMoveRight = col < m_resolutionX - 1;
        bool canMoveUp = row < m_resolutionY - 1;

        float degree = 0;
        float accum1 = 0;
        float accum2 = 0;

        if(canMoveLeft){
            int neighIndex = i - 1;

            //
            float coef = 0.0f;

            if(canMoveUp){
                coef += verticalCotangentWeight;
            }
            if(canMoveDown){
                coef += verticalCotangentWeight;
            }

            coef /= 2.0f;

            //newVector.add((vectorField.at(neighIndex) * coef));
            accum1 += coef * (firstComponent[neighIndex]);
            accum2 += coef * (secondComponent[neighIndex]);

            degree += coef;
        }
        if(canMoveRight){

            int neighIndex = i + 1;

            float coef = 0.0f;
            //

            if(canMoveDown){
                coef += verticalCotangentWeight;
            }
            if(canMoveUp){
                coef += verticalCotangentWeight;
            }

            coef /= 2.0f;

            //newVector.add((vectorField.at(neighIndex) * coef));
            accum1 += coef * (firstComponent[neighIndex]);
            accum2 += coef * (secondComponent[neighIndex]);

            degree += coef;
        }
        if(canMoveDown){

            float coef = 0.0f;

            int neighIndex = i - m_resolutionX;


            if(canMoveLeft){
                coef += horizontalCotangentWeight;
            }
            if(canMoveRight){
                coef += horizontalCotangentWeight;
            }

            coef /= 2.0f;


            //newVector.add((vectorField.at(neighIndex) * coef));

            accum1 += coef * (firstComponent[neighIndex]);
            accum2 += coef * (secondComponent[neighIndex]);

            degree += coef;
        }
        if(canMoveUp){

            int neighIndex = i +  m_resolutionX;

            float coef = 0.0f;

            if(canMoveLeft){
                coef += horizontalCotangentWeight;
            }
            if(canMoveRight){
                coef += horizontalCotangentWeight;
            }

            coef /= 2.0f;

            //newVector.add((vectorField.at(neighIndex) * coef));

            accum1 += coef * (firstComponent[neighIndex]);
            accum2 += coef * (secondComponent[neighIndex]);

            degree += coef;
        }

        //newVector.add(vectorField.at(i) * (-degree));
        newFirstComponent[i] = accum1 - degree * firstComponent[i];
        newSecondComponent[i] = accum2 - degree * secondComponent[i];
    }

    //vectorField.assign(newVectorField.begin(), newVectorField.end());
    firstComponent.setValues(newFirstComponent);    
    secondComponent.setValues(newSecondComponent);    
}

void Grid::multiplyByLaplacian2(Vector &firstComponent) {
    if(m_resolutionX * m_resolutionY != firstComponent.getDimension()) {
        cout << "Error while multiplying grid by vector. Incompatible dimensions." << endl;
        exit(1);
    }

    int numberOfVectors = m_resolutionX * m_resolutionY;
    VECTOR_TYPE newFirstComponent[numberOfVectors];

    float horizontalCotangentWeight = m_delta_x / m_delta_y;
    float verticalCotangentWeight = m_delta_y / m_delta_x;

    //numberOfVectors == numberOfVertices in the grid == resolution * resolution
    for(int i = 0 ; i < numberOfVectors ; ++i){
        int row = i / m_resolutionX;
        int col = i % m_resolutionX;

        bool canMoveLeft = col > 0;
        bool canMoveDown = row > 0;
        bool canMoveRight = col < m_resolutionX - 1;
        bool canMoveUp = row < m_resolutionY - 1;

        float degree = 0;
        float accum1 = 0;

        if(canMoveLeft){
            int neighIndex = i - 1;

            //
            float coef = 0.0;

            if(canMoveUp){
                coef += verticalCotangentWeight;
            }
            if(canMoveDown){
                coef += verticalCotangentWeight;
            }

            coef /= 2.0;

            //newVector.add((vectorField.at(neighIndex) * coef));
            accum1 += coef * (firstComponent[neighIndex]);
            degree += coef;
        }
        if(canMoveRight){

            int neighIndex = i + 1;

            float coef = 0.0;
            //

            if(canMoveDown){
                coef += verticalCotangentWeight;
            }
            if(canMoveUp){
                coef += verticalCotangentWeight;
            }

            coef /= 2.0;

            //newVector.add((vectorField.at(neighIndex) * coef));
            accum1 += coef * (firstComponent[neighIndex]);
            degree += coef;
        }
        if(canMoveDown){

            float coef = 0.0;

            int neighIndex = i - m_resolutionX;


            if(canMoveLeft){
                coef += horizontalCotangentWeight;
            }
            if(canMoveRight){
                coef += horizontalCotangentWeight;
            }

            coef /= 2.0;


            //newVector.add((vectorField.at(neighIndex) * coef));

            accum1 += coef * (firstComponent[neighIndex]);
            degree += coef;
        }
        if(canMoveUp){

            int neighIndex = i +  m_resolutionX;

            float coef = 0.0;

            if(canMoveLeft){
                coef += horizontalCotangentWeight;
            }
            if(canMoveRight){
                coef + = horizontalCotangentWeight;
            }

            coef /= 2.0;

            //newVector.add((vectorField.at(neighIndex) * coef));

            accum1 += coef * (firstComponent[neighIndex]);
            degree += coef;
        }

        //newVector.add(vectorField.at(i) * (-degree));
        newFirstComponent[i] = accum1 - degree * firstComponent[i];
    }

    //vectorField.assign(newVectorField.begin(), newVectorField.end());
    firstComponent.setValues(newFirstComponent);    
}

//#define DEBUG

vector<Grid::Inter> Grid::clipAgainstHorizontalLines(
    const Grid::Inter &g1, const Grid::Inter &g2
) const
{
    /*
        take endpoints of segment,
        return vector of its intersections with horizontal lines
        (for each horizontal lines between the endpoints,
        calculate intersection using linear interpolation)
    */

    vector<Grid::Inter> result;  // result is a vector of intersections

    if (g1.grid_point.Y() == g2.grid_point.Y()) {  // segment is perfectly horizontal
        return result;  // return empty list  // parallel lines do not cross
    }

    float inv_slope =  // inverse of slope = run/rise  // so never is dividing by 0
        (g2.grid_point.X() - g1.grid_point.X()) / (g2.grid_point.Y() - g1.grid_point.Y());

        
    // now that has the slope, can calculate the intersection with each horizontal line using linear interpolation

    if (inv_slope >= 0) 
    {
        float y1 = g1.grid_point.Y(), y2 = g2.grid_point.Y();

        float this_y = y1;
        float next_y = floor(this_y) + 1;

        while (next_y < y2) {
            float u = (next_y - y1) / (y2 - y1);
            Grid::Inter p;
            p.grid_point = Vector2D::lerp(g1.grid_point, g2.grid_point, u);
            // Force correct (integer) coordinate prevent trouble down the line with floating point accuracy.
            p.grid_point.y = next_y;
            p.u = u;
            p.kind = Grid::Inter::Horizontal;
            result.push_back(p);
            this_y = next_y;
            next_y = this_y + 1;
        }
        return result;
    } 

    else  // inv_slope < 0
    {
        /*
        If inv_slope is negative, we transform the problem by flipping the grid upside down.
        Now inv_slope will be positive. 
        Then transform back.
        */

        // flip the grid
        Grid::Inter reverse_g1, reverse_g2;
        float t = (m_resolutionY - 1);
        reverse_g1.grid_point = Vector2D(g1.grid_point.x, t - g1.grid_point.y);
        reverse_g2.grid_point = Vector2D(g2.grid_point.x, t - g2.grid_point.y);
        reverse_g1.u = g1.u;
        reverse_g2.u = g2.u;
        // compute flipped result
        vector<Grid::Inter> reverse_result = clipAgainstHorizontalLines(reverse_g1, reverse_g2);  // call method again with reversed input
        // unflip
        for (size_t i=0; i<reverse_result.size(); ++i) {
            reverse_result[i].grid_point.y = t - reverse_result[i].grid_point.y;
        }
        
        return reverse_result;
    }
}

vector<Grid::Inter> Grid::clipAgainstVerticalLines
    (const Grid::Inter &g1, const Grid::Inter &g2) const
{
    // Flip coordinates so vertical lines become horizontal. Solve the problem,
    // transform back.
    Grid::Inter flipped_g1, flipped_g2;
    flipped_g1.grid_point = Vector2D(g1.grid_point.y, g1.grid_point.x);
    flipped_g2.grid_point = Vector2D(g2.grid_point.y, g2.grid_point.x);
    flipped_g1.u = g1.u;
    flipped_g2.u = g2.u;
    vector<Grid::Inter> flipped_result =
        clipAgainstHorizontalLines(flipped_g1, flipped_g2);
    for (vector<Grid::Inter>::iterator b = flipped_result.begin(),
             e = flipped_result.end(); b!=e; ++b) {
        b->grid_point = Vector2D(b->grid_point.y, b->grid_point.x);
        b->kind = Grid::Inter::Vertical;
    }
    return flipped_result;
}

inline static float get_u_from_points(const Vector2D &v1, const Vector2D &v2,
                                      const Vector2D &u)
{
    if (v2.x != v1.x) {
        return (u.x - v1.x) / (v2.x - v1.x);
    } else {
        return (u.y - v1.y) / (v2.y - v1.y);
    }
}

// http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T>
static int sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

void Grid::clipLine(PolygonalPath &path1) const  // tesselation - find all points where path intersects 
{
    /*
        insert new vertices wherever the path intersects with a grid line

        divide path1 into his segments
        then clipAgainstHorizontalLines and clipAgainstVerticalLines do the actual tesselation
    */
   
    // get segments
    size_t currentVertexIndex = 0;
    while(currentVertexIndex < path1.numberOfPoints() - 1){  // for each segment:
        // position of start and end
        Vector2D from = path1.getPoint(currentVertexIndex).first,  
                   to = path1.getPoint(currentVertexIndex+1).first;
        // time of start and end
        float   tfrom = path1.getPoint(currentVertexIndex).second, 
                  tto = path1.getPoint(currentVertexIndex+1).second;
        // direction
        Vector2D tangent = path1.getTangent(currentVertexIndex);

        // convert cordinates World->Grid
        vector<Grid::Inter> inters;
        Grid::Inter e1, e2;
        e1.grid_point = toGrid(from);
        e2.grid_point = toGrid(to);
        e1.u = 0;  // 0 = start 
        e2.u = 1;  // 1 = end
        e1.kind = e2.kind = Grid::Inter::EndPoint;

        inters.push_back(e1);

        vector<Grid::Inter> horiz = clipAgainstHorizontalLines(e1, e2);
        horiz.push_back(e2);  
        for (size_t i=0; i<horiz.size(); ++i) {  // for each horizontal intersection
            vector<Grid::Inter> vert = clipAgainstVerticalLines(inters.back(), horiz[i]);
            // reported barycentric coords are with respect to clipped lines;
            // make a good u here.
            for (size_t j=0; j<vert.size(); ++j) {
                vert[j].u = get_u_from_points(e1.grid_point, e2.grid_point, vert[j].grid_point);
            }
            copy(vert.begin(), vert.end(), back_inserter(inters));
            inters.push_back(horiz[i]);
        }
        
        // resolve diagonal intersections which remain
        for (size_t i=0; i<inters.size()-1; ++i) {
            int x_square = min(int(inters[i  ].grid_point.x),
                               int(inters[i+1].grid_point.x)),
                y_square = min(int(inters[i  ].grid_point.y),
                               int(inters[i+1].grid_point.y));
            float
                u1 = inters[i  ].grid_point.x - x_square,
                v1 = inters[i  ].grid_point.y - y_square,
                u2 = inters[i+1].grid_point.x - x_square,
                v2 = inters[i+1].grid_point.y - y_square;
            float du = u2 - u1, dv = v2 - v1;
            int s1 = sgn(u1 - v1);
            int s2 = sgn(u2 - v2);

            // if sign of x_i - y_i is different than that of
            // x_i+1 - y_i+1 then there's an intersection with the diagonal.

            if (s1 != s2) {
                // If that's the case, solve the following linear system over point-in-square 
                // coordinates:

                // x = y (diagonal line)
                // (y-y2)/(y2-y1) = (x-x2)/(x2-x1)

                //   or equivalently without divisions:
                // (y-y2)dx = (x-x2)dy

                // The solution is x = y = (v2 du - u2 dv) / (du - dv)
                
                // if du = dv we wouldn't have arrived here, because their signs would
                // have been different
                float x = (v2 * du - u2 * dv) / (du - dv);
                float y = x;
                Grid::Inter new_inter;
                new_inter.grid_point = Vector2D(x_square + x, y_square + y);
                new_inter.u = get_u_from_points(e1.grid_point, e2.grid_point, new_inter.grid_point);
                new_inter.kind = Grid::Inter::Diagonal;
                inters.insert(inters.begin() + (i + 1), new_inter);

                // Since we just inserted a point, and we don't want
                // to check diagonal intersections against it, we increment the index variable.
                ++i;
            }

            // At the end of all this, we insert points 1..end-1, which correspond to
            // new intersections. inters[0] is the origin vertex.
        }
        for (size_t i=1; i<inters.size()-1; ++i) {
            path1.addPoint(++currentVertexIndex,
                           make_pair(toWorld(inters[i].grid_point),
                                     inters[i].u * tto + (1 - inters[i].u) * tfrom),
                           tangent);
        }
        currentVertexIndex++;
    }
}

void update_triplet_matrix
(int index,
 Intersection &cj,
 Intersection &cj_plus_1,
 vector<int> &row_vector,
 vector<int> &col_vector,
 vector<float> &values)
{
    int r1 = index, r2 = index+1;
    
    float
        u1 = 1.0-cj.lambda,        v1 = cj.lambda,
        u2 = 1.0-cj_plus_1.lambda, v2 = cj_plus_1.lambda;

#define ADD_VALUE(r, c, v)                      \
    row_vector.push_back(r);                    \
    col_vector.push_back(c);                    \
    values.push_back(v);

    ADD_VALUE(r1, cj.indexV1, u1 / 3.0);
    ADD_VALUE(r1, cj.indexV2, v1 / 3.0);
    ADD_VALUE(r2, cj.indexV1, u1 / 6.0);
    ADD_VALUE(r2, cj.indexV2, u1 / 6.0);

    ADD_VALUE(r1, cj_plus_1.indexV1, u2 / 6.0);
    ADD_VALUE(r1, cj_plus_1.indexV2, v2 / 6.0);
    ADD_VALUE(r2, cj_plus_1.indexV1, u2 / 3.0);
    ADD_VALUE(r2, cj_plus_1.indexV2, u2 / 3.0);
#undef ADD_VALUE
}

CurveDescription Grid::curve_description(const PolygonalPath &path) const
{
    CurveDescription result;
    Vector &right_hand_side_x = result.rhsx;
    Vector &right_hand_side_y = result.rhsy;
    result.length = 0;
    int numberOfPoints = path.numberOfPoints();
    if(numberOfPoints < 2){
        return result;
    }

    right_hand_side_x = Vector(2*(numberOfPoints-1));
    right_hand_side_y = Vector(2*(numberOfPoints-1));

    for(int i = 0 ; i < numberOfPoints - 1 ; ++i){
        Vector2D current = toGrid(path.getPoint(i).first);
        Vector2D next = toGrid(path.getPoint(i+1).first);
        Vector2D midpoint = (current + next) * 0.5;
        Vector2D desired_tangent = path.getTangent(i);
        PointLocation location;
	//initialize to avoid compilation warnings
	location.barycentric_coords[0] = -1;
	location.barycentric_coords[1] = -1;
	location.barycentric_coords[2] = -1;
        //
	location.face = getFaceWherePointLies(midpoint);
        Segment segment;
        segment.endpoint[0] = location;
        segment.endpoint[1] = location;
        locate_point(segment.endpoint[0], current);
        locate_point(segment.endpoint[1], next);
        segment.time[0] = path.getPoint(i).second;
        segment.time[1] = path.getPoint(i+1).second;
        segment.index = 2 * i;
        result.length += (segment.time[1] - segment.time[0]);
        result.segments.push_back(segment);
        right_hand_side_x[2*i]   = desired_tangent.X();
        right_hand_side_x[2*i+1] = desired_tangent.X();
        right_hand_side_y[2*i]   = desired_tangent.Y();
        right_hand_side_y[2*i+1] = desired_tangent.Y();
    }
    return result;
}