#include "Vector.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdlib>

using namespace std;

Vector::Vector(int dimension):
    dimension(dimension)
{
    values = new VECTOR_TYPE[dimension];
    if(values == NULL){
        cout << "Error while allocating Vector" << endl;
        exit(1);
    }
    setValues(0.0);
}

std::ostream& operator<<(std::ostream& os, const Vector& obj) {
    os << obj.toString();
    return os;
}



Vector::Vector(const Vector& v){
    dimension = v.getDimension();
    values = new VECTOR_TYPE[dimension];
    if(values == NULL){
        cout << "Error while allocating Vector" << endl;
        exit(1);
    }

    setValues(v);
}

Vector::~Vector(){
    if(values != NULL)
        delete[] values;
}

string Vector::toString() const{
    stringstream ss;
    ss << "[ ";
    for(int i = 0 ; i < dimension - 1 ; ++i){

        ss << values[i] << ", ";
        
        if ((i+1)%4==0) ss << "\n";  // keep indentation as python to ease comparission
        
    }
    ss << values[dimension - 1] << "";

    return ss.str();
}
