/****************************** coordinate_system.hpp *****************************

 * Classes for convenient representation of 3D coordinate systems.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 20, 2018

 **********************************************************************************/


#ifndef COORD_H
#define COORD_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <complex>
#include <vector>
#include <map>
#include <stdlib.h>
#include "math.h"


// Useful reference: http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html


/*! \brief Class to store and manipulate points or position vectors in 3D, in Cartesian coordinates.*/
class vect3D
{
public:

	double x = 0.0, y = 0.0, z = 0.0;
	double tol = 1.0e-15;

	// Constructors
	vect3D() {};
    vect3D(std::initializer_list<double> RHS) { x = *RHS.begin(); y = *(RHS.begin()+1); z = *(RHS.begin()+2); };
	vect3D(std::vector<double> RHS) { x = RHS[0]; y = RHS[1]; z = RHS[2]; };
	vect3D(double *RHS) { x = RHS[0]; y = RHS[1]; z = RHS[2]; };

	void SetCustomTolerance(double _tol) { tol = _tol; };

    // Assignment
    vect3D& operator=(const vect3D &RHS) { x = RHS.x; y = RHS.y; z = RHS.z; return *this; };

    // Addition and subtraction
    vect3D& operator+=(const vect3D &RHS) { x += RHS.x; y += RHS.y; z += RHS.z; return *this; };
    vect3D& operator-=(const vect3D &RHS) { x -= RHS.x; y -= RHS.y; z -= RHS.z; return *this; };
    
    vect3D& operator+=(const double &RHS) { x += RHS; y += RHS; z += RHS; return *this; };
    vect3D& operator-=(const double &RHS) { x -= RHS; y -= RHS; z -= RHS; return *this; };

    vect3D operator+(const vect3D &RHS) { return vect3D(*this) += RHS; };
    vect3D operator-(const vect3D &RHS) { return vect3D(*this) -= RHS; };

    vect3D operator+(const double &RHS) { return vect3D(*this) += RHS; };
    vect3D operator-(const double &RHS) { return vect3D(*this) -= RHS; };

    // Scalar product and division
    vect3D& operator*=(const double &RHS) { x *= RHS; y *= RHS; z *= RHS; return *this; };
    vect3D& operator/=(const double &RHS) { x /= RHS; y /= RHS; z /= RHS; return *this; };

    vect3D& operator*(const double &RHS) { return vect3D(*this) *= RHS; };
    vect3D& operator/(const double &RHS) { return vect3D(*this) /= RHS; };

    // Comparisons
    bool operator==(const vect3D &RHS) { return ((x - RHS.x < x*tol) && (y - RHS.y < y*tol) && (z - RHS.z < z*tol)); };
    bool operator!=(const vect3D &RHS) { return !(*this == RHS); };
    
    bool operator>=(const vect3D &RHS) { return ((x >= RHS.x) && (y >= RHS.y) && (z >= RHS.z)); };
    bool operator<=(const vect3D &RHS) { return ((x <= RHS.x) && (y <= RHS.y) && (z <= RHS.z)); };
    bool operator>(const vect3D &RHS) { return !(*this <= RHS); };
    bool operator<(const vect3D &RHS) { return !(*this >= RHS); };


    // Euclidean norm
    double norm2() { return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)); };

    // Dot product
    friend double dot(vect3D const &a, vect3D const &b) { return a.x*b.x + a.y*b.y + a.z*b.z; };

    // Cross product
    friend vect3D cross(vect3D const &a, vect3D const &b) { return {(a.y*b.z - a.z*b.y), (a.z*b.x - a.x*b.z), (a.x*b.y - a.y*b.x)}; };

    // Print to stream
    friend std::ostream& operator<<(std::ostream &stream, vect3D const &p) { return stream << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::flush; };
    

};


#endif
