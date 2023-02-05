/*! \file    CornerGrid.hpp
 *  \brief   Declaration of classes related to the corner grid
 *  \author  Shizhe Li
 *  \date    Nov/16/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __CORNERGRID_HEADER__
#define __CORNERGRID_HEADER__

// Standard header files
#include <math.h>
#include <stdlib.h>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

// Constants used by corner-point grid code
const OCP_DBL SMALL_REAL   = 1E-10; ///< Used for checking determinate of a small matrix
const OCP_DBL TEENY        = 1E-3;  ///< Used for checking distance b/w center to face
const OCP_DBL SMALL        = 1E-3;  ///< Small number as tolerance
const USI     MAX_NEIGHBOR = 80;    ///< Max number of  neighbors allowed

/// A point in 2D.
class Point2D
{
public:
    OCP_DBL x;
    OCP_DBL y;

public:
    Point2D() = default;
    Point2D(OCP_DBL x0, OCP_DBL y0)
        : x(x0)
        , y(y0){};
};

/// A point in 3D.
class Point3D
{
public:
    OCP_DBL x;
    OCP_DBL y;
    OCP_DBL z;

public:
    Point3D() = default;
    Point3D(OCP_DBL x0, OCP_DBL y0, OCP_DBL z0)
        : x(x0)
        , y(y0)
        , z(z0){};

    Point3D& operator=(const Point3D& other);       ///< equal
    Point3D  operator+(const Point3D& other) const; ///< Addition
    Point3D  operator-(const Point3D& other) const; ///< Subtraction
    OCP_DBL  operator*(const Point3D& other) const; ///< Multiplication
    Point3D& operator+=(const Point3D& other);
    Point3D& operator*=(const OCP_DBL& a);
    Point3D& operator/=(const OCP_DBL& a);
    void     Reset()
    {
        x = 0;
        y = 0;
        z = 0;
    };
};

Point3D operator*(const Point3D& p, const OCP_DBL& a);      ///< Point * a
Point3D operator*(const OCP_DBL& a, const Point3D& p);      ///< a * Point
Point3D CrossProduct(const Point3D& p1, const Point3D& p2); ///< Cross product

/// A hexahedron cell.
class Hexahedron
{
public:
    Point3D p0, p1, p2, p3, p4, p5, p6, p7;
};

/// A face of a hexahedron cell.
class HexahedronFace
{
public:
    Point3D p0, p1, p2, p3;
};

/// 3 by 3 matrix.
class Matrix3
{
public:
    OCP_DBL M[3][3];
    Point3D operator*(const Point3D& v) const;
};

/// Get the volume of a hexahedron.
OCP_DBL VolumHexahedron(const Hexahedron& h);

/// Find the center of a hexahedron.
Point3D CenterHexahedron(const Hexahedron& h);

/// Find the normal vector of a face.
Point3D VectorFace(const HexahedronFace& f);

/// Find the center of a face.
Point3D CenterFace(const HexahedronFace& f);

/// ???
Point2D CalCrossingPoint(const Point2D Line1[2], const Point2D Line2[2]);

/// ???
OCP_DBL CalAreaNotQuadr(const HexahedronFace& FACE1, const HexahedronFace& FACE2);

/// ???
class HalfConn
{
public:
    OCP_DBL Ad_dd;
    Point3D d;
    OCP_USI neigh;
    USI     directionType; // 1 - x, 2 - y, 3 - z, 4 - extension
};

/// ???
class ConnGrid
{
public:
    USI              nConn, maxConn;
    vector<HalfConn> halfConn;
    void             Allocate(const USI& max_neighbor);
    void             AddHalfConn(const OCP_USI& n,
                                 const Point3D& area,
                                 const Point3D& d,
                                 const USI&     direction,
                                 const OCP_DBL& flag = 1);
};

/// ???
class GeneralConnect
{
public:
    OCP_USI begin, end;
    USI     directionType; // 1 - x, 2 - y, 3 - z, 4 - extension
    OCP_DBL Ad_dd_begin;
    OCP_DBL Ad_dd_end;
};

/// ???
class OCP_COORD
{
    friend class Grid;

public:
    void     Allocate(const USI& Nx, const USI& Ny, const USI& Nz);
    void     InputData(const vector<OCP_DBL>& coord, const vector<OCP_DBL>& zcorn);
    OCP_BOOL InputCOORDDATA(const vector<OCP_DBL>& coord);
    OCP_BOOL InputZCORNDATA(const vector<OCP_DBL>& zcorn);
    // New version
    void SetupCornerPoints();
    void SetAllFlags(const HexahedronFace& oFace, const HexahedronFace& Face);
    // functions
    OCP_DBL OCP_SIGN(const OCP_DBL& x) { return x >= 0 ? 1 : -1; }

private:
    USI           nx;
    USI           ny;
    USI           nz;
    OCP_DBL***    COORDDATA;
    OCP_DBL****   ZCORNDATA;
    Hexahedron*** cornerPoints;

    OCP_USI         numGrid;
    OCP_USI         numConn;
    OCP_USI         numConnMax;
    vector<OCP_DBL> v;
    vector<OCP_DBL> depth;
    vector<OCP_DBL> dx;
    vector<OCP_DBL> dy;
    vector<OCP_DBL> dz;
    vector<Point3D> center;

    vector<GeneralConnect> connect;

    // Auxiliary variables
    // if the i th point of oFace is deeper than the one of Face, then flagpi = 1;
    // if the i th point of oFace is higher than the one of Face, then flagpi = -1;
    // if the i th point of oFace is very close to the one of Face, then flagpi = 0;
    OCP_INT        flagp0, flagp1, flagp2, flagp3;
    OCP_BOOL       flagQuad;
    OCP_BOOL       upNNC, downNNC;
    OCP_BOOL       flagJump;
    HexahedronFace tmpFace;
    // after the Axes are determined, blocks will be placed along the y+, or along the
    // y- if y+, then flagForward equals 1.0, else -1.0, this relates to calculation of
    // area normal vector
    OCP_DBL flagForward;
};

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/16/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/