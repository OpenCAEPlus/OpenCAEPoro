/*! \file    CornerGrid.cpp
 *  \brief   Declaration of classes related to the corner grid
 *  \author  Shizhe Li
 *  \date    Nov/19/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "CornerGrid.hpp"

Point3D& Point3D::operator=(const Point3D& other)
{
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

Point3D Point3D::operator+(const Point3D& other) const
{
    return Point3D(x + other.x, y + other.y, z + other.z);
}

Point3D Point3D::operator-(const Point3D& other) const
{
    return Point3D(x - other.x, y - other.y, z - other.z);
}

OCP_DBL Point3D::operator*(const Point3D& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

Point3D& Point3D::operator+=(const Point3D& other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Point3D& Point3D::operator*=(const OCP_DBL& a)
{
    x *= a;
    y *= a;
    z *= a;
    return *this;
}

Point3D& Point3D::operator/=(const OCP_DBL& a)
{
    x /= a;
    y /= a;
    z /= a;
    return *this;
}

Point3D operator*(const Point3D& p, const OCP_DBL& a)
{
    return Point3D(a * p.x, a * p.y, a * p.z);
}

Point3D operator*(const OCP_DBL& a, const Point3D& p)
{
    return Point3D(a * p.x, a * p.y, a * p.z);
}

Point3D CrossProduct(const Point3D& p1, const Point3D& p2)
{
    Point3D result;
    result.x = p1.y * p2.z - p1.z * p2.y;
    result.y = p1.z * p2.x - p1.x * p2.z;
    result.z = p1.x * p2.y - p1.y * p2.x;
    return result;
}

Point3D Matrix3::operator*(const Point3D& v) const
{
    Point3D result;
    result.x = M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z;
    result.y = M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z;
    result.z = M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z;
    return result;
}

OCP_DBL VolumHexahedron(const Hexahedron& h)
{
    OCP_DBL result =
        (h.p0.x * (h.p1.y * (-h.p2.z - h.p3.z + h.p4.z + h.p5.z) +
                   h.p2.y * (h.p1.z - h.p3.z) +
                   h.p3.y * (h.p1.z + h.p2.z - h.p4.z - h.p7.z) +
                   h.p4.y * (-h.p1.z + h.p3.z - h.p5.z + h.p7.z) +
                   h.p5.y * (-h.p1.z + h.p4.z) + h.p7.y * (h.p3.z - h.p4.z)) +
         h.p1.x * (h.p0.y * (+h.p2.z + h.p3.z - h.p4.z - h.p5.z) +
                   h.p2.y * (-h.p0.z - h.p3.z + h.p5.z + h.p6.z) +
                   h.p3.y * (-h.p0.z + h.p2.z) + h.p4.y * (h.p0.z - h.p5.z) +
                   h.p5.y * (h.p0.z - h.p2.z + h.p4.z - h.p6.z) +
                   h.p6.y * (-h.p2.z + h.p5.z)) +
         h.p2.x * (h.p0.y * (-h.p1.z + h.p3.z) +
                   h.p1.y * (h.p0.z + h.p3.z - h.p5.z - h.p6.z) +
                   h.p3.y * (-h.p0.z - h.p1.z + h.p6.z + h.p7.z) +
                   h.p5.y * (h.p1.z - h.p6.z) +
                   h.p6.y * (h.p1.z - h.p3.z + h.p5.z - h.p7.z) +
                   h.p7.y * (-h.p3.z + h.p6.z)) +
         h.p3.x * (h.p0.y * (-h.p1.z - h.p2.z + h.p4.z + h.p7.z) +
                   h.p1.y * (h.p0.z - h.p2.z) +
                   h.p2.y * (h.p0.z + h.p1.z - h.p6.z - h.p7.z) +
                   h.p4.y * (-h.p0.z + h.p7.z) + h.p6.y * (h.p2.z - h.p7.z) +
                   h.p7.y * (-h.p0.z + h.p2.z - h.p4.z + h.p6.z)) +
         h.p4.x * (h.p0.y * (h.p1.z - h.p3.z + h.p5.z - h.p7.z) +
                   h.p1.y * (-h.p0.z + h.p5.z) + h.p3.y * (h.p0.z - h.p7.z) +
                   h.p5.y * (-h.p0.z - h.p1.z + h.p6.z + h.p7.z) +
                   h.p6.y * (-h.p5.z + h.p7.z) +
                   h.p7.y * (h.p0.z + h.p3.z - h.p5.z - h.p6.z)) +
         h.p5.x * (h.p0.y * (h.p1.z - h.p4.z) +
                   h.p1.y * (-h.p0.z + h.p2.z - h.p4.z + h.p6.z) +
                   h.p2.y * (-h.p1.z + h.p6.z) +
                   h.p4.y * (h.p0.z + h.p1.z - h.p6.z - h.p7.z) +
                   h.p6.y * (-h.p1.z - h.p2.z + h.p4.z + h.p7.z) +
                   h.p7.y * (h.p4.z - h.p6.z)) +
         h.p6.x * (h.p1.y * (h.p2.z - h.p5.z) +
                   h.p2.y * (-h.p1.z + h.p3.z - h.p5.z + h.p7.z) +
                   h.p3.y * (-h.p2.z + h.p7.z) + h.p4.y * (h.p5.z - h.p7.z) +
                   h.p5.y * (h.p1.z + h.p2.z - h.p4.z - h.p7.z) +
                   h.p7.y * (-h.p2.z - h.p3.z + h.p4.z + h.p5.z)) +
         h.p7.x * (h.p0.y * (-h.p3.z + h.p4.z) + h.p2.y * (h.p3.z - h.p6.z) +
                   h.p3.y * (h.p0.z - h.p2.z + h.p4.z - h.p6.z) +
                   h.p4.y * (-h.p0.z - h.p3.z + h.p5.z + h.p6.z) +
                   h.p5.y * (-h.p4.z + h.p6.z) +
                   h.p6.y * (h.p2.z + h.p3.z - h.p4.z - h.p5.z))) /
        12;
    return result;
}

Point3D CenterHexahedron(const Hexahedron& h)
{
    OCP_DBL r      = 1.0 / 8.0;
    Point3D result = r * (h.p0 + h.p1 + h.p2 + h.p3 + h.p4 + h.p5 + h.p6 + h.p7);
    return result;
}

Point3D VectorFace(const HexahedronFace& f)
{
    Point3D p0, p1, p2, p3;
    p0             = f.p2 - f.p1;
    p1             = f.p0 - f.p1;
    p2             = f.p0 - f.p3;
    p3             = f.p2 - f.p3;
    Point3D result = 0.5 * (CrossProduct(p0, p1) + CrossProduct(p2, p3));
    return result;
}

Point3D CenterFace(const HexahedronFace& f)
{
    OCP_DBL r      = 1.0 / 4.0;
    Point3D result = r * (f.p0 + f.p1 + f.p2 + f.p3);
    return result;
}

Point2D CalCrossingPoint(const Point2D Line1[2], const Point2D Line2[2])
{
    Point2D crosspoint;
    //
    //   LOCALS
    //
    OCP_DBL a11, a12, a21, a22, b1, b2, detA, detX, detY;
    //
    //   assume   x   =   crosspoint.x
    //            y   =   crosspoint.y
    //   calculate x and y with equations in the following
    //
    //    a11 a12     x       b1
    //   [        ] (   ) = (    )
    //    a21 a22     y       b2
    //
    a11  = Line1[1].y - Line1[0].y;
    a12  = Line1[0].x - Line1[1].x;
    a21  = Line2[1].y - Line2[0].y;
    a22  = Line2[0].x - Line2[1].x;
    b1   = a11 * Line1[0].x + a12 * Line1[0].y;
    b2   = a21 * Line2[0].x + a22 * Line2[0].y;
    detA = a11 * a22 - a12 * a21;

    if (fabs(detA) > SMALL_REAL) {
        detX         = b1 * a22 - b2 * a12;
        detY         = a11 * b2 - a21 * b1;
        crosspoint.x = detX / detA;
        crosspoint.y = detY / detA;
    } else {
        crosspoint = Line1[0];
    }
    return crosspoint;
}

OCP_DBL CalAreaNotQuadr(const HexahedronFace& FACE1, const HexahedronFace& FACE2)
{
    // Attention! Only for non quadrilateral!!!  ---- Lishizhe
    //
    // This function calculate the common area of two quadrilaterals FACE1, FACE2.
    //
    // Order of points of Face follows
    //       1 --- 0        0 --- 1
    //       |     |    or  |     |
    //       2 --- 3        3 --- 2
    // p0, p1 are upper, p2, p3 are lower
    // y must be depth!!!
    //
    OCP_DBL CalAreaNotQuadr;
    //
    //   LOCALS
    //
    USI            iret;
    Point2D        crosspoint[4];
    Point2D        Line1[2], Line2[2];
    HexahedronFace FACEtmp1, FACEtmp2;
    Point3D        area, point1, point2, point3;
    //
    CalAreaNotQuadr = 0;
    iret            = 0;
    //
    //   the crossing relations of 4 lines:
    //           Line1 : point0 and point1 of face1
    //           Line2 : point2 and point3 of face1
    //           Line3 : point0 and point1 of face2
    //           Line4 : point2 and point3 of face2
    //
    //   Line1 & Line3
    //
    Line1[0]      = Point2D(FACE1.p0.x, FACE1.p0.y);
    Line1[1]      = Point2D(FACE1.p1.x, FACE1.p1.y);
    Line2[0]      = Point2D(FACE2.p0.x, FACE2.p0.y);
    Line2[1]      = Point2D(FACE2.p1.x, FACE2.p1.y);
    crosspoint[0] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[0].x - Line1[0].x) * (crosspoint[0].x - Line1[1].x) < 0)
        iret = iret + 1;
    //
    //   Line2 & Line3
    //
    Line1[0]      = Point2D(FACE1.p2.x, FACE1.p2.y);
    Line1[1]      = Point2D(FACE1.p3.x, FACE1.p3.y);
    Line2[0]      = Point2D(FACE2.p0.x, FACE2.p0.y);
    Line2[1]      = Point2D(FACE2.p1.x, FACE2.p1.y);
    crosspoint[1] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[1].x - Line1[0].x) * (crosspoint[1].x - Line1[1].x) < 0)
        iret = iret + 2;
    //
    //   Line1 & Line4
    //
    Line1[0]      = Point2D(FACE1.p0.x, FACE1.p0.y);
    Line1[1]      = Point2D(FACE1.p1.x, FACE1.p1.y);
    Line2[0]      = Point2D(FACE2.p2.x, FACE2.p2.y);
    Line2[1]      = Point2D(FACE2.p3.x, FACE2.p3.y);
    crosspoint[2] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[2].x - Line1[0].x) * (crosspoint[2].x - Line1[1].x) < 0)
        iret = iret + 4;
    //
    //   Line2 & Line4
    //
    Line1[0]      = Point2D(FACE1.p2.x, FACE1.p2.y);
    Line1[1]      = Point2D(FACE1.p3.x, FACE1.p3.y);
    Line2[0]      = Point2D(FACE2.p2.x, FACE2.p2.y);
    Line2[1]      = Point2D(FACE2.p3.x, FACE2.p3.y);
    crosspoint[3] = CalCrossingPoint(Line1, Line2);
    if ((crosspoint[3].x - Line1[0].x) * (crosspoint[3].x - Line1[1].x) < 0)
        iret = iret + 8;
    //
    //   consider 12 cases of crossing relation combinations
    //
    switch (iret) {
        case 1:
            //
            //  Line1 & Line3 only
            //
            FACEtmp1.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp2.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);

            if (FACE1.p0.y > FACE2.p0.y) {
                FACEtmp1.p0 = FACE1.p0;
            } else {
                FACEtmp1.p0 = FACE2.p0;
            }

            if (FACE1.p1.y > FACE2.p1.y) {
                FACEtmp2.p1 = FACE1.p1;
            } else {
                FACEtmp2.p1 = FACE2.p1;
            }

            if (FACE1.p3.y > FACE2.p3.y) {
                FACEtmp1.p3 = FACE2.p3;
            } else {
                FACEtmp1.p3 = FACE1.p3;
            }

            if (FACE1.p2.y > FACE2.p2.y) {
                FACEtmp2.p2 = FACE2.p2;
            } else {
                FACEtmp2.p2 = FACE1.p2;
            }

            FACEtmp1.p2 = Point3D(0.5 * (FACEtmp1.p3.x + FACEtmp2.p2.x),
                                  0.5 * (FACEtmp1.p3.y + FACEtmp2.p2.y), 0);
            FACEtmp2.p3 = FACEtmp1.p2;

            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            area            = VectorFace(FACEtmp2);
            CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
            break;
        case 2:
            //
            //  Line2 & Line3 only
            //
            if (FACE1.p3.y > FACE2.p0.y) {
                point1 = FACE1.p3;
                point2 = FACE2.p0;
            } else {
                point1 = FACE1.p2;
                point2 = FACE2.p1;
            }
            point3          = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
            area            = CrossProduct(point1 - point3, point2 - point3);
            CalAreaNotQuadr = fabs(area.z) * 0.5;
            break;
        case 3:
            //
            //  Line1 & Line3
            //  Line2 & Line3
            //
            FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp1.p1 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
            if (FACE1.p0.y < FACE2.p0.y) {
                FACEtmp1.p2 = FACE1.p2;
                FACEtmp1.p3 = FACE1.p1;
            } else {
                FACEtmp1.p2 = FACE1.p3;
                FACEtmp1.p3 = FACE1.p0;
            }
            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            break;
        case 4:
            //
            //  Line1 & Line4 only
            //
            if (FACE1.p0.y < FACE2.p3.y) {
                point1 = FACE1.p0;
                point2 = FACE2.p3;
            } else {
                point1 = FACE1.p1;
                point2 = FACE2.p2;
            }
            point3          = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
            area            = CrossProduct(point1 - point3, point2 - point3);
            CalAreaNotQuadr = fabs(area.z) * 0.5;
            break;
        case 5:
            //
            //  Line1 & Line3
            //  Line1 & Line4
            //
            FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp1.p1 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
            if (FACE2.p3.y > FACE1.p0.y) {
                FACEtmp1.p2 = FACE2.p3;
                FACEtmp1.p3 = FACE2.p0;
            } else {
                FACEtmp1.p2 = FACE2.p2;
                FACEtmp1.p3 = FACE2.p1;
            }
            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            break;
        case 8:
            //
            //  Line2 & Line4 only
            //
            FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

            if (FACE1.p0.y > FACE2.p0.y) {
                FACEtmp1.p0 = FACE1.p0;
            } else {
                FACEtmp1.p0 = FACE2.p0;
            }

            if (FACE1.p1.y > FACE2.p1.y) {
                FACEtmp2.p1 = FACE1.p1;
            } else {
                FACEtmp2.p1 = FACE2.p1;
            }

            if (FACE1.p3.y > FACE2.p3.y) {
                FACEtmp1.p3 = FACE2.p3;
            } else {
                FACEtmp1.p3 = FACE1.p3;
            }

            if (FACE1.p2.y > FACE2.p2.y) {
                FACEtmp2.p2 = FACE2.p2;
            } else {
                FACEtmp2.p2 = FACE1.p2;
            }

            FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
                                  0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
            FACEtmp2.p0 = FACEtmp1.p1;

            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            area            = VectorFace(FACEtmp2);
            CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
            break;
        case 9:
            //
            //  Line1 & Line3
            //  Line2 & Line4
            //
            FACEtmp1.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp2.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

            if (FACE1.p0.y > FACE2.p0.y) {
                FACEtmp1.p0 = FACE1.p0;
            } else {
                FACEtmp1.p0 = FACE2.p0;
            }

            if (FACE1.p1.y > FACE2.p1.y) {
                FACEtmp2.p1 = FACE1.p1;
            } else {
                FACEtmp2.p1 = FACE2.p1;
            }

            if (FACE1.p3.y > FACE2.p3.y) {
                FACEtmp1.p3 = FACE2.p3;
            } else {
                FACEtmp1.p3 = FACE1.p3;
            }

            if (FACE1.p2.y > FACE2.p2.y) {
                FACEtmp2.p2 = FACE2.p2;
            } else {
                FACEtmp2.p2 = FACE1.p2;
            }

            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            area            = VectorFace(FACEtmp2);
            CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
            break;
        case 10:
            //
            //  Line2 & Line3
            //  Line2 & Line4
            //
            FACEtmp1.p0 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
            FACEtmp1.p1 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            if (FACE1.p2.y > FACE2.p1.y) {
                FACEtmp1.p2 = FACE2.p1;
                FACEtmp1.p3 = FACE2.p2;
            } else {
                FACEtmp1.p2 = FACE2.p3;
                FACEtmp1.p3 = FACE2.p0;
            }
            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            break;
        case 11:
            //
            //  Line1 & Line3
            //  Line2 & Line3
            //  Line2 & Line4
            //
            FACEtmp1.p0 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            FACEtmp1.p3 = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
            FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

            if (FACE1.p0.y < FACE2.p0.y) {
                FACEtmp2.p1 = FACE1.p1;
                FACEtmp2.p2 = FACE2.p2;
            } else {
                FACEtmp2.p1 = FACE1.p0;
                FACEtmp2.p2 = FACE2.p3;
            }

            FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
                                  0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
            FACEtmp2.p0 = FACEtmp1.p1;

            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            area            = VectorFace(FACEtmp2);
            CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
            break;
        case 12:
            //
            //  Line1 & Line4
            //  Line2 & Line4
            //
            FACEtmp1.p0 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
            FACEtmp1.p1 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            if (FACE1.p2.y > FACE2.p2.y) {
                FACEtmp1.p2 = FACE1.p3;
                FACEtmp1.p3 = FACE1.p0;
            } else {
                FACEtmp1.p2 = FACE1.p2;
                FACEtmp1.p3 = FACE1.p1;
            }
            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            break;
        case 13:
            //
            //  Line1 & Line3
            //  Line1 & Line4
            //  Line2 & Line4
            //
            FACEtmp1.p2 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            FACEtmp2.p1 = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp2.p2 = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
            FACEtmp2.p3 = Point3D(crosspoint[3].x, crosspoint[3].y, 0);

            if (FACE1.p2.y > FACE2.p2.y) {
                FACEtmp1.p0 = FACE2.p0;
                FACEtmp1.p3 = FACE1.p3;
            } else {
                FACEtmp1.p0 = FACE2.p1;
                FACEtmp1.p3 = FACE1.p2;
            }

            FACEtmp1.p1 = Point3D(0.5 * (FACEtmp1.p0.x + FACEtmp2.p1.x),
                                  0.5 * (FACEtmp1.p0.y + FACEtmp2.p1.y), 0);
            FACEtmp2.p0 = FACEtmp1.p1;

            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            area            = VectorFace(FACEtmp2);
            CalAreaNotQuadr = CalAreaNotQuadr + fabs(area.z);
            break;
        case 15:
            //
            //  Line1 & Line3
            //  Line2 & Line3
            //  Line1 & Line4
            //  Line2 & Line4
            //
            FACEtmp1.p0     = Point3D(crosspoint[0].x, crosspoint[0].y, 0);
            FACEtmp1.p1     = Point3D(crosspoint[2].x, crosspoint[2].y, 0);
            FACEtmp1.p2     = Point3D(crosspoint[3].x, crosspoint[3].y, 0);
            FACEtmp1.p3     = Point3D(crosspoint[1].x, crosspoint[1].y, 0);
            area            = VectorFace(FACEtmp1);
            CalAreaNotQuadr = fabs(area.z);
            break;
        default:
            CalAreaNotQuadr = 0;
            break;
    }
    return CalAreaNotQuadr;
}

void ConnGrid::Allocate(const USI& max_neighbor)
{
    nConn   = 0;
    maxConn = max_neighbor;
    halfConn.resize(maxConn);
}

void ConnGrid::AddHalfConn(const OCP_USI& n,
                           const Point3D& area,
                           const Point3D& d,
                           const USI&     direction,
                           const OCP_DBL& flag)
{
    if (nConn >= maxConn) {
        maxConn *= 2;
        halfConn.resize(maxConn);
        // get larger space
        if (maxConn > MAX_NEIGHBOR) {
            // OCP_ABORT("Too many Neighbors!");
        }
    }
    halfConn[nConn].Ad_dd         = area * d / (d * d) * flag;
    halfConn[nConn].d             = d;
    halfConn[nConn].neigh         = n;
    halfConn[nConn].directionType = direction;
    nConn++;

    // cout << n << "   " << direction << "   " << halfConn[nConn-1].Ad_dd << endl;
}

void OCP_COORD::Allocate(const USI& Nx, const USI& Ny, const USI& Nz)
{
    nx      = Nx;
    ny      = Ny;
    nz      = Nz;
    numGrid = nx * ny * nz;

    COORDDATA = new OCP_DBL**[3];
    for (USI i = 0; i < 3; i++) {
        COORDDATA[i] = new OCP_DBL*[2];
        for (USI j = 0; j < 2; j++) {
            COORDDATA[i][j] = new OCP_DBL[(nx + 1) * (ny + 1)];
        }
    }

    ZCORNDATA = new OCP_DBL***[nx];
    for (USI i = 0; i < nx; i++) {
        ZCORNDATA[i] = new OCP_DBL**[ny];
        for (USI j = 0; j < ny; j++) {
            ZCORNDATA[i][j] = new OCP_DBL*[nz];
            for (USI k = 0; k < nz; k++) {
                ZCORNDATA[i][j][k] = new OCP_DBL[8];
            }
        }
    }

    cornerPoints = new Hexahedron**[nx];
    for (USI i = 0; i < nx; i++) {
        cornerPoints[i] = new Hexahedron*[ny];
        for (USI j = 0; j < ny; j++) {
            cornerPoints[i][j] = new Hexahedron[nz];
        }
    }

    v.resize(numGrid);
    depth.resize(numGrid);
    dx.resize(numGrid);
    dy.resize(numGrid);
    dz.resize(numGrid);
    center.resize(numGrid);
}

void OCP_COORD::InputData(const vector<OCP_DBL>& coord, const vector<OCP_DBL>& zcorn)
{
    if (coord.empty() || !InputCOORDDATA(coord)) {
        OCP_ABORT("ERROR COORD!");
    }
    if (zcorn.empty() || !InputZCORNDATA(zcorn)) {
        OCP_ABORT("ERROR ZCORN!");
    }
}

OCP_BOOL OCP_COORD::InputCOORDDATA(const vector<OCP_DBL>& coord)
{
    // See Eclipse -- COORD
    OCP_BOOL flag = OCP_FALSE;
    OCP_USI  iter = 0;

    for (USI J = 0; J < ny + 1; J++) {
        for (USI I = 0; I < nx + 1; I++) {
            // top
            for (USI i = 0; i < 3; i++) {
                COORDDATA[i][0][J * (nx + 1) + I] = coord[iter];
                iter++;
            }
            // bottom
            for (USI i = 0; i < 3; i++) {
                COORDDATA[i][1][J * (nx + 1) + I] = coord[iter];
                iter++;
            }
        }
    }

    flag = OCP_TRUE;
    return flag;
}

OCP_BOOL OCP_COORD::InputZCORNDATA(const vector<OCP_DBL>& zcorn)
{
    // See Eclipse -- ZCORN
    OCP_BOOL flag = OCP_FALSE;
    OCP_USI  iter = 0;

    for (USI K = 0; K < nz; K++) {
        for (USI J = 0; J < ny; J++) {
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][0] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][1] = zcorn[iter];
                iter++;
            }
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][3] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][2] = zcorn[iter];
                iter++;
            }
        }
        for (USI J = 0; J < ny; J++) {
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][4] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][5] = zcorn[iter];
                iter++;
            }
            for (USI I = 0; I < nx; I++) {
                ZCORNDATA[I][J][K][7] = zcorn[iter];
                iter++;
                ZCORNDATA[I][J][K][6] = zcorn[iter];
                iter++;
            }
        }
    }

    flag = OCP_TRUE;
    return flag;
}

void OCP_COORD::SetAllFlags(const HexahedronFace& oFace, const HexahedronFace& Face)
{
    tmpFace = Face;

    if (oFace.p0.z > Face.p0.z + TEENY) {
        tmpFace.p0 = oFace.p0;
        flagp0     = 1;
    } else if (oFace.p0.z < Face.p0.z - TEENY)
        flagp0 = -1;
    else
        flagp0 = 0;

    if (oFace.p1.z > Face.p1.z + TEENY)
        flagp1 = 1;
    else if (oFace.p1.z < Face.p1.z - TEENY) {
        tmpFace.p1 = oFace.p1;
        flagp1     = -1;
    } else
        flagp1 = 0;

    if (oFace.p2.z > Face.p2.z + TEENY)
        flagp2 = 1;
    else if (oFace.p2.z < Face.p2.z - TEENY) {
        tmpFace.p2 = oFace.p2;
        flagp2     = -1;
    } else
        flagp2 = 0;

    if (oFace.p3.z > Face.p3.z + TEENY) {
        tmpFace.p3 = oFace.p3;
        flagp3     = 1;
    } else if (oFace.p3.z < Face.p3.z - TEENY)
        flagp3 = -1;
    else
        flagp3 = 0;

    // check if interface is empty set
    // check if interface is quadrilateral
    // check if the one contains the other one

    if (((oFace.p1.z <= Face.p0.z) && (oFace.p2.z <= Face.p3.z)) ||
        ((oFace.p0.z >= Face.p1.z) && (oFace.p3.z >= Face.p2.z))) {
        flagJump = OCP_TRUE;
    } else {
        flagJump = OCP_FALSE;
        if ((flagp0 * flagp3 >= 0) && (oFace.p0.z <= Face.p1.z) &&
            (oFace.p3.z <= Face.p2.z) && (flagp1 * flagp2 >= 0) &&
            (oFace.p1.z >= Face.p0.z) && (oFace.p2.z >= Face.p3.z)) {
            flagQuad = OCP_TRUE;
        } else {
            flagQuad = OCP_FALSE;
        }
    }
}

void OCP_COORD::SetupCornerPoints()
{
    OCP_USI cindex, oindex; // current block index and the other block index
    OCP_USI nxny = nx * ny;

    // allocate memoery for connections
    vector<ConnGrid> blockconn(numGrid);
    for (OCP_USI iloop = 0; iloop < numGrid; iloop++) {
        blockconn[iloop].Allocate(10);
    }

    // setup each block including coordinates of points, center, depth, and volume
    OCP_DBL xtop, ytop, ztop, xbottom, ybottom, zbottom, xvalue, yvalue, zvalue;
    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                //
                // corner point 0 and 4
                //
                xtop    = COORDDATA[0][0][j * (nx + 1) + i];
                ytop    = COORDDATA[1][0][j * (nx + 1) + i];
                ztop    = COORDDATA[2][0][j * (nx + 1) + i];
                xbottom = COORDDATA[0][1][j * (nx + 1) + i];
                ybottom = COORDDATA[1][1][j * (nx + 1) + i];
                zbottom = COORDDATA[2][1][j * (nx + 1) + i];

                zvalue = ZCORNDATA[i][j][k][0];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p0 = Point3D(xvalue, yvalue, zvalue);

                zvalue = ZCORNDATA[i][j][k][4];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p4 = Point3D(xvalue, yvalue, zvalue);
                //
                //    corner point 1 and 5
                //
                xtop    = COORDDATA[0][0][j * (nx + 1) + i + 1];
                ytop    = COORDDATA[1][0][j * (nx + 1) + i + 1];
                ztop    = COORDDATA[2][0][j * (nx + 1) + i + 1];
                xbottom = COORDDATA[0][1][j * (nx + 1) + i + 1];
                ybottom = COORDDATA[1][1][j * (nx + 1) + i + 1];
                zbottom = COORDDATA[2][1][j * (nx + 1) + i + 1];

                zvalue = ZCORNDATA[i][j][k][1];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p1 = Point3D(xvalue, yvalue, zvalue);

                zvalue = ZCORNDATA[i][j][k][5];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p5 = Point3D(xvalue, yvalue, zvalue);
                //
                //    corner point 2 and 6
                //
                xtop    = COORDDATA[0][0][(j + 1) * (nx + 1) + i + 1];
                ytop    = COORDDATA[1][0][(j + 1) * (nx + 1) + i + 1];
                ztop    = COORDDATA[2][0][(j + 1) * (nx + 1) + i + 1];
                xbottom = COORDDATA[0][1][(j + 1) * (nx + 1) + i + 1];
                ybottom = COORDDATA[1][1][(j + 1) * (nx + 1) + i + 1];
                zbottom = COORDDATA[2][1][(j + 1) * (nx + 1) + i + 1];

                zvalue = ZCORNDATA[i][j][k][2];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p2 = Point3D(xvalue, yvalue, zvalue);

                zvalue = ZCORNDATA[i][j][k][6];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p6 = Point3D(xvalue, yvalue, zvalue);
                //
                //    corner point 3 and 7
                //
                xtop    = COORDDATA[0][0][(j + 1) * (nx + 1) + i];
                ytop    = COORDDATA[1][0][(j + 1) * (nx + 1) + i];
                ztop    = COORDDATA[2][0][(j + 1) * (nx + 1) + i];
                xbottom = COORDDATA[0][1][(j + 1) * (nx + 1) + i];
                ybottom = COORDDATA[1][1][(j + 1) * (nx + 1) + i];
                zbottom = COORDDATA[2][1][(j + 1) * (nx + 1) + i];

                zvalue = ZCORNDATA[i][j][k][3];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p3 = Point3D(xvalue, yvalue, zvalue);

                zvalue = ZCORNDATA[i][j][k][7];
                xvalue =
                    xbottom - (zbottom - zvalue) / (zbottom - ztop) * (xbottom - xtop);
                yvalue =
                    ybottom - (zbottom - zvalue) / (zbottom - ztop) * (ybottom - ytop);
                cornerPoints[i][j][k].p7 = Point3D(xvalue, yvalue, zvalue);

                //    calculate volumes and pore volumes
                cindex = k * nxny + j * nx + i;
                //
                // NOTE: if there are several points not well ordered, the calculated
                // volume will be negative.
                //
                v[cindex]      = VolumHexahedron(cornerPoints[i][j][k]); // NTG
                v[cindex]      = fabs(v[cindex]);
                center[cindex] = CenterHexahedron(cornerPoints[i][j][k]);
                depth[cindex]  = center[cindex].z;
            }
        }
    }

    // find neighbor and calculate transmissibility
    OCP_USI num_conn = 0; // record the num of connection, a->b & b->a are both included
    Point3D Pcenter, Pface, Pc2f; // center of Hexahedron
    HexahedronFace Face, oFace;   // current face, the other face
    HexahedronFace FaceP, oFaceP; // Projection of Face and the other face
    Point3D        areaV;         // area vector of interface
    OCP_DBL        areaP;         // area of projection of interface
    OCP_INT        iznnc;
    Point3D        dxpoint, dypoint, dzpoint;

    // test
    // cornerPoints[13][1][72].p0; cornerPoints[13][1][72].p1;
    // cornerPoints[13][1][72].p2; cornerPoints[13][1][72].p3;
    // cornerPoints[13][1][72].p4; cornerPoints[13][1][72].p5;
    // cornerPoints[13][1][72].p6; cornerPoints[13][1][72].p7;

    // cornerPoints[13][2][74].p0; cornerPoints[13][2][74].p1;
    // cornerPoints[13][2][74].p2; cornerPoints[13][2][74].p3;
    // cornerPoints[13][2][74].p4; cornerPoints[13][2][74].p5;
    // cornerPoints[13][2][74].p6; cornerPoints[13][2][74].p7;

    /////////////////////////////////////////////////////////////////////
    // Attention that The coordinate axis follows the right-hand rule ! //
    /////////////////////////////////////////////////////////////////////
    //
    //      o----> x
    //     /|
    //    y z
    // For a face, p0 and p3 are the points of upper edge of quadrilateral,
    // p1 and p2 are the points of lower edge
    //       p0 ---- p3
    //        |       |
    //        |       |
    //       p1 ---- p2

    // Determine flagForward
    if (COORDDATA[1][0][nx + 1] > COORDDATA[1][0][0])
        flagForward = 1.0;
    else
        flagForward = -1.0;

    for (USI k = 0; k < nz; k++) {
        for (USI j = 0; j < ny; j++) {
            for (USI i = 0; i < nx; i++) {
                // begin from each block
                const Hexahedron& block = cornerPoints[i][j][k];
                cindex                  = k * nxny + j * nx + i;
                Pcenter                 = center[cindex];

                // cout << "============= " << cindex << " =============" << endl;
                //
                // (x-) direction
                //

                Face.p0 = block.p0;
                Face.p1 = block.p4;
                Face.p2 = block.p7;
                Face.p3 = block.p3;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dxpoint = Pc2f;

                if (i == 0) {
                    // nothing to do
                } else {

                    const Hexahedron& leftblock = cornerPoints[i - 1][j][k];
                    oindex                      = k * nxny + j * nx + i - 1;

                    oFace.p0 = leftblock.p1;
                    oFace.p1 = leftblock.p5;
                    oFace.p2 = leftblock.p6;
                    oFace.p3 = leftblock.p2;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = VectorFace(tmpFace);
                        } else {
                            FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                            FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                            FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                            FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                            oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                            oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                            oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                            oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = VectorFace(Face);
                            // correct
                            if (fabs(areaV.x) < 1E-6) {
                                OCP_WARNING("x is too small");
                            } else {
                                areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                areaV.x = OCP_SIGN(areaV.x) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    // check if upNNC and downNNC exist
                    if ((flagp0 > 0) || (flagp3 > 0))
                        upNNC = OCP_TRUE;
                    else
                        upNNC = OCP_FALSE;
                    if ((flagp1 < 0) || (flagp2 < 0))
                        downNNC = OCP_TRUE;
                    else
                        downNNC = OCP_FALSE;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& leftblock = cornerPoints[i - 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i - 1;
                        oFace.p0 = leftblock.p1;
                        oFace.p1 = leftblock.p5;
                        oFace.p2 = leftblock.p6;
                        oFace.p3 = leftblock.p2;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;
                        if ((flagp0 > 0) || (flagp3 > 0))
                            upNNC = OCP_TRUE;
                        else
                            upNNC = OCP_FALSE;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& leftblock = cornerPoints[i - 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i - 1;
                        oFace.p0 = leftblock.p1;
                        oFace.p1 = leftblock.p5;
                        oFace.p2 = leftblock.p6;
                        oFace.p3 = leftblock.p2;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;

                        if ((flagp1 < 0) || (flagp2 < 0))
                            downNNC = OCP_TRUE;
                        else
                            downNNC = OCP_FALSE;
                    }
                }

                //
                // (x+) direction
                //
                Face.p0 = block.p2;
                Face.p1 = block.p6;
                Face.p2 = block.p5;
                Face.p3 = block.p1;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dxpoint = Pc2f - dxpoint;

                if (i == nx - 1) {
                    // nothing to do
                } else {

                    const Hexahedron& rightblock = cornerPoints[i + 1][j][k];
                    oindex                       = k * nxny + j * nx + i + 1;

                    oFace.p0 = rightblock.p3;
                    oFace.p1 = rightblock.p7;
                    oFace.p2 = rightblock.p4;
                    oFace.p3 = rightblock.p0;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = VectorFace(tmpFace);
                        } else {
                            FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                            FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                            FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                            FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                            oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                            oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                            oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                            oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = VectorFace(Face);
                            // correct
                            if (fabs(areaV.x) < 1E-6) {
                                OCP_WARNING("x is too small");
                            } else {
                                areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                areaV.x = OCP_SIGN(areaV.x) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    if ((flagp0 > 0) || (flagp3 > 0))
                        upNNC = OCP_TRUE;
                    else
                        upNNC = OCP_FALSE;
                    if ((flagp1 < 0) || (flagp2 < 0))
                        downNNC = OCP_TRUE;
                    else
                        downNNC = OCP_FALSE;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& rightblock =
                            cornerPoints[i + 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i + 1;
                        oFace.p0 = rightblock.p3;
                        oFace.p1 = rightblock.p7;
                        oFace.p2 = rightblock.p4;
                        oFace.p3 = rightblock.p0;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;

                        if ((flagp0 > 0) || (flagp3 > 0))
                            upNNC = OCP_TRUE;
                        else
                            upNNC = OCP_FALSE;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& rightblock =
                            cornerPoints[i + 1][j][k + iznnc];
                        oindex   = (k + iznnc) * nxny + j * nx + i + 1;
                        oFace.p0 = rightblock.p3;
                        oFace.p1 = rightblock.p7;
                        oFace.p2 = rightblock.p4;
                        oFace.p3 = rightblock.p0;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p3.y, Face.p3.z, 0);
                                FaceP.p1  = Point3D(Face.p0.y, Face.p0.z, 0);
                                FaceP.p2  = Point3D(Face.p1.y, Face.p1.z, 0);
                                FaceP.p3  = Point3D(Face.p2.y, Face.p2.z, 0);
                                oFaceP.p0 = Point3D(oFace.p3.y, oFace.p3.z, 0);
                                oFaceP.p1 = Point3D(oFace.p0.y, oFace.p0.z, 0);
                                oFaceP.p2 = Point3D(oFace.p1.y, oFace.p1.z, 0);
                                oFaceP.p3 = Point3D(oFace.p2.y, oFace.p2.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.x) < 1E-6) {
                                    OCP_WARNING("x is too small");
                                } else {
                                    areaV.y = areaV.y / fabs(areaV.x) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.x) * areaP;
                                    areaV.x = OCP_SIGN(areaV.x) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 1,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;

                        if ((flagp1 < 0) || (flagp2 < 0))
                            downNNC = OCP_TRUE;
                        else
                            downNNC = OCP_FALSE;
                    }
                }

                //
                // (y-) direction
                //
                Face.p0 = block.p1;
                Face.p1 = block.p5;
                Face.p2 = block.p4;
                Face.p3 = block.p0;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dypoint = Pc2f;

                if (j == 0) {
                    // nothing to do
                } else {

                    const Hexahedron& backblock = cornerPoints[i][j - 1][k];
                    oindex                      = k * nxny + (j - 1) * nx + i;

                    oFace.p0 = backblock.p2;
                    oFace.p1 = backblock.p6;
                    oFace.p2 = backblock.p7;
                    oFace.p3 = backblock.p3;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = VectorFace(tmpFace);
                        } else {
                            FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                            FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                            FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                            FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                            oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                            oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                            oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                            oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = VectorFace(Face);
                            // correct
                            if (fabs(areaV.y) < 1E-6) {
                                OCP_WARNING("y is too small");
                            } else {
                                areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                areaV.y = OCP_SIGN(areaV.y) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    if ((flagp0 > 0) || (flagp3 > 0))
                        upNNC = OCP_TRUE;
                    else
                        upNNC = OCP_FALSE;
                    if ((flagp1 < 0) || (flagp2 < 0))
                        downNNC = OCP_TRUE;
                    else
                        downNNC = OCP_FALSE;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& backblock = cornerPoints[i][j - 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j - 1) * nx + i;
                        oFace.p0 = backblock.p2;
                        oFace.p1 = backblock.p6;
                        oFace.p2 = backblock.p7;
                        oFace.p3 = backblock.p3;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;

                        if ((flagp0 > 0) || (flagp3 > 0))
                            upNNC = OCP_TRUE;
                        else
                            upNNC = OCP_FALSE;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& backblock = cornerPoints[i][j - 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j - 1) * nx + i;
                        oFace.p0 = backblock.p2;
                        oFace.p1 = backblock.p6;
                        oFace.p2 = backblock.p7;
                        oFace.p3 = backblock.p3;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;

                        if ((flagp1 < 0) || (flagp2 < 0))
                            downNNC = OCP_TRUE;
                        else
                            downNNC = OCP_FALSE;
                    }
                }

                //
                // (y+) direction
                //
                Face.p0 = block.p3;
                Face.p1 = block.p7;
                Face.p2 = block.p6;
                Face.p3 = block.p2;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dypoint = Pc2f - dypoint;

                if (j == ny - 1) {
                    // nothing to do
                } else {

                    const Hexahedron& frontblock = cornerPoints[i][j + 1][k];
                    oindex                       = k * nxny + (j + 1) * nx + i;

                    oFace.p0 = frontblock.p0;
                    oFace.p1 = frontblock.p4;
                    oFace.p2 = frontblock.p5;
                    oFace.p3 = frontblock.p1;

                    SetAllFlags(oFace, Face);

                    // calculate the interface of two face
                    if (flagJump) {
                        // nothing to do
                    } else {
                        if (flagQuad) {
                            areaV = VectorFace(tmpFace);
                        } else {
                            FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                            FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                            FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                            FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                            oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                            oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                            oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                            oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                            areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                            // attention the direction of vector
                            areaV = VectorFace(Face);
                            // correct
                            if (fabs(areaV.y) < 1E-6) {
                                OCP_WARNING("y is too small");
                            } else {
                                areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                areaV.y = OCP_SIGN(areaV.y) * areaP;
                            }
                        }
                        blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                      flagForward);
                        num_conn++;
                    }

                    // then find all NNC for current block
                    if ((flagp0 > 0) || (flagp3 > 0))
                        upNNC = OCP_TRUE;
                    else
                        upNNC = OCP_FALSE;
                    if ((flagp1 < 0) || (flagp2 < 0))
                        downNNC = OCP_TRUE;
                    else
                        downNNC = OCP_FALSE;

                    iznnc = -1;
                    while (upNNC) {
                        // if (-iznnc > k) break;
                        if (-iznnc - static_cast<OCP_INT>(k) > 0) break;
                        // find object block
                        const Hexahedron& frontblock =
                            cornerPoints[i][j + 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j + 1) * nx + i;
                        oFace.p0 = frontblock.p0;
                        oFace.p1 = frontblock.p4;
                        oFace.p2 = frontblock.p5;
                        oFace.p3 = frontblock.p1;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc--;

                        if ((flagp0 > 0) || (flagp3 > 0))
                            upNNC = OCP_TRUE;
                        else
                            upNNC = OCP_FALSE;
                    }

                    iznnc = 1;
                    while (downNNC) {
                        if (k + iznnc > nz - 1) break;
                        // find object block
                        const Hexahedron& frontblock =
                            cornerPoints[i][j + 1][k + iznnc];
                        oindex   = (k + iznnc) * nxny + (j + 1) * nx + i;
                        oFace.p0 = frontblock.p0;
                        oFace.p1 = frontblock.p4;
                        oFace.p2 = frontblock.p5;
                        oFace.p3 = frontblock.p1;

                        SetAllFlags(oFace, Face);

                        // calculate the interface of two face
                        if (flagJump) {
                            // nothing to do
                        } else {
                            if (flagQuad) {
                                areaV = VectorFace(tmpFace);
                            } else {
                                FaceP.p0  = Point3D(Face.p0.x, Face.p0.z, 0);
                                FaceP.p1  = Point3D(Face.p3.x, Face.p3.z, 0);
                                FaceP.p2  = Point3D(Face.p2.x, Face.p2.z, 0);
                                FaceP.p3  = Point3D(Face.p1.x, Face.p1.z, 0);
                                oFaceP.p0 = Point3D(oFace.p0.x, oFace.p0.z, 0);
                                oFaceP.p1 = Point3D(oFace.p3.x, oFace.p3.z, 0);
                                oFaceP.p2 = Point3D(oFace.p2.x, oFace.p2.z, 0);
                                oFaceP.p3 = Point3D(oFace.p1.x, oFace.p1.z, 0);
                                areaP     = CalAreaNotQuadr(FaceP, oFaceP);
                                // attention the direction of vector
                                areaV = VectorFace(Face);
                                // correct
                                if (fabs(areaV.y) < 1E-6) {
                                    OCP_WARNING("y is too small");
                                } else {
                                    areaV.x = areaV.x / fabs(areaV.y) * areaP;
                                    areaV.z = areaV.z / fabs(areaV.y) * areaP;
                                    areaV.y = OCP_SIGN(areaV.y) * areaP;
                                }
                            }
                            blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 2,
                                                          flagForward);
                            num_conn++;
                        }
                        iznnc++;

                        if ((flagp1 < 0) || (flagp2 < 0))
                            downNNC = OCP_TRUE;
                        else
                            downNNC = OCP_FALSE;
                    }
                }

                //
                // (z-) direction
                //
                Face.p0 = block.p0;
                Face.p1 = block.p3;
                Face.p2 = block.p2;
                Face.p3 = block.p1;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dzpoint = Pc2f;
                if (k == 0) {
                    // nothing to do
                } else {
                    // upblock
                    oindex = (k - 1) * nxny + j * nx + i;

                    tmpFace = Face;
                    areaV   = VectorFace(tmpFace);
                    blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 3, flagForward);
                    num_conn++;
                }

                //
                // (z+) direction
                //
                Face.p0 = block.p5;
                Face.p1 = block.p6;
                Face.p2 = block.p7;
                Face.p3 = block.p4;
                Pface   = CenterFace(Face);
                Pc2f    = Pface - Pcenter;
                dzpoint = Pc2f - dzpoint;

                if (k == nz - 1) {
                    // nothing to do
                } else {
                    // downblock
                    oindex = (k + 1) * nxny + j * nx + i;

                    tmpFace = Face;
                    areaV   = VectorFace(tmpFace);
                    blockconn[cindex].AddHalfConn(oindex, areaV, Pc2f, 3, flagForward);
                    num_conn++;
                }

                // calculate dx,dy,dz
                dx[cindex] = sqrt(dxpoint.x * dxpoint.x + dxpoint.y * dxpoint.y +
                                  dxpoint.z * dxpoint.z);
                dy[cindex] = sqrt(dypoint.x * dypoint.x + dypoint.y * dypoint.y +
                                  dypoint.z * dypoint.z);
                dz[cindex] = sqrt(dzpoint.x * dzpoint.x + dzpoint.y * dzpoint.y +
                                  dzpoint.z * dzpoint.z);

                OCP_ASSERT(isfinite(dx[cindex]), "Wrong dx!");
                OCP_ASSERT(isfinite(dy[cindex]), "Wrong dy!");
                OCP_ASSERT(isfinite(dz[cindex]), "Wrong dz!");
            }
        }
    }

    OCP_ASSERT(num_conn % 2 == 0, "Wrong Conn!");
    numConnMax = num_conn / 2;
    connect.resize(numConnMax);
    //
    //    calculate the x,y,z direction transmissibilities of each block and save them
    //
    // make the connections
    OCP_USI iter_conn = 0;
    for (OCP_USI n = 0; n < numGrid; n++) {
        for (USI j = 0; j < blockconn[n].nConn; j++) {
            OCP_USI nn = blockconn[n].halfConn[j].neigh;
            if (nn < n) continue;
            USI jj;
            for (jj = 0; jj < blockconn[nn].nConn; jj++) {
                if (blockconn[nn].halfConn[jj].neigh == n) {
                    break;
                }
            }
            if (jj == blockconn[nn].nConn) {
                continue;
            }
            if (blockconn[n].halfConn[j].Ad_dd <= 0 ||
                blockconn[nn].halfConn[jj].Ad_dd <= 0) {
                // OCP_FALSE connection
                continue;
            }

            //
            // now, blockconn[n].halfConn[j]
            //     blockconn[nn].halfConn[jj]
            //     are a pair of connections
            connect[iter_conn].begin         = n;
            connect[iter_conn].Ad_dd_begin   = blockconn[n].halfConn[j].Ad_dd;
            connect[iter_conn].end           = nn;
            connect[iter_conn].Ad_dd_end     = blockconn[nn].halfConn[jj].Ad_dd;
            connect[iter_conn].directionType = blockconn[n].halfConn[j].directionType;
            iter_conn++;
        }
    }
    numConn = iter_conn;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/19/2021      Create file                          */
/*----------------------------------------------------------------------------*/