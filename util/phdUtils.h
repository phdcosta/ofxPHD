#pragma once

#include "ofMain.h"
//#include "phdDefinitions.h"

//-------------------------------------------------------------------------------------------------
enum phdStretchMode { psmNone, psmCentered, psmStretched, psmFitted, psmFitCentered };

//-------------------------------------------------------------------------------------------------
void printDebug(string & _str);

// rotate a point
//-------------------------------------------------------------------------------------------------
ofVec2f rotatePoint(float x, float y, float cx, float cy, float ccx, float ccy, float sine, float cosine);

// calculate tangents between two circles
//-------------------------------------------------------------------------------------------------
void circlesTangent(float _x1, float _y1, float _r1, float _x2, float _y2, float _r2, vector<ofVec2f> & _tangents);

// calculate intersection points between two circles
//-------------------------------------------------------------------------------------------------
void circlesIntersection(double x0, double y0, double r0, double x1, double y1, double r1, vector<ofVec2f> & _points);

// calculate triangle points from segment sizes
//-------------------------------------------------------------------------------------------------
bool getTriangleFromSegmentSize(double _szA, double _szB, double _szC, double x0, double y0, double _angle, double _scale, ofVec2f & _ptA, ofVec2f & _ptB, ofVec2f & _ptC);

// compute the area of a triangle using Heron's formula
//-------------------------------------------------------------------------------------------------
double triareaHeron(double a, double b, double c);

// compute the distance between two 2d points
//-------------------------------------------------------------------------------------------------
double dist(double x0, double y0, double z0, double x1, double y1, double z1);

// calculate barycentric coordinates
// triangle 1st vertex: x0,y0,z0
// triangle 2nd vertex: x1,y1,z1
//  triangle 3rd vertex: x2,y2,z2
// point inside triangle: vx, vy,vz
// *u,*v,*w are the coordinates returned
void barycent(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double vx, double vy, double vz, double *u, double *v, double *w);

void setIdentityMatrix(ofMatrix4x4 & _m);
void setRotateMatrixX(ofMatrix4x4 & _m, double angle);
void setRotateMatrixY(ofMatrix4x4 & _m, double angle);
void setRotateMatrixZ(ofMatrix4x4 & _m, double angle);
void setRotateMatrixXYZ(ofMatrix4x4 & _m, double angleX, double angleY, double angleZ);
void setScaleMatrix(ofMatrix4x4 & _m, double scaleX, double scaleY, double scaleZ);
void setTranslateMatrix(ofMatrix4x4 & _m, double moveX, double moveY, double moveZ);

void setTransformationMatrix(ofMatrix4x4 & _m, ofVec2f & _curPos, float _curAngle, ofVec2f & _curScale, ofVec2f & _newPos, float _newAngle, ofVec2f & _newScale);

//--------------------------------------------------------------
// Code by Theo from URL above
//--------------------------------------------------------------
// http://www.openframeworks.cc/forum/viewtopic.php?f=9&t=1443
//--------------------------------------------------------------
bool phdPointInsidePolygon(vector<ofVec2f> & _points, float px, float py);

bool phdPointInsideTriangle(float _x, float _y, float _xa, float _ya, float _xb, float _yb, float _xc, float _yc);

double phdDistancePointLineSquared(double x, double y, double x1, double y1, double x2, double y2);

bool phdPointNearEdge(double x1, double y1, double x2, double y2, double px, double py, double _dist);

int phdPointNearEdge(vector<ofVec2f> & _points, double px, double py, double _dist, bool _closed);

int phdPointOverVertex(vector<ofVec2f> & _points, double px, double py);

double phdGetCentroid(vector<ofVec2f> & _points, ofVec2f & _centroid);

void getTransformedPolyline(vector<ofVec2f> & _src, vector<ofVec2f> & _dst, ofMatrix4x4 _mat);

// monotone implementation of convex hull contour
//--------------------------------------------------------------
typedef long long CoordType;

struct phdConvexPoint {
	CoordType x, y;
 	bool operator <(const phdConvexPoint &p) const {
		return x < p.x || (x == p.x && y < p.y);
	}
};
 
void contourToConvexHull(vector<ofVec2f> &src, vector<ofVec2f> &dst);
//--------------------------------------------------------------

// area of triangle from its vertices
double triareaVertices(double _ax, double _ay, double _bx, double _by, double _cx, double _cy);

void drawEdgeWithSizeLabel(double _x1, double _y1, double _x2, double _y2, string _label, bool _edge);

void drawGrid(float _xGap, float _yGap, ofTrueTypeFont * _font);

void drawFilledBorderRectangle(float _x1, float _y1, float _w, float _h, ofColor & _fill, ofColor & _border);

void drawBorderRectangle(float _x1, float _y1, float _w, float _h, ofColor & _border);

void drawBorderRectangle(float _x1, float _y1, float _w, float _h);
//--------------------------------------------------------------------------------------------------------------
bool isClickEqualsDrag(ofMouseEventArgs & onClick, ofMouseEventArgs & onDrag);

void bezierSplinePoints(vector<ofVec2f> pnts, int count, int segments, vector<ofVec2f> & points);

//--------------------------------------------------------------------------------------------------------------
float getBestHeightbyFixedW(float _width, float _height, float _minW);
float getMinScale(float _dstW, float _dstH, float _srcW, float _srcH);
void getBestSizeWH(float _dstW, float _dstH, float _srcW, float _srcH, float & _scaledW, float & _scaledH);
void getStrecthModeArea(ofRectangle & area, float srcW, float srcH, float dstW, float dstH, phdStretchMode stretchMode);
ofRectangle getStrecthModeArea(float _x, float _y, float _w, float _h, float srcW, float srcH, phdStretchMode stretchMode);

enum phdTriangleDirection { ptdLeft, ptdRight, ptdTop, ptdBottom };
void drawBorderTriangle(float _x, float _y, float _w, float _h, phdTriangleDirection _dir = ptdTop);
void drawFilledTriangle(float _x, float _y, float _w, float _h, phdTriangleDirection _dir);

void genColorsFromBase(ofColor _base, ofColor & _normal, ofColor & _focused, ofColor & _selected, float _fNormal = 1.0, float _fFocused = 0.75, float _fSelected = 1.25, bool _alpha = false);

//--------------------------------------------------------------------------------------------------------------
void rotSliderV(float _x, float _y, float _w, float _h, float _value, ofColor _color, bool _filled = true);
void rotSliderV2(float _x, float _y, float _w, float _h, float _value, ofColor _clrA, ofColor _clrB, float _aIni = 135.0, float _aSize = 270.0);

double rotAngleToValue(float _angle, float _aIni = 135.0, float _aSize = 270.0, float _mark = 90.0);
void rotSliderA(float _x, float _y, float _w, float _h, float _angle, ofColor _clrA, ofColor _clrB);

double atanPH(double _x, double _y, double _cx, double _cy);

float mapEaseQuartic(float value, float sizeIn, float sizeOut);

//--------------------------------------------------------------------------------------------------------------
void genTriangle(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices);
void genDiamond(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices);
void genQuad(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices);

void adjustRectangle(ofPoint & _center, ofPoint & _maxSize, ofRectangle & _resultArea);

//---------------------------------------------------------------------------------------------
// return projection v1 on to v2
ofVec2f projection(const ofVec2f & v1, const ofVec2f & v2);

//---------------------------------------------------------------------------------------------
void phdMakeTransformationMatrix(ofMatrix4x4 & _mat, float _cx, float _cy, float _angFrom, float _angTo, float _dx, float _dy, float _sx, float _sy);

//---------------------------------------------------------------------------------------------
void getTransformationSinCos(float _aFrom, float _aTo, float & _sinF, float & _cosF, float & _sinT, float & _cosT);

//---------------------------------------------------------------------------------------------
void transformPoint(ofVec2f & _pt, double _cx, double _cy, double _sinF, double _cosF, double _sinT, double _cosT, double _dx, double _dy, double _sx, double _sy);

//---------------------------------------------------------------------------------------------
void fill4(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border4(float _x, float _y, float _w, float _h);
void fill4d(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border4d(float _x, float _y, float _w, float _h);
void fill3u(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border3u(float _x, float _y, float _w, float _h);
void fill3d(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border3d(float _x, float _y, float _w, float _h);
void fill3r(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border3r(float _x, float _y, float _w, float _h);
void fill3l(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border3l(float _x, float _y, float _w, float _h);
void fill8(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src);
void border8(float _x, float _y, float _w, float _h);
