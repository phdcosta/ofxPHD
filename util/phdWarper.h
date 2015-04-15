#pragma once

#include "ofMain.h"

//--------------------------------------------------------------------
// extracted from Dan Wilcox work
// available at https://github.com/danomatika/ofxAppUtils
//--------------------------------------------------------------------

/*
Copyright (c) 2011-2012 Dan Wilcox <danomatika@gmail.com> All rights reserved.

The following terms (the "Standard Improved BSD License") apply to all files
associated with the software unless explicitly disclaimed in individual files:

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
   and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//------------------------------------------------------------------
// 3x3 matrix manipulation routines for
// projection warping from a modified Theo example on the OF forums:
// http://threeblindmiceandamonkey.com/?p=31
//------------------------------------------------------------------

// multiply matrix: c = a * b
inline void multiplyMatrix(double a[3][3], double b[3][3], double c[3][3]) {
    c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
    c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
    c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
    c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
    c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
    c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
    c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
    c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
    c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
}

// determinant of a 2x2 matrix
inline double det2(double a, double b, double c, double d) { return (a*d - b*c); }

// adjoint matrix: b = adjoint(a); returns determinant(a)
inline double adjointMatrix(double a[3][3], double b[3][3]) {
    b[0][0] = det2(a[1][1], a[1][2], a[2][1], a[2][2]);
    b[1][0] = det2(a[1][2], a[1][0], a[2][2], a[2][0]);
    b[2][0] = det2(a[1][0], a[1][1], a[2][0], a[2][1]);
    b[0][1] = det2(a[2][1], a[2][2], a[0][1], a[0][2]);
    b[1][1] = det2(a[2][2], a[2][0], a[0][2], a[0][0]);
    b[2][1] = det2(a[2][0], a[2][1], a[0][0], a[0][1]);
    b[0][2] = det2(a[0][1], a[0][2], a[1][1], a[1][2]);
    b[1][2] = det2(a[0][2], a[0][0], a[1][2], a[1][0]);
    b[2][2] = det2(a[0][0], a[0][1], a[1][0], a[1][1]);
    return a[0][0]*b[0][0] + a[0][1]*b[0][1] + a[0][2]*b[0][2];
}

#define MATRIX_TOLERANCE 1e-13
#define MATRIX_ZERO(x) ((x)<MATRIX_TOLERANCE && (x)>-MATRIX_TOLERANCE)

// calculate matrix for unit square to quad mapping
inline void mapSquareToQuad(double quad[4][2],  // vertices of quadrilateral
                     double SQ[3][3])    // square->quad transform
{
    double px, py;

    px = quad[0][0]-quad[1][0]+quad[2][0]-quad[3][0];
    py = quad[0][1]-quad[1][1]+quad[2][1]-quad[3][1];

    if (MATRIX_ZERO(px) && MATRIX_ZERO(py)) {
        SQ[0][0] = quad[1][0]-quad[0][0];
        SQ[1][0] = quad[2][0]-quad[1][0];
        SQ[2][0] = quad[0][0];
        SQ[0][1] = quad[1][1]-quad[0][1];
        SQ[1][1] = quad[2][1]-quad[1][1];
        SQ[2][1] = quad[0][1];
        SQ[0][2] = 0.;
        SQ[1][2] = 0.;
        SQ[2][2] = 1.;
        return;
    } else {
        double dx1, dx2, dy1, dy2, del;

        dx1 = quad[1][0]-quad[2][0];
        dx2 = quad[3][0]-quad[2][0];
        dy1 = quad[1][1]-quad[2][1];
        dy2 = quad[3][1]-quad[2][1];
        del = det2(dx1,dx2, dy1,dy2);

        SQ[0][2] = det2(px,dx2, py,dy2)/del;
        SQ[1][2] = det2(dx1,px, dy1,py)/del;
        SQ[2][2] = 1.;
        SQ[0][0] = quad[1][0]-quad[0][0]+SQ[0][2]*quad[1][0];
        SQ[1][0] = quad[3][0]-quad[0][0]+SQ[1][2]*quad[3][0];
        SQ[2][0] = quad[0][0];
        SQ[0][1] = quad[1][1]-quad[0][1]+SQ[0][2]*quad[1][1];
        SQ[1][1] = quad[3][1]-quad[0][1]+SQ[1][2]*quad[3][1];
        SQ[2][1] = quad[0][1];
    }
}

// calculate matrix for general quad to quad mapping
inline void mapQuadToQuad( double in[4][2],    // starting quad
                    double out[4][2],   // target quad
                    double ST[3][3])    // the matrix (returned)
{
    double quad[4][2], MS[3][3];
    double SM[3][3], MT[3][3];

    quad[0][0] = in[0][0]; quad[0][1] = in[0][1];
    quad[1][0] = in[1][0]; quad[1][1] = in[1][1];
    quad[2][0] = in[2][0]; quad[2][1] = in[2][1];
    quad[3][0] = in[3][0]; quad[3][1] = in[3][1];
    mapSquareToQuad(quad, MS);
    adjointMatrix(MS, SM);

    quad[0][0] = out[0][0] ; quad[0][1] = out[0][1] ;
    quad[1][0] = out[1][0] ; quad[1][1] = out[1][1] ;
    quad[2][0] = out[2][0] ; quad[2][1] = out[2][1] ;
    quad[3][0] = out[3][0] ; quad[3][1] = out[3][1] ;
    mapSquareToQuad(quad, MT);

    multiplyMatrix(SM, MT, ST);
}

//---------------------------------------------------------------------------------------------------------------
// calculate transformation matrix
//---------------------------------------------------------------------------------------------------------------
// projection warping from a modified Theo example on the OF forums:
// http://threeblindmiceandamonkey.com/?p=31
//---------------------------------------------------------------------------------------------------------------

inline void makeIdentity(GLfloat glWarpMatrix[16]) {
	// we set it to the default - 0 translation
	// and 1.0 scale for x y z and w
	for(int i = 0; i < 16; i++) {
		if(i % 5 != 0) {
			glWarpMatrix[i] = 0.0;
		} else {
			glWarpMatrix[i] = 1.0;
		}
	}
}

inline void calcWarpMatrix(double quadIn[4][2], double quadOut[4][2], GLfloat glWarpMatrix[16]) {

	// projection warping matrice
	double	_warpMatrix[3][3];

	// makes a identity matrix
	makeIdentity(glWarpMatrix);

	// perform the warp calculation SRC to DST
	mapQuadToQuad(quadIn, quadOut, _warpMatrix);

	// copy the values
	glWarpMatrix[0]	= _warpMatrix[0][0];
	glWarpMatrix[1]	= _warpMatrix[0][1];
	glWarpMatrix[3]	= _warpMatrix[0][2];

	glWarpMatrix[4]	= _warpMatrix[1][0];
	glWarpMatrix[5]	= _warpMatrix[1][1];
	glWarpMatrix[7]	= _warpMatrix[1][2];

	glWarpMatrix[12] = _warpMatrix[2][0];
	glWarpMatrix[13] = _warpMatrix[2][1];
	glWarpMatrix[15] = _warpMatrix[2][2];
}

//---------------------------------------------------------------------------------------------------------------
// here begins PHD code
//---------------------------------------------------------------------------------------------------------------
inline void getWarpedPoint(float inX, float inY, float & outX, float & outY, GLfloat glWarpMatrix[16]) {

	GLfloat z = 0.0;
	GLfloat r0 = 0.0; GLfloat r1 = 0.0; GLfloat r2 = 0.0; GLfloat r3 = 0.0;
	r0 = (float)(inX * glWarpMatrix[0] + inY*glWarpMatrix[4] + z*glWarpMatrix[8]  + 1.0*glWarpMatrix[12]);
	r1 = (float)(inX * glWarpMatrix[1] + inY*glWarpMatrix[5] + z*glWarpMatrix[9]  + 1.0*glWarpMatrix[13]);
	r2 = (float)(inX * glWarpMatrix[2] + inY*glWarpMatrix[6] + z*glWarpMatrix[10] + 1.0*glWarpMatrix[14]);
	r3 = (float)(inX * glWarpMatrix[3] + inY*glWarpMatrix[7] + z*glWarpMatrix[11] + 1.0*glWarpMatrix[15]);

	outX = outY = 0.0;

	if(r3 != 0.0) {
		outX = (float)(r0 / r3);
		outY = (float)(r1 / r3);
	}
}

//---------------------------------------------------------------------------------------------------------------
enum phdCurveInterpolatorMode { pcimLinear, pcimBSpline, pcimCatmullRom };

inline void getInterpolatedPoints(vector<ofVec2f> & ctrls, vector<ofVec2f> & curve, int resolution = 2) {
	float _res = 1.0 / resolution;
	for(int i = 1; i < ctrls.size(); i++) {
		for(int j = 0; j < resolution; j++) {
			float x = ctrls[i-1].x + (ctrls[i].x - ctrls[i-1].x) * ((float)j * _res);
			float y = ctrls[i-1].y + (ctrls[i].y - ctrls[i-1].y) * ((float)j * _res);
			curve.push_back(ofVec2f(x,y));
		}
		if(i == ctrls.size()-1) {
			curve.push_back(ctrls[i]); // last point
		}
	}
}

//---------------------------------------------------------------------------------------------------------------
float linear(float mu, float v0, float v1) {
	return (v1-v0) * mu + v0;
}

//---------------------------------------------------------------------------------------------------------------
float cubic(float mu, float v0, float v1, float v2, float v3) {
	float mu2 = mu*mu;
	float mu3 = mu2*mu;
	float a0 = v3 - v2 - v0 + v1;
	float a1 = v0 - v1 - a0;
	float a2 = v2 - v0;
	float a3 = v1;
	return (a0 * mu3 + a1 * mu2 + a2 * mu + a3);
}

//---------------------------------------------------------------------------------------------------------------
inline float bSpline(float mu, float p0, float p1, float p2, float p3) {
	double mu2 = mu * mu;
	double mu3 = mu * mu2;
	return ((1.0/6.0)*
			(mu3*(-p0+3.0*p1-3.0*p2+p3)+
			 mu2*(3.0*p0-6.0*p1+3.0*p2)+
			 mu* (-3.0*p0+3.0*p2)+(p0+4.0*p1+p2)));
}

inline void getBSplineCurvePoints(vector<ofVec2f> & ctrls, vector<ofVec2f> & curve, int resolution = 2) {

	if(ctrls.size() < 2) return;

	vector<ofVec2f> ctrlPts = ctrls;

	// Phantom Points
	ofVec2f pha, phb;
	pha.set(2.0*ctrls[0].x-ctrls[1].x, 2.0*ctrls[0].y-ctrls[1].y);
	phb.set(2.0*ctrls[ctrls.size()-1].x-ctrls[ctrls.size()-2].x, 2.0*ctrls[ctrls.size()-1].y-ctrls[ctrls.size()-2].y);
	ctrlPts.insert(ctrlPts.begin() + 0, pha);
	ctrlPts.push_back(phb);

	double mudelta = 1.0 / resolution;
	int count = ctrlPts.size();
	int point = 0;

	curve.resize((ctrlPts.size() - 3) * resolution + 1);

	for(int n = 3; n < count; n++) {
		for(int j = 0; j < resolution; j++) {
			float mu = mudelta * (float)j;
			float x = bSpline(mu, ctrlPts[n-3].x, ctrlPts[n-2].x, ctrlPts[n-1].x, ctrlPts[n].x);
			float y = bSpline(mu, ctrlPts[n-3].y, ctrlPts[n-2].y, ctrlPts[n-1].y, ctrlPts[n].y);
			curve[point].set(x, y);
			point += 1;
		}
		if(n == count-1) {
			float x = bSpline(1.0, ctrlPts[n-3].x, ctrlPts[n-2].x, ctrlPts[n-1].x, ctrlPts[n].x);
			float y = bSpline(1.0, ctrlPts[n-3].y, ctrlPts[n-2].y, ctrlPts[n-1].y, ctrlPts[n].y);
			curve[point].set(x, y);
			point += 1;
		}
	}
}

//---------------------------------------------------------------------------------------------------------------
inline float catmullRom( float u, float x0, float x1, float x2, float x3 ) {
	float u2 = u * u;
	float u3 = u2 * u;
	return ((2.0 * x1) + 
		   (-x0 + x2) * u + 
		   (2.0*x0 - 5.0*x1 + 4.0*x2 - x3) * u2 + 
		   (-x0 + 3.0*x1 - 3.0*x2 + x3) * u3) * 0.5f;
}

inline void getCatmullCurvePoints(vector<ofVec2f> & ctrls, vector<ofVec2f> & curve, int resolution = 2) {

	if(ctrls.size() < 2) return;

	vector<ofVec2f> ctrlPts = ctrls;

	// Phantom Points
	ofVec2f pha, phb;
	pha.set(2.0*ctrls[0].x-ctrls[1].x, 2.0*ctrls[0].y-ctrls[1].y);
	phb.set(2.0*ctrls[ctrls.size()-1].x-ctrls[ctrls.size()-2].x, 2.0*ctrls[ctrls.size()-1].y-ctrls[ctrls.size()-2].y);
	ctrlPts.insert(ctrlPts.begin() + 0, pha);
	ctrlPts.push_back(phb);

	double mudelta = 1.0 / resolution;
	int count = ctrlPts.size();
	int point = 0;

	curve.resize((ctrlPts.size() - 3) * resolution + 1);

	for(int n = 3; n < count; n++) {
		for(int j = 0; j < resolution; j++) {
			float mu = mudelta * (float)j;
			float x = catmullRom(mu, ctrlPts[n-3].x, ctrlPts[n-2].x, ctrlPts[n-1].x, ctrlPts[n].x);
			float y = catmullRom(mu, ctrlPts[n-3].y, ctrlPts[n-2].y, ctrlPts[n-1].y, ctrlPts[n].y);
			curve[point].set(x, y);
			point += 1;
		}
		if(n == count-1) {
			float x = catmullRom(1.0, ctrlPts[n-3].x, ctrlPts[n-2].x, ctrlPts[n-1].x, ctrlPts[n].x);
			float y = catmullRom(1.0, ctrlPts[n-3].y, ctrlPts[n-2].y, ctrlPts[n-1].y, ctrlPts[n].y);
			curve[point].set(x, y);
			point += 1;
		}
	}
}

//---------------------------------------------------------------------------------------------------------------
inline void getCurvePoints(vector<ofVec2f> & ctrls, vector<ofVec2f> & curve, int resolution = 2, phdCurveInterpolatorMode mode = pcimCatmullRom) {

	if(ctrls.size() < 2) return;

	if(mode == pcimLinear) {
		getInterpolatedPoints(ctrls, curve, resolution);
	} else if(mode == pcimBSpline) {
		getBSplineCurvePoints(ctrls, curve, resolution);
	} else if(mode == pcimCatmullRom) {
		getCatmullCurvePoints(ctrls, curve, resolution);
	}
}

//---------------------------------------------------------------------------------------------------------------
class phdQuad {
protected:
	float A, B, D, E, G, H;
	bool boundsUpdated, solveUpdated;
	ofRectangle bounds;

	void solvePerspective() {

		if(solveUpdated) return;

		float T = (pts[2][0] - pts[1][0]) * (pts[2][1] - pts[3][1]) - (pts[2][0] - pts[3][0]) * (pts[2][1] - pts[1][1]);

        G = ((pts[2][0] - pts[0][0]) * (pts[2][1] - pts[3][1]) - (pts[2][0] - pts[3][0]) * (pts[2][1] - pts[0][1])) / T;
        H = ((pts[2][0] - pts[1][0]) * (pts[2][1] - pts[0][1]) - (pts[2][0] - pts[0][0]) * (pts[2][1] - pts[1][1])) / T;

        A = G * (pts[1][0] - pts[0][0]);
        D = G * (pts[1][1] - pts[0][1]);
        B = H * (pts[3][0] - pts[0][0]);
        E = H * (pts[3][1] - pts[0][1]);

        G -= 1.0;
        H -= 1.0;

		solveUpdated = true;
	}

public:
	double pts[4][2];

	phdQuad(float _x = 0.0, float _y = 0.0, float _w = 1.0, float _h = 1.0) {
		pts[0][0] = _x;		pts[0][1] = _y;
		pts[1][0] = _x+_w;	pts[1][1] = _y;
		pts[2][0] = _x+_w;	pts[2][1] = _y+_w;
		pts[3][0] = _x;		pts[3][1] = _y+_w;
		solveUpdated = boundsUpdated = false;
	}

	phdQuad(float _x0, float _y0, float _x1, float _y1, float _x2, float _y2, float _x3, float _y3) {
		pts[0][0] = _x0; pts[0][1] = _y0;
		pts[1][0] = _x1; pts[1][1] = _y1;
		pts[2][0] = _x2; pts[2][1] = _y2;
		pts[3][0] = _x3; pts[3][1] = _y3;
		solveUpdated = boundsUpdated = false;
	}

	void operator=(phdQuad & _value) {
		setVertex(0, _value.getX(0), _value.getY(0));
		setVertex(1, _value.getX(1), _value.getY(1));
		setVertex(2, _value.getX(2), _value.getY(2));
		setVertex(3, _value.getX(3), _value.getY(3));
	}

	bool operator==(phdQuad & _value) {
		return
		pts[0][0] == _value.pts[0][0] && pts[0][1] == _value.pts[0][1] &&
		pts[1][0] == _value.pts[1][0] && pts[1][1] == _value.pts[1][1] &&
		pts[2][0] == _value.pts[2][0] && pts[2][1] == _value.pts[2][1] &&
		pts[3][0] == _value.pts[3][0] && pts[3][1] == _value.pts[3][1];
	}

	bool operator!=(phdQuad & _value) { return !(_value == *this); }

	inline bool isChanged() { return !(solveUpdated); }

	inline double getX(unsigned int index) { return index < 4 ? pts[index][0] : -1; }
	inline double getY(unsigned int index) { return index < 4 ? pts[index][1] : -1; }

	inline ofVec2f getVertex(unsigned int index) {
		return (index < 4) ? ofVec2f(pts[index][0], pts[index][1]) : ofVec2f(-1,-1);
	}

	bool setVertex(unsigned int index, double x, double y) {
		if(index < 4) {
			if(pts[index][0] != x || pts[index][1] != y) {
				pts[index][0] = x;
				pts[index][1] = y;
				solveUpdated = boundsUpdated = false;
				return true;
			}
		}
		return false;
	}
	inline bool setVertex(unsigned int index, ofVec2f & vertex) { return setVertex(index, vertex.x, vertex.y); }

	inline void setVertices(float _x, float _y, float _w, float _h) {
		pts[0][0] = _x;		pts[0][1] = _y;
		pts[1][0] = _x+_w;	pts[1][1] = _y;
		pts[2][0] = _x+_w;	pts[2][1] = _y+_h;
		pts[3][0] = _x;		pts[3][1] = _y+_h;
		solveUpdated = boundsUpdated = false;
	}

	inline void setVertices(float _x0, float _y0, float _x1, float _y1, float _x2, float _y2, float _x3, float _y3) {
		pts[0][0] = _x0; pts[0][1] = _y0;
		pts[1][0] = _x1; pts[1][1] = _y1;
		pts[2][0] = _x2; pts[2][1] = _y2;
		pts[3][0] = _x3; pts[3][1] = _y3;
		solveUpdated = boundsUpdated = false;
	}

	inline void getVertices(float _vertices[8]) {
		_vertices[0] = pts[0][0]; _vertices[1] = pts[0][1];
		_vertices[2] = pts[1][0]; _vertices[3] = pts[1][1];
		_vertices[4] = pts[2][0]; _vertices[5] = pts[2][1];
		_vertices[6] = pts[3][0]; _vertices[7] = pts[3][1];
	}

	ofVec2f getPerspectiveVertex(float u, float v) {

		solvePerspective();

		float T = G * u + H * v + 1.0;
		return ofVec2f((A * u + B * v) / T + pts[0][0], (D * u + E * v) / T + pts[0][1]);
	}

	ofVec2f getBilinearVertex(float u, float v) {

		ofVec2f vr, vs;

        vr.x = (1.0 - u) * pts[0][0] + u * pts[1][0];
        vr.y = (1.0 - u) * pts[0][1] + u * pts[1][1];
        vs.x = (1.0 - u) * pts[3][0] + u * pts[2][0];
        vs.y = (1.0 - u) * pts[3][1] + u * pts[2][1];

        return ofVec2f((1.0 - v) * vr.x + v * vs.x, (1.0 - v) * vr.y + v * vs.y);
	}

	ofRectangle & getBounds() {

		if(boundsUpdated) return bounds;

		ofVec2f tl = ofVec2f(pts[0][0], pts[0][1]);
		ofVec2f br = ofVec2f(pts[0][0], pts[0][1]);
		for(int i = 1; i < 4; i++) {
			if(tl.x > pts[i][0]) tl.x = pts[i][0];
			if(tl.y > pts[i][1]) tl.y = pts[i][1];
			if(br.x < pts[i][0]) br.x = pts[i][0];
			if(br.y < pts[i][0]) br.y = pts[i][1];
		}
		bounds.set(tl.x, tl.y, br.x-tl.x, br.y-tl.y);

		boundsUpdated = true;

		return bounds;
	}

	void draw() {
		glBegin(GL_LINES);
			glVertex2f(pts[0][0], pts[0][1]); glVertex2f(pts[1][0], pts[1][1]);
			glVertex2f(pts[1][0], pts[1][1]); glVertex2f(pts[2][0], pts[2][1]);
			glVertex2f(pts[2][0], pts[2][1]); glVertex2f(pts[3][0], pts[3][1]);
			glVertex2f(pts[3][0], pts[3][1]); glVertex2f(pts[0][0], pts[0][1]);
		glEnd();
	}
};
//---------------------------------------------------------------------------------------------------------------
enum phdWarpFaceType { wftQuadDLD, wftQuadDLU, wftTriTL, wftTriTR, wftTriBR, wftTriBL };

//---------------------------------------------------------------------------------------------------------------
class phdWarper {
protected:
	bool changed;
	phdQuad quadIn;
	phdQuad quadOut;
	phdWarpFaceType faceType;

public:
	GLfloat	glWarpMatrix[16];
	GLfloat glInvWarpMatrix[16];

	phdWarper() {
		makeIdentity(glWarpMatrix);
		makeIdentity(glInvWarpMatrix);
		faceType = wftQuadDLD;
		changed = true;
	}

	inline void setInputVertex(unsigned int vertex, float x, float y) { 
		if(quadIn.setVertex(vertex, x, y)) changed = true;
	}
	inline void setOutputVertex(unsigned int vertex, float x, float y) {
		if(quadOut.setVertex(vertex, x, y)) changed = true;
	}

	inline phdQuad & getQuadIn()  { return quadIn; }
	inline void setQuadIn(phdQuad & _value) {
		setInputVertex(0, _value.getX(0), _value.getY(0));
		setInputVertex(1, _value.getX(1), _value.getY(1));
		setInputVertex(2, _value.getX(2), _value.getY(2));
		setInputVertex(3, _value.getX(3), _value.getY(3));
	}

	inline phdQuad & getQuadOut() { return quadOut; }
	inline void setQuadOut(phdQuad & _value) {
		setOutputVertex(0, _value.getX(0), _value.getY(0));
		setOutputVertex(1, _value.getX(1), _value.getY(1));
		setOutputVertex(2, _value.getX(2), _value.getY(2));
		setOutputVertex(3, _value.getX(3), _value.getY(3));
	}

	void updateWarpMatrix() {

		if(!changed) return; // dont need any update

		calcWarpMatrix(quadIn.pts, quadOut.pts, glWarpMatrix);

		calcWarpMatrix(quadOut.pts, quadIn.pts, glInvWarpMatrix);

		changed = false;
	}

	inline void begin() {
		updateWarpMatrix();
		glPushMatrix();
		glMultMatrixf(glWarpMatrix);
	}
	inline void end() { glPopMatrix(); }

	inline void invBegin() {
		updateWarpMatrix();
		glPushMatrix();
		glMultMatrixf(glInvWarpMatrix);
	}
	inline void invEnd() { glPopMatrix(); }
};

//--------------------------------------------------------------------------------------------------------------
enum phdCropAreaMode { pcrmNone, pcrmRect, pcrmQuad };

//--------------------------------------------------------------------------------------------------------------
class phdCropArea {
protected:
	phdCropAreaMode cropMode;
	ofRectangle cropRect;
	phdQuad cropQuad;

public:
	phdCropArea() {
		cropMode = pcrmNone;
		cropRect.set(0.0,0.0,1.0,1.0);
		cropQuad.setVertices(0.0,0.0,1.0,1.0);
	}
	virtual ~phdCropArea() { }

	inline void setCropArea(float _x, float _y, float _w, float _h) {
		cropRect.set(_x,_y,_w,_h);
		cropRect.standardize();
		if(cropMode == pcrmQuad) cropQuad.setVertices(cropRect.x, cropRect.y, cropRect.width, cropRect.height);
	}

	inline void setCropArea(float _x0, float _y0, float _x1, float _y1, float _x2, float _y2, float _x3, float _y3) {
		cropQuad.setVertices(_x0,_y0,_x1,_y1,_x2,_y2,_x3,_y3);
		if(cropMode == pcrmRect) {
			cropRect.set(_x0, _y0, _x2-_x0, _y2-_y0);
			cropRect.standardize();
		}
	}

	inline void fillQuadFromCropResizedByWH(phdQuad & _quad, float _w, float _h) {
		if(cropMode == pcrmNone) {
			_quad.setVertices(0.0, 0.0, _w, _h);
		} else if(cropMode == pcrmRect) {
			_quad.setVertices(cropRect.x * _w, cropRect.y * _h, (cropRect.x + cropRect.width) * _w, (cropRect.y + cropRect.height) * _h);
		} else if(cropMode == pcrmQuad) {
			_quad.setVertex(0, cropQuad.getX(0) * _w, cropQuad.getY(0) * _h);
			_quad.setVertex(1, cropQuad.getX(1) * _w, cropQuad.getY(1) * _h);
			_quad.setVertex(2, cropQuad.getX(2) * _w, cropQuad.getY(2) * _h);
			_quad.setVertex(3, cropQuad.getX(3) * _w, cropQuad.getY(3) * _h);
		}
	}
};
