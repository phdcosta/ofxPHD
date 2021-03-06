#pragma once

#include "phdUtils.h"

void printDebug(string & _str) { printf((_str + "\n").c_str()); }

ofVec2f rotatePoint(float x, float y, float cx, float cy, float ccx, float ccy, float sine, float cosine) {
	float ax = x - cx;
	float ay = y - cy;
	return ofVec2f(ccx + ax*cosine + ay*sine, ccy + ax*-sine + ay*cosine);
}

void circlesTangent(float _x1, float _y1, float _r1, float _x2, float _y2, float _r2, vector<ofVec2f> & _tangents) {
	double x1 = 200.0;
	double y1 = 200.0;
	double r1 = 50;

	double x2 = 250;
	double y2 = 200;
	double r2 = 30;

	double d_sq = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
	if (d_sq <= (r1-r2)*(r1-r2)) return;
 
	double d = sqrt(d_sq);
	double vx = (x2 - x1) / d;
	double vy = (y2 - y1) / d;
 
	double res[4][4];
	int i = 0;

	ofNoFill();
	ofSetColor(255,255,255,255);
	ofCircle(x1,y1,r1);
	ofCircle(x2,y2,r2);
	ofFill();

	for (int sign1 = +1; sign1 >= -1; sign1 -= 2) {
		double c = (r1 - sign1 * r2) / d;
		if(c*c > 1.0) continue;
		double h = sqrt(max(0.0, 1.0 - c*c));
		for (int sign2 = +1; sign2 >= -1; sign2 -= 2) {
			double nx = vx * c - sign2 * h * vy;
			double ny = vy * c + sign2 * h * vx;
			double a[4];
			a[0] = x1 + r1 * nx;
			a[1] = y1 + r1 * ny;
			a[2] = x2 + sign1 * r2 * nx;
			a[3] = y2 + sign1 * r2 * ny;

			_tangents.push_back(ofVec2f(a[0], a[1]));
			_tangents.push_back(ofVec2f(a[2], a[3]));
		}
	}
}

void circlesIntersection(double x0, double y0, double r0, double x1, double y1, double r1, vector<ofVec2f> & _points) {

	double xi, yi, xi_prime, yi_prime;
	double a, dx, dy, d, h, rx, ry;
	double x2, y2;

	dx = x1 - x0;
	dy = y1 - y0;

	d = hypot(dx,dy);

	if (d > (r0 + r1)) { return; }
	if (d < fabs(r0 - r1)) { return; }

	a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

	x2 = x0 + (dx * a/d);
	y2 = y0 + (dy * a/d);

	h = sqrt((r0*r0) - (a*a));

	rx = -dy * (h/d);
	ry = dx * (h/d);

	xi = x2 + rx;
	xi_prime = x2 - rx;
	yi = y2 + ry;
	yi_prime = y2 - ry;

	_points.push_back(ofVec2f(xi, yi));
	_points.push_back(ofVec2f(xi_prime, yi_prime));
}

bool getTriangleFromSegmentSize(double _szA, double _szB, double _szC, double x0, double y0, double _angle, double _scale, ofVec2f & _ptA, ofVec2f & _ptB, ofVec2f & _ptC) {

	double r0 = _szB * _scale;
	double x1 = 0 + _szA * _scale;
	double y1 = 0;
	double r1 = _szC * _scale;

	vector<ofVec2f> _points;

	circlesIntersection(0, 0, r0, x1, y1, r1, _points);

	if(_points.size() == 0) return false; // dont intersect --- sides dont make a triangle
	
	_ptA.set(0, 0);
	_ptB.set(x1, y1);
	_ptC.set(_points[1].x, _points[1].y);

	_ptC.x = _ptB.x + (_ptA.x - _ptC.x);
	_ptC.y = _ptB.y + (_ptA.y - _ptC.y);

	double _sin = sin(ofDegToRad(_angle));
	double _cos = cos(ofDegToRad(_angle));

	ofVec2f _pivot = (_ptA + _ptB + _ptC) / 3.0;
	
	_ptA -= _pivot;
	_ptB -= _pivot;
	_ptC -= _pivot;

	_ptA = rotatePoint(_ptA.x, _ptA.y, 0, 0, 0, 0, _sin, _cos);
	_ptB = rotatePoint(_ptB.x, _ptB.y, 0, 0, 0, 0, _sin, _cos);
	_ptC = rotatePoint(_ptC.x, _ptC.y, 0, 0, 0, 0, _sin, _cos);

	_ptA += ofVec2f(x0,y0);
	_ptB += ofVec2f(x0,y0);
	_ptC += ofVec2f(x0,y0);

	return true;
}

// compute the area of a triangle using Heron's formula
double triareaHeron(double a, double b, double c)
{
    double s = (a + b + c)/2.0;
    double area=sqrt(fabs(s*(s-a)*(s-b)*(s-c)));
    return area;     
}

// compute the distance between two 2d points
double dist(double x0, double y0, double z0, double x1, double y1, double z1)
{
    double a = x1 - x0;   
    double b = y1 - y0;
    double c = z1 - z0;
    return sqrt(a*a + b*b + c*c);
}

// calculate barycentric coordinates
// triangle 1st vertex: x0,y0,z0
// triangle 2nd vertex: x1,y1,z1
//  triangle 3rd vertex: x2,y2,z2
// point inside triangle: vx, vy,vz
// *u,*v,*w are the coordinates returned
void barycent(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2,
                         double vx, double vy, double vz,
                         double *u, double *v, double *w)
{
    // compute the area of the big triangle
    double a = dist(x0, y0, z0, x1, y1, z1);
    double b = dist(x1, y1, z1, x2, y2, z2);
    double c = dist(x2, y2, z2, x0, y0, z0);
    double totalarea = triareaHeron(a, b, c);
        
    // compute the distances from the outer vertices to the inner vertex
    double length0 = dist(x0, y0, z0, vx, vy, vz);        
    double length1 = dist(x1, y1, z1, vx, vy, vz);        
    double length2 = dist(x2, y2, z2, vx, vy, vz);        
    
    // divide the area of each small triangle by the area of the big triangle
    *u = triareaHeron(b, length1, length2)/totalarea;
    *v = triareaHeron(c, length0, length2)/totalarea;
    *w = triareaHeron(a, length0, length1)/totalarea;          
}

void setIdentityMatrix(ofMatrix4x4 & _m) {
	_m(0,0) = 1.0; _m(0,1) = 0.0; _m(0,2) = 0.0; _m(0,3) = 0.0;
	_m(1,0) = 0.0; _m(1,1) = 1.0; _m(1,2) = 0.0; _m(1,3) = 0.0;
	_m(2,0) = 0.0; _m(2,1) = 0.0; _m(2,2) = 1.0; _m(2,3) = 0.0;
	_m(3,0) = 0.0; _m(3,1) = 0.0; _m(3,2) = 0.0; _m(3,3) = 1.0;
};

void setRotateMatrixX(ofMatrix4x4 & _m, double angle) {
	double sine = sin(angle);
	double cosine = cos(angle);
	setIdentityMatrix(_m);
	_m(1,1) = cosine; _m(1,2) = -sine;
	_m(2,1) = sine;   _m(2,2) = cosine;
};

void setRotateMatrixY(ofMatrix4x4 & _m, double angle) {
	double sine = sin(angle);
	double cosine = cos(angle);
	setIdentityMatrix(_m);
	_m(0,0) = cosine; _m(0,2) = sine;
	_m(2,0) = -sine;  _m(2,2) = cosine;
};

void setRotateMatrixZ(ofMatrix4x4 & _m, double angle) {
	double sine = sin(angle);
	double cosine = cos(angle);
	setIdentityMatrix(_m);
	_m(0,0) = cosine; _m(0,1) = -sine;
	_m(1,0) = sine;   _m(1,1) = cosine;
};

void setRotateMatrixXYZ(ofMatrix4x4 & _m, double angleX, double angleY, double angleZ) {
	ofMatrix4x4 _mx; setRotateMatrixX(_mx, angleX);
	ofMatrix4x4 _my; setRotateMatrixX(_my, angleY);
	ofMatrix4x4 _mz; setRotateMatrixX(_mz, angleZ);
	_m = _mx * _my * _mz;
};

void setScaleMatrix(ofMatrix4x4 & _m, double scaleX, double scaleY, double scaleZ) {
	setIdentityMatrix(_m);
	_m(0,0) = scaleX;
	_m(1,1) = scaleY;
	_m(2,2) = scaleZ;
};

void setTranslateMatrix(ofMatrix4x4 & _m, double moveX, double moveY, double moveZ) {
	setIdentityMatrix(_m);
	_m(3,0) = moveX;
	_m(3,1) = moveY;
	_m(3,2) = moveZ;
};

void setTransformationMatrix(ofMatrix4x4 & _m, ofVec2f & _curPos, float _curAngle, ofVec2f & _curScale, ofVec2f & _newPos, float _newAngle, ofVec2f & _newScale) {

	ofMatrix4x4 _aux;

	_m.makeTranslationMatrix(-_curPos);
	
	_aux.makeRotationMatrix(-_curAngle, 0.0, 0.0, 1.0);
	_m = _m * _aux;

	_aux.makeScaleMatrix(_newScale.x / _curScale.x, _newScale.y / _curScale.y, 1.0);
	_m = _m * _aux;

	_aux.makeRotationMatrix(_newAngle, 0.0, 0.0, 1.0);
	_m = _m * _aux;

	_aux.makeTranslationMatrix(_newPos);
	_m = _m * _aux;
}

//--------------------------------------------------------------
// Code by Theo from URL above
//--------------------------------------------------------------
// http://www.openframeworks.cc/forum/viewtopic.php?f=9&t=1443
//--------------------------------------------------------------
bool phdPointInsidePolygon(vector<ofVec2f> & _points, float px, float py) {
	int N = _points.size();
	if (N == 0) return false;
	int counter = 0;
	int i;
	double xinters;
	ofVec2f p1, p2;

	p1 = _points[0];
	for (int i = 1; i <= N; i++) {
		p2 = _points[i % N];
		if (py > MIN(p1.y,p2.y)) {
			if (py <= MAX(p1.y,p2.y)) {
				if (px <= MAX(p1.x,p2.x)) {
					if (p1.y != p2.y) {
						xinters = (py-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
						if (p1.x == p2.x || px <= xinters) counter++;
					}
				}
			}
		}
		p1 = p2;
	}

	if (counter % 2 == 0)
		return false;
	else
		return true;
}

bool phdPointInsideTriangle(float _x, float _y, float _xa, float _ya, float _xb, float _yb, float _xc, float _yc) {

    float s = _ya * _xc - _xa * _yc + (_yc - _ya) * _x + (_xa - _xc) * _y;
    float t = _xa * _yb - _ya * _xb + (_ya - _yb) * _x + (_xb - _xa) * _y;

    if ((s < 0) != (t < 0)) return false;

    float A = -_yb * _xc + _ya * (_xc - _xb) + _xa * (_yb - _yc) + _xb * _yc;

    if (A < 0.0) { s = -s; t = -t; A = -A; }

    return s > 0 && t > 0 && (s + t) < A;
}

double phdDistancePointLineSquared(double x, double y, double x1, double y1, double x2, double y2) {
	double A = x - x1;
	double B = y - y1;
	double C = x2 - x1;
	double D = y2 - y1;
	double dot = A * C + B * D;
	double len_sq = C * C + D * D;
	double param = dot / len_sq;
	double xx,yy;
	if(param < 0) {
		xx = x1;
		yy = y1;
	} else if(param > 1) {
		xx = x2;
		yy = y2;
	} else {
		xx = x1 + param * C;
		yy = y1 + param * D;
	}
	return ((x - xx) * (x - xx)) + ((y - yy) * (y - yy));
}

bool phdPointNearEdge(double x1, double y1, double x2, double y2, double px, double py, double _dist) {
	return (phdDistancePointLineSquared(px,py,x1,y1,x2,y2) < _dist*_dist);
}

int phdPointNearEdge(vector<ofVec2f> & _points, double px, double py, double _dist, bool _closed) {

	int nPoints = _points.size();

	if(nPoints < 2) return -1;

	for(int i = 0; i < nPoints-1; i++) {
		if(phdDistancePointLineSquared(px, py, _points[i].x, _points[i].y, _points[i+1].x, _points[i+1].y) < _dist*_dist) return i;
	}
	if(_closed && nPoints > 2) {
		if(phdDistancePointLineSquared(px, py, _points[0].x, _points[0].y, _points[nPoints-1].x, _points[nPoints-1].y) < _dist*_dist) return nPoints-1;
	}

	return -1;
}

int phdPointOverVertex(vector<ofVec2f> & _points, double px, double py) {
	for(int i = 0; i < _points.size(); i++) {
		if(fabs(_points[i].x - px) < 5 && fabs(_points[i].y - py) < 5) {
			return i;
		}
	}
	return -1;
}

// calculates centroid and return polygon area
double phdGetCentroid(vector<ofVec2f> & _points, ofVec2f & _centroid) {

	int i,j;
	double ai, atmp, xtmp, ytmp;

	int _sz = _points.size();

	if(_sz == 0) { // no points
		_centroid.set(0.0,0.0);
		return 0.0;
	}

	if(_sz == 1) { // just one point
		_centroid = _points[0];
		return 0.0;
	}

	if(_sz == 2) { // a line
		_centroid.set((_points[0].x + _points[1].x) * 0.5, (_points[0].y + _points[1].y) * 0.5);
		return 0.0;
	}

	atmp = xtmp = ytmp = 0.0;

	i = _sz-1;

	for(int j = 0; j < _sz; j++) {
		ai = _points[i].x * _points[j].y - _points[j].x * _points[i].y;
		atmp = atmp + ai;
		xtmp = xtmp + (_points[j].x + _points[i].x) * ai;
		ytmp = ytmp + (_points[j].y + _points[i].y) * ai;
		i = j;
	}

	if(atmp != 0.0) {
		_centroid.x = xtmp / (3.0 * atmp);
		_centroid.y = ytmp / (3.0 * atmp);
	}

	return atmp * 0.5;
}

void getTransformedPolyline(vector<ofVec2f> & _src, vector<ofVec2f> & _dst, ofMatrix4x4 _mat) {
	_dst = _src;
	for(int i = 0; i < _dst.size(); i++) {
		_dst[i] = ofVec3f(_dst[i]) * _mat;
	}
}

// 2D cross product.
// Return a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
static CoordType cross(const phdConvexPoint &O, const phdConvexPoint &A, const phdConvexPoint &B)
{
	return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

void contourToConvexHull(vector<ofVec2f> &src, vector<ofVec2f> &dst) {

	dst.clear();

	vector<phdConvexPoint> P(src.size());
	for(int i = 0; i < src.size(); i++) {
		P[i].x = src[i].x;
		P[i].y = src[i].y;
	}

	int n = src.size(), k = 0;
	vector<phdConvexPoint> H(2*n);
 
	// Sort points lexicographically
	sort(P.begin(), P.end());
 
	// Build lower hull
	for (int i = 0; i < n; i++) {
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
 
	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
 
	H.resize(k);

	for(int i = 0; i < H.size(); i++) { dst.push_back(ofVec2f(H[i].x + 500, H[i].y)); }
}

double triareaVertices(double _ax, double _ay, double _bx, double _by, double _cx, double _cy) {
	return fabs(_ax*(_by-_cy) + _bx*(_cy-_ay) + _cx*(_ay-_by)) / 2.0;
}

void drawEdgeWithSizeLabel(double _x1, double _y1, double _x2, double _y2, string _label, bool _edge) {
	double _cx = (_x1+_x2) / 2.0;
	double _cy = (_y1+_y2) / 2.0;
	double _a  = atan2(_y2-_y1,_x2-_x1);

	if(_a >= 0.0 && _a < PI*0.5) {  // x-  y+
		ofDrawBitmapString(_label, _cx - 80, _cy + 15);
	} else if(_a >= PI*0.5 && _a < PI) {  // x+  y-
		ofDrawBitmapString(_label, _cx - 80, _cy - 20);
	} else if(_a >= -PI && _a < -PI*0.5) {  // x-  y-
		ofDrawBitmapString(_label, _cx + 10, _cy - 20);
	} else if(_a >= -PI*0.5 && _a < 0.0) {  // x-  y+
		ofDrawBitmapString(_label, _cx + 10, _cy + 15);
	}
	if(_edge) ofLine(_x1,_y1,_x2,_y2);
}

void drawGrid(float _xGap, float _yGap, ofTrueTypeFont * _font) {

	ofSetLineWidth(3);
	ofSetColor(255,255,255,128);
	if(_font != NULL) {
		for(int x = 0; x <= 1024; x += _xGap * 5) {
			for(int y = 0; y <= 768; y += _yGap * 5) {
				if(_font != NULL) _font->drawString(ofToString(y/64,0)+":"+ofToString(x/64,0), x+16, y+40);
			}
		}
	}

	ofSetLineWidth(1.5);
	ofSetColor(255,255,255,96);
	for(int x = 0; x <= 1024; x += _xGap) { ofLine(x, 0, x, 768); }
	for(int y = 0; y <= 768;  y += _yGap) { ofLine(0, y, 1024, y); }

	ofSetColor(255,255,255,255);
	for(int x = 0; x <= 1024; x += _xGap * 5) { ofLine(x, 0, x, 768); }
	for(int y = 0; y <= 768;  y += _yGap * 5) { ofLine(0, y, 1024, y); }

	ofSetColor(255,255,255,255);
	ofSetLineWidth(3);
	ofLine(0,0,1024,0);
	ofLine(1024,0,1024,768);
	ofLine(1024,768,0,768);
	ofLine(0,768,0,0);
	ofSetLineWidth(1);
}

void drawFilledBorderRectangle(float _x1, float _y1, float _w, float _h, ofColor & _fill, ofColor & _border) {
	ofFill();
	ofSetColor(_fill);
	ofRect(_x1,_y1,_w,_h);
	ofSetColor(_border);
	glBegin(GL_LINES);
		glVertex2f(_x1,_y1);		glVertex2f(_x1+_w,_y1);
		glVertex2f(_x1+_w,_y1);		glVertex2f(_x1+_w,_y1+_h);
		glVertex2f(_x1+_w,_y1+_h);	glVertex2f(_x1,_y1+_h);
		glVertex2f(_x1,_y1+_h);		glVertex2f(_x1,_y1);
	glEnd();
}

void drawBorderRectangle(float _x1, float _y1, float _w, float _h, ofColor & _border) {
	ofSetColor(_border);
	glBegin(GL_LINES);
		glVertex2f(_x1,_y1);		glVertex2f(_x1+_w,_y1);
		glVertex2f(_x1+_w,_y1);		glVertex2f(_x1+_w,_y1+_h);
		glVertex2f(_x1+_w,_y1+_h);	glVertex2f(_x1,_y1+_h);
		glVertex2f(_x1,_y1+_h);		glVertex2f(_x1,_y1);
	glEnd();
}

void drawBorderRectangle(float _x1, float _y1, float _w, float _h) {
	glBegin(GL_LINES);
		glVertex2f(_x1,_y1);		glVertex2f(_x1+_w,_y1);
		glVertex2f(_x1+_w,_y1);		glVertex2f(_x1+_w,_y1+_h);
		glVertex2f(_x1+_w,_y1+_h);	glVertex2f(_x1,_y1+_h);
		glVertex2f(_x1,_y1+_h);		glVertex2f(_x1,_y1);
	glEnd();
}

//--------------------------------------------------------------------------------------------------------------
bool isClickEqualsDrag(ofMouseEventArgs & onClick, ofMouseEventArgs & onDrag) { // onClick == onDrag when was dragging something
	return (onClick.x == onDrag.x && onClick.y == onDrag.y && onClick.button == onDrag.button);
}


//--------------------------------------------------------------------------------------------------------------
// Draw a BezierSpline approximating a curve defined by the array of points.
// Beware! A B-Spline generaly does not normally pass through the points
// defining it !
//--------------------------------------------------------------------------------------------------------------
double Calculate(double mu, double p0, double p1, double p2, double p3) {
	double mu2, mu3;
	mu2 = mu * mu;
	mu3 = mu * mu2;
	return (double)((1.0/6.0)*(	mu3*(-p0+3.0*p1-3.0*p2+p3)+
								mu2*(3.0*p0-6.0*p1+3.0*p2)+
  								mu *(-3.0*p0+3.0*p2)+(p0+4.0*p1+p2)));
}

void bezierSplinePoints(vector<ofVec2f> pnts, int count, int segments, vector<ofVec2f> & points) {
	double mu, mudelta;
	int x1, y1, x2, y2, n, h;
	ofVec2f pha, phb;

	if(count < 4 || count > 16383) return;

	// Phantom Points
	pha = ofVec2f(2.0*pnts[0].x-pnts[1].x, 2.0*pnts[0].y-pnts[1].y);
	phb = ofVec2f(2.0*pnts[count-1].x-pnts[count-2].x, 2.0*pnts[count-1].y-pnts[count-2].y);

	mudelta = 1.0 / segments;

	for(n = 2; n < count; n++) {

		mu = 0;

		if(n == 2) {
			x1 = Calculate(mu,pha.x,pnts[n-2].x,pnts[n-1].x,pnts[n].x);
			y1 = Calculate(mu,pha.y,pnts[n-2].y,pnts[n-1].y,pnts[n].y);
		} else if(n == count) {
			x1 = Calculate(mu,pnts[n-3].x,pnts[n-2].x,pnts[n-1].x,phb.x);
			y1 = Calculate(mu,pnts[n-3].y,pnts[n-2].y,pnts[n-1].y,phb.y);
		} else {
			x1 = Calculate(mu,pnts[n-3].x,pnts[n-2].x,pnts[n-1].x,pnts[n].x);
			y1 = Calculate(mu,pnts[n-3].y,pnts[n-2].y,pnts[n-1].y,pnts[n].y);
		}

		points.push_back(ofVec2f(x1, y1));

		mu = mu + mudelta;

		for(h = 1; h < segments; h++) {

			if(n == 2) {
				x2 = Calculate(mu,pha.x,pnts[n-2].x,pnts[n-1].x,pnts[n].x);
				y2 = Calculate(mu,pha.y,pnts[n-2].y,pnts[n-1].y,pnts[n].y);
			} else if(n == count) {
				x2 = Calculate(mu,pnts[n-3].x,pnts[n-2].x,pnts[n-1].x,phb.x);
				y2 = Calculate(mu,pnts[n-3].y,pnts[n-2].y,pnts[n-1].y,phb.y);
			} else {
				x2 = Calculate(mu,pnts[n-3].x,pnts[n-2].x,pnts[n-1].x,pnts[n].x);
				y2 = Calculate(mu,pnts[n-3].y,pnts[n-2].y,pnts[n-1].y,pnts[n].y);
			}

			points.push_back(ofVec2f(x2, y2));

			mu = mu + mudelta;
		}
	}
}


/*		float dx = _x - mx;
		float dy = _y - my;

		float _x1 = mx;
		float _y1 = my;
		float _x2, _y2;

		ofSetColor(255,255,255,255);

		glBegin(GL_LINES);
		glVertex2d(_x1,_y1);
		for(int i = 1; i <= 10; i++) {
			_x2 = mx + fabs(mapEaseQuartic(i, 10, distance));
			_y2 = my + mapLinear(i, 10, distance);
			glVertex2d(_x2,_y2);
			_x1 = _x2; _y1 = _y2;
			if(i < 10) glVertex2d(_x1,_y1);
		}
		glEnd();

		_x1 = mx; _y1 = my;

		glBegin(GL_LINES);
		glVertex2d(_x1,_y1);
		for(int i = 1; i <= 10; i++) {
			_x2 = mx + fabs(mapEaseQuartic(i, 10, distance));
			_y2 = my - mapLinear(i, 10, distance);
			glVertex2d(_x2,_y2);
			_x1 = _x2; _y1 = _y2;
			if(i < 10) glVertex2d(_x1,_y1);
		}
		glEnd();
		ofCircle(mx + fabs(mapEaseQuartic(dy, distance, distance)), my + mapLinear(dy, distance, distance), 5);
*/

//--------------------------------------------------------------------------------------------------------------
float getBestHeightbyFixedW(float _width, float _height, float _minW) {
	return _minW / _width * _height;
}

//--------------------------------------------------------------------------------------------------------------
float getMinScale(float _dstW, float _dstH, float _srcW, float _srcH) {
	return MIN(_dstW / _srcW, _dstH / _srcH);
}

//--------------------------------------------------------------------------------------------------------------
void getBestSizeWH(float _dstW, float _dstH, float _srcW, float _srcH, float & _scaledW, float & _scaledH) {

	float _sw = _dstW / _srcW;
	float _sh = _dstH / _srcH;

	if(_sw < _sh) {
		_scaledW = _srcW * _sw;
		_scaledH = _srcH * _sw;
	} else {
		_scaledW = _srcW * _sh;
		_scaledH = _srcH * _sh;
	}
}

//--------------------------------------------------------------------------------------------------------------
void getStrecthModeArea(ofRectangle & area, float srcW, float srcH, float dstW, float dstH, phdStretchMode stretchMode) {
	if(stretchMode == psmNone) { // no resize --> at topleft
		area.x		= 0.0;
		area.y		= 0.0;
		area.width	= srcW;
		area.height	= srcH;
	} else if(stretchMode == psmCentered) { // no resize --> at center
		float _ccx = dstW * 0.5;
		float _ccy = dstH * 0.5;
		area.x		= _ccx-srcW*0.5;
		area.y		= _ccy-srcH*0.5;
		area.width	= srcW;
		area.height	= srcH;
	} else if(stretchMode == psmStretched) { // resized to destination size
		area.x		= 0.0;
		area.y		= 0.0;
		area.width	= dstW;
		area.height	= dstH;
	} else if(stretchMode == psmFitted) { // resized to fit proportionally to destination size at topleft
		float _scale = getMinScale(dstW, dstH, srcW, srcH);
		area.x		= 0.0;
		area.y		= 0.0;
		area.width	= srcW*_scale;
		area.height	= srcH*_scale;
	} else if(stretchMode == psmFitCentered) { // resized to fit proportionally to destination size at center
		float _ccx = dstW * 0.5;
		float _ccy = dstH * 0.5;
		float _scale = getMinScale(dstW, dstH, srcW, srcH);
		area.x		= _ccx-srcW*0.5*_scale;
		area.y		= _ccy-srcH*0.5*_scale;
		area.width	= srcW*_scale;
		area.height	= srcH*_scale;
	}
}

//--------------------------------------------------------------------------------------------------------------
ofRectangle getStrecthModeArea(float _x, float _y, float _w, float _h, float srcW, float srcH, phdStretchMode stretchMode) {
	ofRectangle area;
	if(stretchMode == psmNone) { // no resize --> at topleft
		area.x		= _x;
		area.y		= _y;
		area.width	= srcW;
		area.height	= srcH;
	} else if(stretchMode == psmCentered) { // no resize --> at center
		float _ccx = _w * 0.5;
		float _ccy = _h * 0.5;
		area.x		= _x + _ccx-srcW*0.5;
		area.y		= _y + _ccy-srcH*0.5;
		area.width	= srcW;
		area.height	= srcH;
	} else if(stretchMode == psmFitted) { // resized to fit proportionally to destination size at topleft
		float _scale = getMinScale(_w, _h, srcW, srcH);
		area.x		= _x;
		area.y		= _y;
		area.width	= srcW*_scale;
		area.height	= srcH*_scale;
	} else if(stretchMode == psmFitCentered) { // resized to fit proportionally to destination size at center
		float _ccx = _w * 0.5;
		float _ccy = _h * 0.5;
		float _scale = getMinScale(_w, _h, srcW, srcH);
		area.x		= _x + _ccx-srcW*0.5*_scale;
		area.y		= _y + _ccy-srcH*0.5*_scale;
		area.width	= srcW*_scale;
		area.height	= srcH*_scale;
	} else {//if(stretchMode == psmStretched) { // resized to destination size
		area.x		= _x;
		area.y		= _y;
		area.width	= _w;
		area.height	= _h;
	}
	return area;
}

//--------------------------------------------------------------------------------------------------------------
void drawBorderTriangle(float _x, float _y, float _w, float _h, phdTriangleDirection _dir) {

	float _h2 = _h / 2.0;
	float _w2 = _w / 2.0;

	if(_dir == ptdTop) {
		ofLine(_x + _w2, _y, _x + _w, _y + _h);
		ofLine(_x + _w, _y + _h, _x, _y + _h);
		ofLine(_x, _y + _h, _x + _w2, _y);
	} else if(_dir == ptdBottom) {
		ofLine(_x, _y, _x + _w, _y);
		ofLine(_x + _w, _y, _x + _w2, _y + _h);
		ofLine(_x + _w2, _y + _h, _x, _y);
	} else if(_dir == ptdLeft) {
		ofLine(_x + _w, _y, _x + _w, _y + _h);
		ofLine(_x + _w, _y + _h, _x, _y + _h2);
		ofLine(_x, _y + _h2, _x + _w, _y);
	} else if(_dir == ptdRight) {
		ofLine(_x, _y, _x + _w, _y + _h2);
		ofLine(_x + _w, _y + _h2, _x, _y + _h);
		ofLine(_x, _y + _h, _x, _y);
	}
}

//--------------------------------------------------------------------------------------------------------------
void drawFilledTriangle(float _x, float _y, float _w, float _h, phdTriangleDirection _dir) {

	float _h2 = _h / 2.0;
	float _w2 = _w / 2.0;

	if(_dir == ptdTop) {
		ofTriangle(_x + _w2, _y, _x + _w, _y + _h, _x, _y + _h);
	} else if(_dir == ptdBottom) {
		ofTriangle(_x, _y, _x + _w, _y, _x + _w2, _y + _h);
	} else if(_dir == ptdLeft) {
		ofTriangle(_x + _w, _y, _x + _w, _y + _h, _x, _y + _h2);
	} else if(_dir == ptdRight) {
		ofTriangle(_x, _y, _x + _w, _y + _h2, _x, _y + _h);
	}
}

//--------------------------------------------------------------------------------------------------------------
void genColorsFromBase(ofColor _base, ofColor & _normal, ofColor & _focused, ofColor & _selected, float _fNormal, float _fFocused, float _fSelected, bool _alpha) {
	if(_alpha) {
		_normal = ofColor(
			ofClamp(_base.r*_fNormal, 0, 255),
			ofClamp(_base.g*_fNormal, 0, 255),
			ofClamp(_base.b*_fNormal, 0, 255),
			ofClamp(_base.a*_fNormal, 0, 255));
		_focused = ofColor(
			ofClamp(_base.r*_fFocused, 0, 255),
			ofClamp(_base.g*_fFocused, 0, 255),
			ofClamp(_base.b*_fFocused, 0, 255),
			ofClamp(_base.a*_fFocused, 0, 255));
		_selected = ofColor(
			ofClamp(_base.r*_fSelected, 0, 255),
			ofClamp(_base.g*_fSelected, 0, 255),
			ofClamp(_base.b*_fSelected, 0, 255),
			ofClamp(_base.a*_fSelected, 0, 255));
	} else {
		_normal = ofColor(
			ofClamp(_base.r*_fNormal, 0, 255),
			ofClamp(_base.g*_fNormal, 0, 255),
			ofClamp(_base.b*_fNormal, 0, 255),
			_base.a);
		_focused = ofColor(
			ofClamp(_base.r*_fFocused, 0, 255),
			ofClamp(_base.g*_fFocused, 0, 255),
			ofClamp(_base.b*_fFocused, 0, 255),
			_base.a);
		_selected = ofColor(
			ofClamp(_base.r*_fSelected, 0, 255),
			ofClamp(_base.g*_fSelected, 0, 255),
			ofClamp(_base.b*_fSelected, 0, 255),
			_base.a);
	}
}

//--------------------------------------------------------------------------------------------------------------
void rotSliderV(float _x, float _y, float _w, float _h, float _value, ofColor _color, bool _filled) {
	ofPoint p(_x,_y,0.0);
	ofPath curve;
	curve.arc(p, _w, _h, 135, 135 + 270*_value);
	curve.arcNegative(p, _w*0.5, _h*0.5, 135 + 270*_value, 135);
	curve.close();
	curve.setArcResolution(36);
	curve.setFillColor(_color);
	curve.setFilled(_filled);
	curve.draw();
}

//--------------------------------------------------------------------------------------------------------------
void rotSliderV2(float _x, float _y, float _w, float _h, float _value, ofColor _clrA, ofColor _clrB, float _aIni, float _aSize) {
	ofPoint p(_x,_y,0.0);
	ofPath curve;
	if(_value != 0.0) {
		curve.arc(p, _w, _h, _aIni, _aIni + _aSize*_value);
		curve.arcNegative(p, _w*0.5, _h*0.5, _aIni + _aSize*_value, _aIni);
		curve.close();
		curve.setArcResolution(36);
		curve.setFillColor(_clrA);
		curve.setFilled(true);
		curve.draw();
	}
	if(_value != 1.0) {
		curve.clear();
		curve.arc(p, _w, _h, _aIni + _aSize*_value, _aIni+_aSize);
		curve.arcNegative(p, _w*0.5, _h*0.5, _aIni+_aSize, _aIni + _aSize*_value);
		curve.close();
		curve.setArcResolution(36);
		curve.setFillColor(_clrB);
		curve.setFilled(true);
		curve.draw();
	}
}

//--------------------------------------------------------------------------------------------------------------
double rotAngleToValue(float _angle, float _aIni, float _aSize, float _mark) {
	if(_angle < _mark) _angle += 360;
	return ofMap(_angle, _aIni, _aIni+_aSize, 0.0, 1.0, true);
}

//--------------------------------------------------------------------------------------------------------------
void rotSliderA(float _x, float _y, float _w, float _h, float _angle, ofColor _clrA, ofColor _clrB) {
	if(_angle < 90) _angle += 360;
	float _value = ofMap(_angle, 135, 405, 0.0, 1.0, true);
	rotSliderV2(_x,_y,_w,_h,_value,_clrA,_clrB);
}

//--------------------------------------------------------------------------------------------------------------
double atanPH(double _x, double _y, double _cx, double _cy) {
	double _a = atan2(_y-_cy,_x-_cx);
	if(_a < 0.0) return ofMap(_a, -3.14, 0.0, 180.0, 360.0, true);
	return ofMap(_a, 0.0, 3.14, 0.0, 180.0, true);
}

//--------------------------------------------------------------------------------------------------------------
float mapEaseQuartic(float value, float sizeIn, float sizeOut = 1.0) {
	float r = value / sizeIn; r = (r*r*r*r);
	if(value < 0.0) r = -r;
	return MIN(sizeOut, MAX(-sizeOut, r * sizeOut));
}

//--------------------------------------------------------------------------------------------------------------
void genTriangle(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices) {
	float cosine = cos(_angle);
	float sine = sin(_angle);
	float cx = _x+_w*0.5;
	float cy = _y+_h*0.5;
	_vertices.clear();
	_vertices.addVertex(rotatePoint(_x,    _y,        cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(_x+_w, _y+_h*0.5, cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(_x,    _y+_h,     cx, cy, cx, cy, sine, cosine));
	_vertices.setClosed(true);
}

void genDiamond(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices) {
	float cosine = cos(_angle);
	float sine = sin(_angle);
	float cx = _x+_w*0.5;
	float cy = _y+_h*0.5;
	_vertices.clear();
	_vertices.addVertex(rotatePoint(cx,			cy-_h*0.5,	cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(cx+_w*0.5,	cy,			cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(cx,			cy+_h*0.5,	cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(cx-_w*0.5,	cy,			cx, cy, cx, cy, sine, cosine));
	_vertices.setClosed(true);
}

void genQuad(float _x, float _y, float _w, float _h, float _angle, ofPolyline & _vertices) {
	float cosine = cos(_angle);
	float sine = sin(_angle);
	float cx = _x+_w*0.5;
	float cy = _y+_h*0.5;
	_vertices.clear();
	_vertices.addVertex(rotatePoint(_x,		_y,		cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(_x+_w,	_y,		cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(_x+_w,	_y+_h,	cx, cy, cx, cy, sine, cosine));
	_vertices.addVertex(rotatePoint(_x,		_y+_h,	cx, cy, cx, cy, sine, cosine));
	_vertices.setClosed(true);
}

//---------------------------------------------------------------------------------------------
// return projection v1 on to v2
ofVec2f projection(const ofVec2f & v1, const ofVec2f & v2) {
	float v2_ls = v2.lengthSquared();
	if( v2_ls > 0.00000000000001 ) return v2 * ( v2.dot(v1)/v2_ls );
	return ofVec2f(0.0, 0.0);
}

//---------------------------------------------------------------------------------------------
void phdMakeTransformationMatrix(ofMatrix4x4 & _mat, float _cx, float _cy, float _angFrom, float _angTo, float _dx, float _dy, float _sx, float _sy) {

	_mat.makeIdentityMatrix();

	ofMatrix4x4 _aux;
		
	_aux.makeTranslationMatrix(-_cx, -_cy, 0.0);			_mat = _mat * _aux;
	_aux.makeRotationMatrix(-_angFrom, 0.0, 0.0, 1.0);		_mat = _mat * _aux;
	_aux.makeScaleMatrix(_sx, _sy, 1.0);					_mat = _mat * _aux;
	_aux.makeRotationMatrix(_angTo, 0.0, 0.0, 1.0);			_mat = _mat * _aux;
	_aux.makeTranslationMatrix(_cx + _dx, _cy + _dy, 0.0);	_mat = _mat * _aux;
}

//---------------------------------------------------------------------------------------------
void getTransformationSinCos(float _aFrom, float _aTo, float & _sinF, float & _cosF, float & _sinT, float & _cosT) {
	_sinF = sinf(ofDegToRad(_aFrom));
	_cosF = cosf(ofDegToRad(_aFrom));
	_sinT = sinf(ofDegToRad(_aTo));
	_cosT = cosf(ofDegToRad(_aTo));
}

//---------------------------------------------------------------------------------------------
void transformPoint(ofVec2f & _pt, double _cx, double _cy, double _sinF, double _cosF, double _sinT, double _cosT, double _dx, double _dy, double _sx, double _sy) {

	double _x = _pt.x;
	double _y = _pt.y;

	double _u, _v;

	_x -= _cx; // translate to 0
	_y -= _cy;

	_u = _x * _cosF + _y * _sinF; // rotate to 0
	_v = _x *-_sinF + _y * _cosF;

	_x = _u * _sx; // scale
	_y = _v * _sy;

	_u = _x * _cosT + _y *-_sinT; // rotate new angle
	_v = _x * _sinT + _y * _cosT;

	_x = _u + _cx + _dx; // translate back + desloc
	_y = _v + _cy + _dy;

	_pt.set(_x,_y);
}

//----------------------------------------------------------------------------------------------------
//------------- rectangle
void fill4(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_QUADS);
		glTexCoord2f(_r.x,            _r.y);             glVertex2f( _x,      _y);
		glTexCoord2f(_r.x + _r.width, _r.y);             glVertex2f( _x + _w, _y);
		glTexCoord2f(_r.x + _r.width, _r.y + _r.height); glVertex2f( _x + _w, _y + _h);
		glTexCoord2f(_r.x,            _r.y + _r.height); glVertex2f( _x,      _y + _h);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border4(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x,      _y);
		glVertex2f( _x + _w, _y);
		glVertex2f( _x + _w, _y + _h);
		glVertex2f( _x,      _y + _h);
		glVertex2f( _x,      _y);
	glEnd();
}

//------------- rectangle
void fill4d(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_QUADS);
		glTexCoord2f(_r.x + _r.width*0.5, _r.y);             glVertex2f( _x + _w*0.5,_y);
		glTexCoord2f(_r.x + _r.width, _r.y + _r.height*0.5); glVertex2f( _x + _w, _y + _h*0.5);
		glTexCoord2f(_r.x + _r.width*0.5, _r.y + _r.height); glVertex2f( _x + _w*0.5,_y + _h);
		glTexCoord2f(_r.x,            _r.y + _r.height*0.5); glVertex2f( _x,      _y + _h*0.5);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border4d(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x + _w*0.5, _y);
		glVertex2f( _x + _w,     _y + _h*0.5);
		glVertex2f( _x + _w*0.5, _y + _h);
		glVertex2f( _x,          _y + _h*0.5);
		glVertex2f( _x + _w*0.5, _y);
	glEnd();
}
//------------- triangle upper
void fill3u(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_TRIANGLES);
		glTexCoord2f(_r.x + _r.width*0.5, _r.y);             glVertex2f( _x + _w*0.5,_y);
		glTexCoord2f(_r.x + _r.width,     _r.y + _r.height); glVertex2f( _x + _w,    _y + _h);
		glTexCoord2f(_r.x,                _r.y + _r.height); glVertex2f( _x,         _y + _h);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border3u(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x + _w*0.5,_y);
		glVertex2f( _x + _w,    _y + _h);
		glVertex2f( _x,         _y + _h);
		glVertex2f( _x + _w*0.5,_y);
	glEnd();
}
//------------- triangle down
void fill3d(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_TRIANGLES);
		glTexCoord2f(_r.x + _r.width*0.5, _r.y + _r.height); glVertex2f( _x + _w*0.5,_y + _h);
		glTexCoord2f(_r.x,                _r.y);             glVertex2f( _x,         _y);
		glTexCoord2f(_r.x + _r.width,     _r.y);             glVertex2f( _x + _w,    _y);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border3d(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x + _w*0.5,_y + _h);
		glVertex2f( _x,         _y);
		glVertex2f( _x + _w,    _y);
		glVertex2f( _x + _w*0.5,_y + _h);
	glEnd();
}
//-------------------------
void fill3r(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_TRIANGLES);
		glTexCoord2f(_r.x, _r.y);                            glVertex2f( _x,_y);
		glTexCoord2f(_r.x + _r.width, _r.y + _r.height*0.5); glVertex2f( _x + _w, _y + _h*0.5);
		glTexCoord2f(_r.x,            _r.y + _r.height);     glVertex2f( _x,      _y + _h);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border3r(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x,_y);
		glVertex2f( _x + _w, _y + _h*0.5);
		glVertex2f( _x,      _y + _h);
		glVertex2f( _x,_y);
	glEnd();
}
//-------------------------
void fill3l(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_TRIANGLES);
		glTexCoord2f(_r.x + _r.width, _r.y);                 glVertex2f( _x + _w, _y);
		glTexCoord2f(_r.x,            _r.y + _r.height*0.5); glVertex2f( _x,      _y + _h*0.5);
		glTexCoord2f(_r.x + _r.width, _r.y + _r.height);     glVertex2f( _x + _w, _y + _h);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border3l(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x + _w, _y);
		glVertex2f( _x,      _y + _h*0.5);
		glVertex2f( _x + _w, _y + _h);
		glVertex2f( _x + _w, _y);
	glEnd();
}
//-------------------------
void fill8(float _x, float _y, float _w, float _h, ofTexture * _tex, ofRectangle * _src) {
	ofRectangle _r(0,0,1,1); if(_src != NULL) _r = *_src;
	if(_tex != NULL) _tex->bind();
	glBegin(GL_POLYGON);
		glTexCoord2f(_r.x + _r.width*0.33, _r.y);                  glVertex2f( _x + _w*0.33, _y);
		glTexCoord2f(_r.x + _r.width*0.66, _r.y);                  glVertex2f( _x + _w*0.66, _y);
		glTexCoord2f(_r.x + _r.width,      _r.y + _r.height*0.33); glVertex2f( _x + _w,      _y + _h*0.33);
		glTexCoord2f(_r.x + _r.width,      _r.y + _r.height*0.66); glVertex2f( _x + _w,      _y + _h*0.66);
		glTexCoord2f(_r.x + _r.width*0.66, _r.y + _r.height);      glVertex2f( _x + _w*0.66, _y + _h);
		glTexCoord2f(_r.x + _r.width*0.33, _r.y + _r.height);      glVertex2f( _x + _w*0.33, _y + _h);
		glTexCoord2f(_r.x,                 _r.y + _r.height*0.66); glVertex2f( _x,           _y + _h*0.66);
		glTexCoord2f(_r.x,                 _r.y + _r.height*0.33); glVertex2f( _x,           _y + _h*0.33);
	glEnd();
	if(_tex != NULL) _tex->unbind();
}

void border8(float _x, float _y, float _w, float _h) {
	glBegin(GL_LINE_STRIP);
		glVertex2f( _x + _w*0.33, _y);
		glVertex2f( _x + _w*0.66, _y);
		glVertex2f( _x + _w,      _y + _h*0.33);
		glVertex2f( _x + _w,      _y + _h*0.66);
		glVertex2f( _x + _w*0.66, _y + _h);
		glVertex2f( _x + _w*0.33, _y + _h);
		glVertex2f( _x,           _y + _h*0.66);
		glVertex2f( _x,           _y + _h*0.33);
		glVertex2f( _x + _w*0.33, _y);
	glEnd();
}
