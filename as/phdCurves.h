#pragma once

#include "ofMain.h"

#define TWEEN_EPSILON 0.00000001

// invert a value inside a range if _I condition is true
#define INVIF(_I, _R, _V) ( ((_I) ? _R-_V : _V) )

//---------------------------------------------------------------------------------------------------------------
class phdCurveCalculator {
protected:
	float pa, pb;

public:
	phdCurveCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	inline void set(float _a, float _b) { pa = _a; pb = _b; }

	virtual float calculate(float x) { return x; } // linear

	virtual float calculateMap(float x, float vA, float vB) { return ofMap(calculate(x), 0.0, 1.0, vA, vB); }

	virtual float calculateMed(float x, float vA, float vB, float vM) { return (x <= 0.5) ? calculateMap(x*2.0,vA,vM) : calculateMap((1.0-x)*2.0,vB,vM); }

	virtual void draw(float x, float y, float w, float h, float res = 5.0, bool inv = false) {

		ofNoFill();
		ofRect(x, y, w, h);

		glBegin(GL_LINE_STRIP);
			double _sp = 1.0 / (w / res); // line segment each RES px
			double _cp = 0.0;
			while(_cp < 1.0) {
				glVertex2d(x + INVIF(inv, w, (w * _cp)), y + h - (calculate(_cp) * h));
				_cp += _sp;
			}
			glVertex2d(x + INVIF(inv, w, (w * 1.0)), y + h - (calculate(1.0) * h));
		glEnd();
		ofCircle(x + INVIF(inv, w, 0), y + h - (calculate(0.0) * h), 2);
		ofCircle(x + INVIF(inv, w, w), y + h - (calculate(1.0) * h), 2);
	}

	void mapCurvePosToArea(double _pos, float x, float y, float w, float h, float & u, float & v, float inv = false) {
		u = x + ofMap(_pos, 0.0, 1.0, (inv?w:0.0), (inv?0.0:w), true);
		v = y + ofMap(calculate(_pos), 0.0, 1.0, h, 0.0, true);
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdExponentialEasingCalculator : public phdCurveCalculator {
public:
	phdExponentialEasingCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
  
		// emphasis
		if (a < 0.5) return powf(x, 2.0*(a));

		// de-emphasis
		return powf(x, 1.0/(1.0-(2.0*(a-0.5))));
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdExponentialSeatCalculator : public phdCurveCalculator  {
public:
	phdExponentialSeatCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));

		if (x <= 0.5) return (pow(2.0*x, 1.0-a))/2.0;

		return 1.0 - (pow(2.0*(1.0-x), 1.0-a))/2.0;
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdExponentialSigmoidCalculator : public phdCurveCalculator {
public:
	phdExponentialSigmoidCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));

		a = 1.0-a; // for sensible results
  
		if (x <= 0.5) return (pow(2.0*x, 1.0/a))/2.0;

		return 1.0 - (pow(2.0*(1.0-x), 1.0/a))/2.0;
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdExponentialLogisticCalculator : public phdCurveCalculator {
public:
	phdExponentialLogisticCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) { // n.b.: this Logistic Sigmoid has been normalized.

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		a = (1.0/(1.0-a) - 1.0);

		float A = 1.0 / (1.0 + exp(0 -((x-0.5)*a*2.0)));
		float B = 1.0 / (1.0 + exp(a));
		float C = 1.0 / (1.0 + exp(0-a)); 
		return (A-B)/(C-B);
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdCircularEaseInCalculator : public phdCurveCalculator {
public:
	phdCircularEaseInCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {
		return 1.0 - sqrt(1.0 - x*x);
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdCircularEaseOutCalculator : public phdCurveCalculator {
public:
	phdCircularEaseOutCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {
		return sqrt(1.0 - (1.0-x)*(1.0-x));
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdDoubleCircleSeatCalculator : public phdCurveCalculator {
public:
	phdDoubleCircleSeatCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }
	virtual float calculate(float x) {
		float a = MAX(0.0, MIN(1.0, pa)); 
		if (x <= a) {
			float y = (x == 0.0 ? 0.0 : (a*a) - (x-a)*(x-a));
			return (y == 0.0 ? 0.0 : sqrt(y));
		} else {
			float y = (x == 1.0 ? 0.0 : (1.0-a)*(1.0-a) - (x-a)*(x-a));
			return 1.0 - (y == 0.0 ? 0.0 : sqrt(y));
		}
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdDoubleCircleSigmoidCalculator : public phdCurveCalculator {
public:
	phdDoubleCircleSigmoidCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }
	virtual float calculate(float x) {
		float a = MAX(0.0, MIN(1.0, pa)); 
		if (x <= a){
			return (a - sqrt(a*a - x*x));
		} else {
			return (a + sqrt((1.0-a)*(1.0-a) - (x-1.0)*(x-1.0)));
		}
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdDoubleEllipticSeatCalculator : public phdCurveCalculator {
public:
	phdDoubleEllipticSeatCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		float b = MAX(0.0, MIN(1.0, pb)); 

		if (x <= a) return (b/a) * sqrt((a*a) - (x-a)*(x-a));

		return 1.0 - ((1.0-b)/(1.0-a))*sqrt((1-a)*(1-a) - (x-a)*(x-a));
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdDoubleEllipticSigmoidCalculator : public phdCurveCalculator {
public:
	phdDoubleEllipticSigmoidCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		float b = MAX(0.0, MIN(1.0, pb)); 
 
		if (x <=a ) return b * (1.0 - (sqrt((a*a) - (x*x))/a));

		return b + ((1.0-b)/(1.0-a))*sqrt((1.0-a)*(1.0-a) - (x-1.0)*(x-1.0));
	}
};


//---------------------------------------------------------------------------------------------------------------
class phdDoubleCubicSeatCalculator : public phdCurveCalculator {
public:
	phdDoubleCubicSeatCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		float b = MAX(0.0, MIN(1.0, pb)); 
  
		if (x <= a) return b - b*pow(1.0-x/a, 3.0);

		return b + (1.0-b)*pow((x-a)/(1.0-a), 3.0);
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdDoubleCubicSeatLinearBlendCalculator : public phdCurveCalculator {
public:
	phdDoubleCubicSeatLinearBlendCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		float b = MAX(0.0, MIN(1.0, pb)); 
		b = 1.0 - b; //reverse for intelligibility.
  
		if (x <= a) return b*x + (1.0-b)*a*(1.0-pow(1.0-x/a, 3.0));

		return b*x + (1.0-b)*(a + (1.0-a)*pow((x-a)/(1.0-a), 3.0));
	}
};

//---------------------------------------------------------------------------------------------------------------
class phdQuadraticThroughGivenPointCalculator : public phdCurveCalculator {
public:
	phdQuadraticThroughGivenPointCalculator(float _a = 0.5, float _b = 0.5) { pa = _a; pb = _b; }

	virtual float calculate(float x) {

		float a = MAX(TWEEN_EPSILON, MIN(1.0-TWEEN_EPSILON, pa));
		float b = MAX(0.0, MIN(1.0, pb)); 
  
		float A = (1.0-b)/(1.0-a) - (b/a);
		float B = (A*(a*a)-b)/a;
		float y = A*(x*x) - B*(x);

		return MAX(0.0, MIN(1.0, y)); 
	}
};

typedef std::map<string, phdCurveCalculator*> phdMapCurveDef;
typedef phdMapCurveDef::iterator phdMapCurveItemDef;

//---------------------------------------------------------------------------------------------------------------
class phdCurveSystem {
protected:
	phdMapCurveDef items;
	vector<phdCurveCalculator*> curves;

public:
	phdCurveSystem() {
		addCurveCalculator("LINEAR",			new phdCurveCalculator());
		addCurveCalculator("EXP_EASING_IN",		new phdExponentialEasingCalculator(0.867));
		addCurveCalculator("EXP_EASING_OUT",	new phdExponentialEasingCalculator(0.220));
		addCurveCalculator("EXP_SEAT_A",		new phdExponentialSeatCalculator(0.247)); 
		addCurveCalculator("EXP_SEAT_B",		new phdExponentialSeatCalculator(0.607));
		addCurveCalculator("EXP_SEAT_C",		new phdExponentialSeatCalculator(0.807)); 
		addCurveCalculator("EXP_SIGMOID_A",		new phdExponentialSigmoidCalculator(0.367)); 
		addCurveCalculator("EXP_SIGMOID_B",		new phdExponentialSigmoidCalculator(0.727));
		addCurveCalculator("EXP_SIGMOID_C",		new phdExponentialSigmoidCalculator(0.887)); 
		addCurveCalculator("EXP_LOGISTIC_A",	new phdExponentialLogisticCalculator(0.657)); 
		addCurveCalculator("EXP_LOGISIIC_B",	new phdExponentialLogisticCalculator(0.787)); 
		addCurveCalculator("EXP_LOGISIIC_C",	new phdExponentialLogisticCalculator(0.920)); 
		addCurveCalculator("CIRCLE_EASING_IN",	new phdCircularEaseInCalculator()); 
		addCurveCalculator("CIRCLE_EASING_OUT",	new phdCircularEaseOutCalculator()); 
		addCurveCalculator("CIRCLE_SEAT_A",		new phdDoubleCircleSeatCalculator(0.293));
		addCurveCalculator("CIRCLE_SEAT_B",		new phdDoubleCircleSeatCalculator(0.500));
		addCurveCalculator("CIRCLE_SEAT_C",		new phdDoubleCircleSeatCalculator(0.770)); 
		addCurveCalculator("CIRCLE_SIGMOID_A",	new phdDoubleCircleSigmoidCalculator(0.293)); 
		addCurveCalculator("CIRCLE_SIGMOID_B",	new phdDoubleCircleSigmoidCalculator(0.500)); 
		addCurveCalculator("CIRCLE_SIGMOID_C",	new phdDoubleCircleSigmoidCalculator(0.747)); 
		addCurveCalculator("ELLIPTIC_SEAT_A",	new phdDoubleEllipticSeatCalculator(0.373, 0.747)); 
		addCurveCalculator("ELLIPTIC_SEAT_B",	new phdDoubleEllipticSeatCalculator(0.693, 0.187)); 
		addCurveCalculator("ELLIPTIC_SIGMOID_A",new phdDoubleEllipticSigmoidCalculator(0.693, 0.507)); 
		addCurveCalculator("ELLIPTIC_SIGMOID_B",new phdDoubleEllipticSigmoidCalculator(0.253, 0.513)); 
		addCurveCalculator("CUBIC_SEAT_A",		new phdDoubleCubicSeatCalculator(0.407, 0.720)); 
		addCurveCalculator("CUBIC_SEAT_B",		new phdDoubleCubicSeatCalculator(0.607, 0.247));
		addCurveCalculator("CUBIC_SEAT_LB_A",	new phdDoubleCubicSeatLinearBlendCalculator(0.640, 0.827)); 
		addCurveCalculator("CUBIC_SEAT_LB_B",	new phdDoubleCubicSeatLinearBlendCalculator(0.347, 0.887)); 
		addCurveCalculator("QUADRATIC_POINT_A",	new phdQuadraticThroughGivenPointCalculator(0.350, 0.650));
		addCurveCalculator("QUADRATIC_POINT_B",	new phdQuadraticThroughGivenPointCalculator(0.650, 0.350)); 
	}

	~phdCurveSystem() {
		for(phdMapCurveItemDef it = items.begin(); it != items.end(); ++it) { it->second = NULL; }
		for(int i = 0; i < curves.size(); i++) { delete curves[i]; curves[i] = NULL; }
		curves.clear();
	}

	inline int size() { return curves.size(); }

	phdCurveCalculator & operator[](unsigned int _index) { return *getCurveByIndex(_index); }
	phdCurveCalculator & get(unsigned int _index)        { return *getCurveByIndex(_index); }

	phdCurveCalculator & operator[](string _name) { return *getCurveByName(_name); }

	inline phdCurveCalculator * getCurveByIndex(unsigned int _index) { return curves[MAX(0, MIN(curves.size()-1, _index))]; }

	inline phdCurveCalculator * getCurveByName(string _name) { return getItemByName(_name)->second; }

	inline phdMapCurveItemDef getItemByName(string _name) {
		phdMapCurveItemDef it = items.find(_name);
		return (it == items.end() ? ++items.begin() : it);
	}

	inline phdMapCurveItemDef getItemByIndex(unsigned int _index) {
		phdMapCurveItemDef it = items.begin();
		std::advance(it, MAX(0, MIN(items.size(), _index)));
		return it;
	}

	inline phdMapCurveItemDef getItemByCurveIndex(unsigned int _index) {
		for(phdMapCurveItemDef it = items.begin(); it != items.end(); ++it) {
			if(it->second == curves[_index]) return it;
		}
		return --items.end();
	}

	inline int indexByName(string _name) {
		phdMapCurveItemDef it = items.find(_name);
		if(it == items.end()) return -1; // name doesnt exist
		for(int i = 0; i < curves.size(); i++) {
			if(it->second == curves[i]) return i;
		}
		return -1;
	}

	inline const string nameByIndex(unsigned int _index) {
		if(_index >= curves.size()) return "";
		for(phdMapCurveItemDef it = items.begin(); it != items.end(); ++it) {
			if(it->second == curves[_index]) return it->first;
		}
		return "";
	}

	void addCurveCalculator(string curveName, phdCurveCalculator * curveCalculator) {
		curves.push_back(curveCalculator);
		items.insert(make_pair(curveName, curves[curves.size()-1]));
	}

	// calculates the bounds of a specific curve inside a big area used for all curves
	inline bool getCurveBoundsInsideDrawingArea(unsigned int _index, float x, float y, float w, float h, ofRectangle & _bounds) {

		if(_index >= curves.size()) return false;

		float _gap = 15;
		float _x = x, _y = y, _w = (w-_gap*5.0) / 6.0, _h = (h-_gap*4.0) / 5.0;

		for(int _curve = 0; _curve < _index; _curve++) {
			_x += _w + _gap;
			if(_x > w) { _x = x; _y += _h + _gap; }
		}

		_bounds.set(_x,_y,_w,_h);

		return true;
	}

	void drawCurves(float x, float y, float w, float h) {

		float _gap = 15;
		float _x = x, _y = y, _w = (w-_gap*5.0) / 6.0, _h = (h-_gap*4.0) / 5.0;

		for(int _curve = 0; _curve < curves.size(); _curve++) {

			curves[_curve]->draw(_x, _y, _w, _h);

			_x += _w + _gap;
			if(_x > w) { _x = x; _y += _h + _gap; }
		}
	}
};
//---------------------------------------------------------------------------------------------------------------
// global curve system
//---------------------------------------------------------------------------------------------------------------
static phdCurveSystem phdMainCurves;

phdCurveSystem & phdMainCurveSystem();

//---------------------------------------------------------------------------------------------------------------
// when try to add points after one
// pcaReposAll --> repos all points to new scaled size
//		if you add after 1.0 the curve grows to new pos and is scaled back to 1.0
//		if you delete the last point, the curve shrinks and is scaled back to n-1 point pos
// pcaIgnoreOp
//		ignores add after 1.0 pos
//		dont changes the last point pos

enum phdCurveEditLastPointMode { pcaAddReposAll, pcaIgnoreAdd };

//---------------------------------------------------------------------------------------------------------------
class phdMultPointCurveCalculator : public phdCurveCalculator {
protected:
	struct phdPosPoint {
		double value;
		double pos;
		phdCurveCalculator * curve;
		phdPosPoint(double _pos = 0.0, double _value = 0.0, phdCurveCalculator * _curve = NULL) : pos(_pos), value(_value), curve(_curve) { }
	};

	struct phdSegment {
		double va, vb;
		float  ta, tb;
		phdCurveCalculator * curve;
		phdSegment(double _va = 0.0, double _vb = 0.0, float _ta = 0.0, float _tb = 0.0, phdCurveCalculator * _curve = NULL) : va(_va), vb(_vb), ta(_ta), tb(_tb), curve(_curve) { }
	};

	vector<phdPosPoint> points;
	phdCurveEditLastPointMode pointEditMode;

public:
	phdMultPointCurveCalculator(float _startValue = 0, float _endValue = 1.0, int _curveIndex = 0);

	virtual float calculate(float x) { return getValueInPos(x); }

	inline const phdPosPoint & getPoint(unsigned int _index) const { return points[_index]; }
	inline int getPointCount() { return points.size(); }

	void setPointPos(unsigned int _index, double _pos);
	void setPointValue(unsigned int _index, double _value);

	inline int getSegmentCount() { return points.size()-1; }

	void setSegmentCurveMode(unsigned int _index, int _curveIndex);
	void setSegmentCurveMode(unsigned int _index, string _curveName);

	int segmentIndexByPos(double _pos);
	bool getSegment(unsigned int _index, phdSegment & _segment);

	double getValueInPos(double _pos);
	double getValueInSegmentByPos(phdSegment & _seg, double _pos);

	void addPoint(double _pos, double _value, int _curveIndex);
	void delPoint(unsigned int _index);

	void normalizePoints();

	virtual void draw(float x, float y, float w, float h, float res = 5.0, bool inv = false);
};
