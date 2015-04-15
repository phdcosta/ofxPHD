#pragma once

#include "ofMain.h"

//------------------------------------------------
template<class T>
class phdInterpolator {
public:
	inline virtual bool useSubKnots() { return false; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) { T _result; return _result + mu; }
};
//------------------------------------------------
template<class T>
class phdLinearInterpolator : public phdInterpolator<T> {
public:
	inline virtual bool useSubKnots() { return false; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
		return (v2-v1) * mu + v1;
	}
};

//------------------------------------------------
template<class T>
class phdCosineInterpolator : public phdInterpolator<T> {
public:
	inline virtual bool useSubKnots() { return false; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
	   float mu2 = (1.0 - cos(mu*PI)) / 2.0;
	   return (v1 * (1.0-mu2) + v2 * mu2);
	}
};

//------------------------------------------------
template<class T>
class phdCubicInterpolator : public phdInterpolator<T> {
public:
	inline virtual bool useSubKnots() { return true; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
		float mu2 = mu * mu;
		float mu3 = mu * mu2;
		return((v3 - v2 - v0 + v1) * mu3 + (v0 - v1 - (v3 - v2 - v0 + v1)) * mu2 + (v2 - v0) * mu + v1);
	}
};

//------------------------------------------------
template<class T>
class phdCatmullromInterpolator : public phdInterpolator<T> {
public:
	inline virtual bool useSubKnots() { return true; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
		float mu2 = mu * mu;
		float mu3 = mu2 * mu;
		return ((2.0 * v1) + 
			   (-v0 + v2) * mu + 
			   (2.0*v0 - 5.0*v1 + 4.0*v2 - v3) * mu2 + 
			   (-v0 + 3.0*v1 - 3.0*v2 + v3) * mu3) * 0.5f;
	}
};

//------------------------------------------------
template<class T>
class phdBezierInterpolator : public phdInterpolator<T> {
public:
	inline virtual bool useSubKnots() { return true; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
		float mu2 = mu * mu;
		float mu3 = mu * mu2;
		return ((1.0/6.0)*
				(mu3*(-v0+3.0*v1-3.0*v2+v3)+
				 mu2*(3.0*v0-6.0*v1+3.0*v2)+
				 mu* (-3.0*v0+3.0*v2)+(v0+4.0*v1+v2)));


		/*return 1.0/6.0 *
			( ( v0 + 4.0*v1 + v2 )
			  + mu * ( ( -3.0*v0 + 3.0*v2 )
			  + mu * ( ( 3.0*v0 - 6.0*v1 + 3.0*v2 )
			  + mu * ( -v0 + 3.0*v1 - 3.0*v2 + v3 ) ) )
			);*/
	}
};

//------------------------------------------------
template<class T>
class phdHermiteInterpolator : public phdInterpolator<T> {
protected:
	float bias; //Bias: (0) is even, (+) is towards first segment, (-) towards the other
	float tension; // Tension: 1 is high, 0 normal, -1 is low
public:
	phdHermiteInterpolator() { bias = 1.0; tension = 0.0; }
	inline virtual bool useSubKnots() { return true; } // this interpolator uses subKnots to get point in perimeter
	inline virtual T calculate(float mu, T & v0, T & v1, T & v2, T & v3) {
		T m0,m1;
		float mu2,mu3;
		float a0,a1,a2,a3;
		mu2 = mu * mu;
		mu3 = mu2 * mu;
		m0 =  (v1-v0)*(1.0+bias)*(1.0-tension)/2.0;
		m0 =  m0 + (v2-v1)*(1.0-bias)*(1.0-tension)/2.0;
		m1 =  (v2-v1)*(1.0+bias)*(1.0-tension)/2.0;
		m1 =  m1 + (v3-v2)*(1.0-bias)*(1.0-tension)/2.0;
		a0 =  2.0*mu3 - 3.0*mu2 + 1.0;
		a1 =      mu3 - 2.0*mu2 + mu;
		a2 =      mu3 -     mu2;
		a3 = -2.0*mu3 + 3.0*mu2;
		return (a0*v1 + a1*m0 + a2*m1 + a3*v2);
	}
};

//------------------------------------------------
enum phdInterpolatorType { pcitLinear, pcitCubic, pcitCatmullrom, pcitBezier, pcitHermite, pcitCosine };

//------------------------------------------------
template<class T>
class phdInterpolatedCurve {
protected:
	class segment {
	protected:
		float len_;
		T * A_;
		T * B_;
		vector<T> subKnots_;
		vector<float> distances_;
	public:
		segment() { A_ = NULL; B_ = NULL; len_ = 0.0; }
		virtual ~segment() { A_ = NULL; B_ = NULL; len_ = 0.0; }
		inline T & A() { return *A_; }
		inline T & B() { return *B_; }
		inline float len() { return len_; }
		inline void setLen(float _len) { len_ = _len; }
		inline void set(T * _a, T * _b, float _len) { A_ = _a; B_ = _b; len_ = _len; }
		inline vector<T> & subKnots() { return subKnots_; }
		inline vector<float> & distances() { return distances_; }

		inline const T getAt(float _position) { return getAtDistance(len_ * MAX( 0.0, MIN(1.0, _position)) ); }

		inline const T getAtDistance(float _distance) {
			float _d = MAX(0.0, MIN(len_, _distance));
			int   _i = 0;
			float _len = 0.0;
			while(_len + distances_[_i] < _d && _i < distances_.size()-1) {
				_len += distances_[_i];
				_i   += 1;
			}
			float _pos = ((_d-_len)/distances_[_i]); // pos inside subSegment

			// linear interpolation inside the subSegment
			return T(subKnots_[_i] + ((subKnots_[_i+1] - subKnots_[_i]) * _pos));
		}
	};
	vector<segment> segs_;

	float precision_;  // precision indicates distance between values of base curve
	float tolerance_;  // tolerance of base curve

	float curveLength_;
	vector<T*> knots_;

	T phantomA, phantomB; // phantom points

	inline void freeItems() {
		for(int i = knots_.size()-1; i > -1; i--) del(i);
		knots_.clear();
		onUpdateCurve();
	}

	inline void computeBaseCurveAndLength() {
		
		curveLength_ = 0.0;
		if(knots_.size() < 2) { segs_.clear(); return; }

		float _segLength = 0.0;

		segs_.front().setLen(0.0); // first seg is a phantom
		segs_.front().distances().resize(0);
		segs_.front().subKnots().resize(0);

		if(interpolator_->useSubKnots()) { // this interpolator need subKnots?
			for(int i = 1; i < segs_.size()-1; i++) { // ignore first and last segs

				T _a, _b, _e;

				vector<T> _ks;

				// working with square of distance to avoid squareRoot calculation of distances
				float _maxDistance = pow(precision_ + tolerance_, 2);
				float _minDistance = pow(precision_ - tolerance_, 2);

				// get begining of segment
				_a = interpolator_->calculate(0.0, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());

				// get end of segment
				_e = interpolator_->calculate(1.0, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());

				_ks.push_back(_a);

				float _d2End = _ks.back().distanceSquared(_e); // distance to the End of segment

				float _muF = 0.0; // mu where the last point was founded
				float _muC = 0.5; // current mu for interpolator

				while(_d2End > _maxDistance) { // can we interpolate more?

					bool _foundedInsidePrecision = false;
					bool _foundedByMinDifference = false;

					while(!_foundedInsidePrecision && !_foundedByMinDifference) { // lets try n times

						_b = interpolator_->calculate(_muC, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());

						float _dab = _a.distanceSquared(_b); // distance between _a and _b

						if(_dab < _minDistance) {
							_muC = _muC + ((_muC - _muF) * 0.5);
						} else if(_dab > _maxDistance) {
							_muC = _muC - ((_muC - _muF) * 0.5);
						} else if((_muC-_muF) < 0.01) { // minimum mu = 1/100th of curve
							_foundedByMinDifference = true;
						} else {
							_foundedInsidePrecision = true;
						}
					}
					_muF = _muC; // founded at that mu -- saves it

					_ks.push_back(_b); // add founded value
					_a = _b; // next time start from founded value

					_d2End = _ks.back().distanceSquared(_e); // distance to the End of segment

					_muC = _muC + (1.0 - _muF) * 0.5;
				}

				// make room for subknots and distances of this segment
				segs_[i].subKnots().resize(_ks.size() + 1);
				segs_[i].distances().resize(_ks.size());

				// add founded values to segment subKnots
				for(int k = 0; k < _ks.size(); k++) { segs_[i].subKnots()[k] = _ks[k]; }

				// calculates and add last subKnot == using mu == 1.0
				_b = interpolator_->calculate(1.0, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());
				segs_[i].subKnots().back() = _b; // last subKnot is at mu == 1.0

				// calculate distances between subKnots
				_segLength = 0.0; // length of this segment
				for(int k = 0; k < segs_[i].subKnots().size()-1; k++) {
					segs_[i].distances()[k] = segs_[i].subKnots()[k].distance(segs_[i].subKnots()[k+1]);
					_segLength += segs_[i].distances()[k];
				}

				segs_[i].setLen(_segLength);
				curveLength_ += _segLength;
			}
		} else { // if this kind of interpolator dont need subKnots, the segment length is distance between A and B
			for(int i = 1; i < segs_.size()-1; i++) { // ignore first and last segs
				
				segs_[i].subKnots().resize(2);  // make room for subknots and distances of this segment
				segs_[i].distances().resize(1);

				segs_[i].subKnots().front() = interpolator_->calculate(0.0, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());
					//segs_[i].A(); // first subKnot == A()
				segs_[i].subKnots().back()  = interpolator_->calculate(1.0, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B());
					//segs_[i].B(); // last  subKnot == B()

				segs_[i].distances().front() = segs_[i].subKnots().front().distance(segs_[i].subKnots().back()); // gets distance front->back

				segs_[i].setLen(segs_[i].distances().front()); // distance between front->back

				curveLength_ += segs_[i].len();
			}
		}

		segs_.back().setLen(0.0); // last seg is a phantom
		segs_.front().distances().resize(0);
		segs_.front().subKnots().resize(0);
	}

	inline void onUpdateCurve() { // on update curve recalculate segments

		curveLength_ = 0.0;
		if(knots_.size() < 2) { segs_.clear(); return; }

		phantomA = (2.0 * (*knots_.front()) - (*knots_[1]));
		phantomB = (2.0 * (*knots_.back()) - (*knots_[knots_.size()-2]));

		segs_.resize(knots_.size() + 1); // make room for phantom segments
		segs_.front().set(&phantomA, knots_.front(), 0.0);
		for(int i = 1; i < knots_.size(); i++) { segs_[i].set(knots_[i-1], knots_[i], 0.0); }
		segs_.back().set(knots_.back(), &phantomB, 0.0);

		computeBaseCurveAndLength();
	}

	phdInterpolator<T> * interpolator_;
	phdInterpolatorType type_;

public:
	phdInterpolatedCurve() {
		precision_ = 40.0; // distance between values on the base curve
		tolerance_ = 1.0; // values will be inside (_precision +/- _tolerance)
		type_ = pcitLinear; // default is an linear curve
		interpolator_ = new phdLinearInterpolator<T>; // creates a linear interpolator
	}
	virtual ~phdInterpolatedCurve() { freeItems(); delete interpolator_; }

	inline int size() { return knots_.size(); }

	inline T & get(unsigned int _index) { return *knots_[_index]; }
	inline void set(unsigned int _index, const T _value) { knots_[_index]->set(T); onUpdateCurve(); }

	inline phdInterpolator<T> & interpolator() { return *interpolator_; }

	inline vector<segment> & segs() { return segs_; }

	inline float curveLength() { return curveLength_; }

	inline void setup(phdInterpolatorType _type, float _precision, float _tolerance) {
		precision_ = _precision;
		tolerance_ = _tolerance;
		setType(_type);
	}

	inline float precision() { return precision_; }
	inline void setPrecision(float _value) { precision_ = _value; computeBaseCurveAndLength(); }

	inline float tolerance() { return tolerance_; }
	inline void setTolerance(float _value) { tolerance_ = _value; computeBaseCurveAndLength(); }

	inline phdInterpolatorType type() { return type_; }

	inline void setType(phdInterpolatorType _value) {
		if(_value != type_) {
			delete interpolator_;
			type_ = _value;
			switch(type_) {
				case pcitCubic      : interpolator_ = new phdCubicInterpolator<T>; break;
				case pcitCatmullrom : interpolator_ = new phdCatmullromInterpolator<T>; break;
				case pcitBezier     : interpolator_ = new phdBezierInterpolator<T>; break;
				case pcitHermite    : interpolator_ = new phdHermiteInterpolator<T>; break;
				case pcitCosine     : interpolator_ = new phdCosineInterpolator<T>; break;
				default             : interpolator_ = new phdLinearInterpolator<T>; type_ = pcitLinear; break;
			}
			computeBaseCurveAndLength();
		}
	}

	void add(const T _value) {
		knots_.push_back(new T(_value));
		onUpdateCurve();
	}

	void del(unsigned int _index) {
		if(_index < knots_.size()) {
			delete knots_[_index]; knots_[_index] = NULL;
			knots_.erase(knots_.begin() + _index);
			onUpdateCurve();
		}
	}

	void getBaseCurve(vector<T> & _output) {
		for(int i = 1; i < segs_.size()-1; i++) {
			for(int k = 0; k < segs_[i].subKnots().size()-1; k++) { // dont include last subknot
				_output.push_back(segs_[i].subKnots()[k]); // subknots
			}
		}
		_output.push_back(segs_[segs_.size()-1].A()); // last subknot == last knot
	}

	void computeCurve(vector<T> & _output, int _dotsPerSegment = 5) {
		
		if(knots_.size() == 0) { _output.clear(); return; } // no points

		if(knots_.size() == 1) { _output.clear(); _output.push_back(*(knots_[0])); return; } // just one point

		// at least 2 points means at least 4 segments with two phantom that A == B
		float _step = 1.0 / (float)_dotsPerSegment;

		for(int i = 1; i < segs_.size()-1; i++) {
			for(float k = 0.0; k < 1.0; k += _step) {
				_output.push_back(interpolator_->calculate(k, segs_[i-1].A(), segs_[i].A(), segs_[i].B(), segs_[i+1].B()));
			}
		}
		int _s = segs_.size()-1;
		_output.push_back(interpolator_->calculate(1.0, segs_[_s-2].A(), segs_[_s-1].A(), segs_[_s-1].B(), segs_[_s].B()));

		//_output.push_back(T(*knots_.back()));
	}

	void computeEquidistanteCurve(vector<T> & _output, int _dotsPerCurve = 25) {
		if(knots_.size() == 0) { _output.clear(); return; } // no points
		if(knots_.size() == 1) { _output.clear(); _output.push_back(T(*(knots_[0]))); return; } // just one point
		float _step = curveLength_ / ((float)_dotsPerCurve-1.0);
		_output.push_back(T(*knots_.front()));
		for(float i = _step; i < curveLength_; i += _step) _output.push_back(getAtDistance(i));
		_output.push_back(T(*knots_.back()));
	}

	const T getAt(float _position) { return getAtDistance(MAX(0.0, MIN(1.0, _position)) * curveLength_); }

	const T getAtDistance(float _distance) {

		if(knots_.size() == 0) return T(0);
		if(knots_.size() == 1) return T(*knots_[0]);

		if(_distance < 0.0) return T(*knots_.front());
		if(_distance >= curveLength_) return T(*knots_.back());

		int   _i = 1;
		float _len = 0.0;
		while(_len + segs_[_i].len() < _distance && _i < segs_.size()-2) {
			_len += segs_[_i].len();
			_i   += 1;
		}

		return T(segs_[_i].getAtDistance(_distance-_len)); // get using desired distance inside segment
	}
};

typedef phdInterpolatedCurve<ofVec2f> phdPath2d;
typedef phdInterpolatedCurve<ofVec3f> phdPath3d;