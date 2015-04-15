#pragma once

#include "ofMain.h"

#include "phdCurves.h"

//----------------------------------------
enum phdTickerLoopMode   { palmNormal = 0, palmPingPong };
enum phdTickerStartMode  { pasmStartForward = 0, pasmStartBackward };
enum phdTickerStatusType { psmdArmed, psmdStopped, psmdPaused, psmdPlaying };
enum phdTickerEventType  { paetOnPlay, paetOnStop, paetOnPause, paetOnCycleCompleted, paetOnActivate, paetOnDeactivate };

//----------------------------------------
#define FRACINT(_V)  ( (_V < 0) ? (_V +  (int) (-_V)) : (_V -  (int) _V) )
#define FRACLONG(_V) ( (_V < 0) ? (_V + (long) (-_V)) : (_V - (long) _V) )

enum phdLoopBoundsMode      { plbmLoopOver, plbmLoopBelow, plbmLoopForever };
enum phdLoopModeType        { plmtNone, plmtNormal, plmtPingPong };
enum phdAnimationStatusType { pastStopped, pastPaused, pastPlaying };

//----------------------------------------
class phdTheClock;

//----------------------------------------
class phdTheClockEventData {
protected:
	phdTheClock * clock_;
	phdTickerEventType type_;
public:
	phdTheClockEventData(phdTheClock * _clock, phdTickerEventType _type) { clock_ = _clock; type_ = _type; }
	inline phdTickerEventType type() { return type_; }
	inline phdTheClock * clock() { return clock_; }
};

enum phdTheClockMode { ptcmOFTimef, ptcmOFFrame, ptcmUserData };

//----------------------------------------
class phdTheClock {
protected:
	double internalClock_;
	double velocity_;
	double lastUpdate_;
	phdTickerStatusType status_;

	ofEvent<phdTheClockEventData> events_;
	inline void notify(phdTickerEventType _type) {
		phdTheClockEventData _data(this, _type);
		ofNotifyEvent(events_, _data);
	}

	bool usingEvents;
	phdTheClockMode mode_;

	double userData_;

	inline double getAutoUpdateValue() {
		switch(mode_) {
			case ptcmOFTimef  : return ofGetElapsedTimef(); break;
			case ptcmOFFrame  : return ofGetFrameNum(); break;
			default           : return userData_; break;
		}
	}

	inline void updateClock(bool _force = false) {
		bool _wasArmed = (status_ == psmdArmed);
		if(_force || _wasArmed || isPlaying()) {
			double _now = getAutoUpdateValue();
			internalClock_ = internalClock_ + ( (_now - lastUpdate_) * velocity_ );
			lastUpdate_ = _now;
			if(_wasArmed) { status_ = psmdPlaying; notify(paetOnPlay); }
			//if(internalClock_ <= 0.0) stop();
		}
	}

public:
	phdTheClock() {
		usingEvents = false;
		velocity_ = 1.0;
		internalClock_ = lastUpdate_ = userData_ = 0.0;
		status_ = psmdStopped;
		mode_ = ptcmUserData; setMode(ptcmOFTimef); 
	}
	virtual ~phdTheClock() { disableEvents(); }

	ofEvent<phdTheClockEventData> events() { return events_; }

	inline void setMode(phdTheClockMode _value) {
		if(_value != mode_) {
			if(_value == ptcmUserData) { disableEvents(); } else { enableEvents(); }
			mode_ = _value;
		}
	}

	inline double get() { return internalClock_; }

	inline void set(double _value) {
		if(mode_ == ptcmUserData) {
			userData_ = internalClock_ = lastUpdate_ = _value;
		} else {
			internalClock_ = _value;
			userData_ = lastUpdate_ = getAutoUpdateValue();
		}
	}

	inline void onOFupdate(ofEventArgs & _args) { updateClock(false); }

	inline bool isPaused()   { return status_ == psmdPaused; }
	inline bool isPlaying()  { return status_ == psmdPlaying; }
	inline bool isAutoMode() { return mode_   != ptcmUserData; }

	inline double velocity() { return velocity_; }
	inline void setVelocity(double _value) { velocity_ = _value; }

	inline virtual void play() {
		if(!isPlaying()) {
			lastUpdate_ = getAutoUpdateValue();
			status_ = psmdArmed;
		}
	}
	inline virtual void stop() {
		bool _wasPlaying = isPlaying();
		internalClock_ = 0.0; status_ = psmdStopped;
		if(_wasPlaying) notify(paetOnStop);
	}
	inline virtual void pause() {
		bool _wasPlaying = isPlaying();
		status_ = psmdPaused;
		if(_wasPlaying) notify(paetOnPause);
	}

	void enableEvents() { if(!usingEvents) { ofAddListener(ofEvents().update, this, &phdTheClock::onOFupdate); usingEvents = true; } }
	void disableEvents() { if(usingEvents) { ofRemoveListener(ofEvents().update, this, &phdTheClock::onOFupdate); usingEvents = false; } }

};

//----------------------------------------
static phdTheClock phdMasterClock;

phdTheClock & phdMainClock();

//----------------------------------------
class phdTicker;

//----------------------------------------
class phdTickerEventData {
protected:
	phdTicker * ticker_;
	phdTickerEventType type_;
public:
	phdTickerEventData(phdTicker * _ticker, phdTickerEventType _type) {
		ticker_ = _ticker; type_ = _type;
	}
	inline phdTickerEventType type() { return type_; }
	inline phdTicker * ticker() { return ticker_; }
};

//template<typedef T>
class phdTicker {
protected:
	double duration_;    // scalar value in time units
	double executions_;  // scalar number of executions greater than 0
	double velocity_;    // scalar velocity

	double lastUpdate_;  // when was last update
	double stepSize_;    // how many units each deltaClock

	double curClock_;    // current scalar that control this anim inside its own track
	double curExec_;     // executions == 0.0 ? INFINITE else RANGE(0, executions_)
	double curPos_;      // normalized inside duration_
	double curVel_;      // real velocity of this
	double dirExec_;      // just to control direction on ping pong loops

	ofVec2f curValue_;   // which object is varying by that animation

	unsigned int curveIndex_; // which curve to use from phdMainCurveSystem;  0 == Linear

	phdTheClock * theClock_;
	phdTheClock & clock() { return (theClock_ == NULL ? phdMainClock() : *theClock_); }

	phdTicker * blendTo_;

	bool active_;
	phdLoopModeType loopMode_;
	phdLoopBoundsMode boundsMode_;
	phdAnimationStatusType status_;

	ofEvent<phdTickerEventData> events_;
	inline void notify(phdTickerEventType _type) {
		phdTickerEventData _data(this, _type);
		ofNotifyEvent(events_, _data);
	}

public:
	phdTicker(phdTheClock * _clock = NULL) {
		setDuration(0.0);   // to force calculation of stepSize_
		setExecutions(1.0); // to keep executions != 0.0
		velocity_ = dirExec_ = 1.0;
		curExec_ = curPos_ = curVel_ = curClock_ = 0.0;
		blendTo_ = NULL;
		theClock_ = _clock;
		loopMode_ = plmtNormal;
		boundsMode_ = plbmLoopForever;
		status_ = pastStopped;
		active_ = false;
		lastUpdate_ = clock().get();
	}

	phdTicker(double _duration, double _executions, double _velocity, phdTheClock * _clock = NULL) {
		setDuration(_duration);     // to force calculation of stepSize_
		setExecutions(_executions); // to keep executions != 0.0
		velocity_ = _velocity;
		curExec_ = curPos_ = curVel_ = curClock_ = 0.0;
		dirExec_ = 1.0;
		blendTo_ = NULL;
		theClock_ = _clock;
		loopMode_ = plmtNormal;
		boundsMode_ = plbmLoopForever;
		status_ = pastStopped;
		active_ = false;
		lastUpdate_ = clock().get();
	}

	ofEvent<phdTickerEventData> & events() { return events_; }

	inline double curExec() { return curExec_; }
	inline double curPos()  { return curPos_; }
	inline double curVal()  { return phdMainCurveSystem().get(curveIndex_).calculate(curPos_); }

	inline double duration()   { return duration_; }
	inline double velocity()   { return velocity_; }
	inline double executions() { return executions_; }

	inline void setDuration(double _value)   { duration_ = _value; stepSize_ = ( duration_ == 0.0 ? 1.0 : 1.0 / duration_ ); }
	inline void setExecutions(double _value) { executions_ = ( _value == 0.0 ? 1.0 : _value ); }
	inline void setVelocity(double _value)   { velocity_ = curVel_ = _value; }
	inline void setCurExec(double _value)    { curExec_ = MAX(0.0, MIN(executions_, _value)); setPosition(curExec_); }

	inline void setPosition(double _value)   { // works like a move to
		double _desloc = _value - curPos_;
		update(_desloc);
	}

	inline void setLoopMode(phdLoopModeType _value, double _dir) {
		loopMode_ = _value;
		dirExec_  = _dir; // indicates direction when starting execution (pingpong inverts it every cycle --- normal keep it constant)
		update(0);
	}
	inline void setLoopBoundsMode(phdLoopBoundsMode _value) { boundsMode_ = _value; update(0); }

	inline void setup(double _duration, double _executions, double _velocity, phdLoopModeType _loopMode, phdLoopBoundsMode _boundsMode, bool _startForward = true) {
		setDuration(_duration);     // force calculation of stepSize_
		setExecutions(_executions); // keep executions != 0.0
		setVelocity(_velocity);     // force change curVel_
		boundsMode_ = _boundsMode;
		loopMode_ = _loopMode;
		dirExec_ = (_startForward ? 1.0 : -1.0); // direction of execution
		curPos_ = curExec_ = curClock_ = 0.0;
		update(0);
	}

	inline ofVec2f & getValue() { return curValue_ * curVal(); }
	inline void setValue(const ofVec2f _value) { curValue_ = _value; }
	inline void setValue(const ofVec2f _value, phdTicker & _animation) { }

	inline void pause() { status_ = pastPaused; curVel_ = 0.0; notify(paetOnPause); }

	// stop reset everything --> back to beginning
	inline void stop()  {
		status_ = pastStopped;
		curVel_ = curPos_ = curExec_ = curClock_ = 0.0;
		lastUpdate_ = clock().get();
		notify(paetOnStop);
	}

	inline void play()  {
		status_ = pastPlaying; 
		curVel_ = velocity_;
		lastUpdate_ = clock().get();
		notify(paetOnPlay);
	}

	bool isActive()  {
		if( (curClock_ >= 0 && curClock_ <= executions_) ) return true; // inside range
		if( (loopMode_ != plmtNone) && (boundsMode_ == plbmLoopForever) ) return true; // looping forever
		if( (loopMode_ != plmtNone) && (boundsMode_ == plbmLoopOver)   && (curClock_ >= 0.0) ) return true; // looping after range
		if( (loopMode_ != plmtNone) && (boundsMode_ == plbmLoopBelow)  && (curClock_ <= executions_) ) return true; // looping before range
		return false;
	}
	bool isPaused()  { return (status_ == pastPaused); }
	bool isStopped() { return (status_ == pastStopped); }
	bool isPlaying() { return (status_ == pastPlaying && isActive()); }

	void update(double _deltaClock) {

		double _oldClock = curClock_;

		curClock_ += _deltaClock * stepSize_ * curVel_;// * dirExec_;

		lastUpdate_ = clock().get(); // save current clock position

		bool _wasActive = active_;

		active_ = isActive();

		if(!_wasActive &&  active_) notify(paetOnActivate);

		if(active_) {

			double _delta = (curClock_ - _oldClock) * dirExec_; // adjusted by direction of execution

			curExec_ += _delta; // curExec is a clock shrinked inside range
			//curPos_  += _delta;

			double _startPos = 0.0;
			double _endPos = 1.0;

			bool _overBounds = false;
			/*while(curPos_ < _startPos || curPos_ > _endPos) { // out of execution range
				switch(loopMode_) {
					case plmtNormal   :	if(curPos_ < _startPos) {
											curPos_ = _endPos - ( _startPos - (curPos_) ); _overBounds = true;
										} else if(curPos_ > _endPos) {
											curPos_ = _startPos + ( curPos_ - (_endPos) ); _overBounds = true;
										} break;

					case plmtPingPong :	if(curPos_ < _startPos) {
											curPos_ = _startPos + ( _startPos - (curPos_) );
											if(_delta < 0.0) dirExec_ *= -1.0; // execution going backward
											_overBounds = true;
										} else if(curPos_ > _endPos) {
											curPos_ = _endPos - ( curPos_ - (_endPos) );
											if(_delta > 0.0) dirExec_ *= -1.0; // execution going forward
											_overBounds = true;
										} break;

					default           : curPos_ = MAX(0.0, MIN(_endPos, curPos_)); break;
				}
			}*/

			_startPos = 0.0;
			_endPos = _startPos + executions_;
			while(curExec_ < _startPos || curExec_ > _endPos) { // out of execution range
				switch(loopMode_) {
					case plmtNormal   :	if(curExec_ < _startPos) {
											curExec_ = _endPos - ( _startPos - (curExec_) ); _overBounds = true;
										} else if(curExec_ > _endPos) {
											curExec_ = _startPos + ( curExec_ - (_endPos) ); _overBounds = true;
										} break;

					case plmtPingPong :	if(curExec_ < _startPos) {
											curExec_ = _startPos + ( _startPos - (curExec_) );
											if(_delta < 0.0) dirExec_ *= -1.0; // execution going backward
											_overBounds = true;
										} else if(curExec_ > _endPos) {
											curExec_ = _endPos - ( curExec_ - (_endPos) );
											if(_delta > 0.0) dirExec_ *= -1.0; // execution going forward
											_overBounds = true;
										} break;

					default           : curExec_ = MAX(0.0, MIN(_endPos, curExec_)); break;
				}
			}

			curPos_ = FRACLONG(curExec_);  // normalized position
			if(curPos_ == 0.0 && curExec_ != 0.0) curPos_ = 1.0; // at end of execution position = 1.0

			if(_overBounds) notify(paetOnCycleCompleted);

		} else {
			if(_wasActive && !active_) {
				curExec_ = MAX(0.0, MIN(executions_, curClock_));
				curPos_ = FRACLONG(curExec_);  // normalized position
				if(curPos_ == 0.0 && curExec_ != 0.0) curPos_ = 1.0; // at end of execution position = 1.0
				notify(paetOnDeactivate);
			}
		}
	}

	inline void update() { update(clock().get() - lastUpdate_); }

	void saveToXML(ofXml & _xml) {
		_xml.addChild("PHDTICKER"); _xml.setTo("PHDTICKER");
			_xml.setValue("DURATION",   ofToString(duration()));
			_xml.setValue("EXECUTIONS", ofToString(executions()));
			_xml.setValue("VELOCITY",   ofToString(velocity()));
			_xml.setValue("BOUNDSMODE", ofToString((int)boundsMode_));
			_xml.setValue("LOOPMODE",   ofToString((int)loopMode_));
			_xml.setValue("DIREXEC",    ofToString((int)dirExec_));
		_xml.setToParent(1); // PHDTICKER
	}

	void loadFromXML(ofXml & _xml) {
		_xml.setTo("PHDTICKER");
			setExecutions(ofToFloat(_xml.getAttribute("EXECUTIONS")));
			setVelocity(ofToFloat(_xml.getAttribute("VELOCITY")));
			setLoopBoundsMode((phdLoopBoundsMode) ofToInt(_xml.getAttribute("BOUNDSMODE")));
			dirExec_ = ofToInt(_xml.getAttribute("DIREXEC"));
			setLoopMode((phdLoopModeType) ofToInt(_xml.getAttribute("LOOPMODE")), dirExec_);
		_xml.setToParent(1); // PHDTICKER
	}

	void draw(float _x, float _y, float _w, float _h) {
		ofPushStyle();
			ofFill();

			ofColor & _c = ofGetStyle().color;
			ofSetColor((192-_c.r)*0.25, (192-_c.g)*0.25, (192-_c.b)*0.25, 255);
			ofRect(_x,_y,_w,_h);

			ofSetColor(_c.r,_c.g,_c.b,255);
			ofRect(_x,_y+_h*0.8,_w,_h*0.2);

			float _ecx = _x + _w*0.5; float _eRay = MIN(_w,_h*0.75);
			float _ecy = _y + _h*0.4;
			ofEllipse(_ecx, _ecy, _eRay, _eRay);

			ofSetColor(255-_c.r, 255-_c.g, 255-_c.b, 255);
			ofSetColor(0,0,0,255);
			ofLine(_ecx, _ecy, _ecx + _eRay*0.5*cos(TWO_PI*curPos_), _ecy + _eRay*0.5*sin(TWO_PI*curPos_));
			ofDrawBitmapString(ofToString(curPos_, 3), _ecx - 30, _ecy);
			ofDrawBitmapString(ofToString(velocity_,2) + "x", _ecx - 30, _ecy + 16);

			ofSetColor(0,0,0,255);
			ofLine(_x + _w*curExec_/executions_, _y+_h*0.8+2, _x + _w*curExec_/executions_, _y+_h-2);
			ofDrawBitmapString(ofToString(curExec_, 3), _ecx - 30, _y+_h*0.8+15);
		ofPopStyle();
	}

	void draw(float _x, float _y, float _w, float _h, double _from, double _to) {
		ofPushStyle();
			ofFill();
			//ofSetColor(192,90,32,255);
			ofColor _cA = ofGetStyle().color; // track color
			ofColor _cB = ofColor(ABS(192-_cA.r)*0.75, ABS(192-_cA.g)*0.75, ABS(192-_cA.b)*0.75); // executions color
			ofColor _cC = ofColor(ABS(192-_cA.r)*0.50, ABS(192-_cA.g)*0.50, ABS(192-_cA.b)*0.55); // position color
			ofColor _cD = ofColor(ABS(192-_cA.r)*0.25, ABS(192-_cA.g)*0.25, ABS(192-_cA.b)*0.25); // cursor color

			ofSetColor(_cA, 255);
			ofRect(_x,_y,_w,_h);

			double _min = MIN(curClock_, _from);
			double _max = MAX(curClock_, _to);

			float _x0 = ofMap(0.0,       _min, _max, _x, _x+_w, true);
			float _xc = ofMap(curClock_, _min, _max, _x, _x+_w, true);
			float _xe = ofMap(curExec_,  _min, _max, _x, _x+_w, true);
			float _xp = ofMap(curPos_,   _min, _max, _x, _x+_w, true);

			float _vw = ofMap(1.0,  _min, _max, _x, _x+_w, true);

			ofSetColor(_cB, 160); // curPos 0.0 - 1.0
			ofRect(_x0, _y+_h*0.0, _vw-_x0, _h*0.20);

			if(loopMode_ != plmtNone) {
				ofSetColor(_cC,160); // executions
				_vw = ofMap(executions_, _min, _max, _x, _x+_w, true);
				ofRect(_x0, _y+_h*0.2, _vw-_x0, _h*0.8);
				ofSetColor(_cD,160); // curExec
				switch(boundsMode_) {
					case plbmLoopBelow    :
						ofRect( _x, _y+_h*0.2, _x0-_x, _h*0.8);
						ofSetColor(_cC,160); // executions
						ofRect(MIN(_x0,_xc), _y+_h*0.6, _vw-MIN(_x0,_xc), _h*0.4);
						break;
					case plbmLoopOver     :
						ofRect(_x0, _y+_h*0.2, _w-(_x0-_x), _h*0.8);
						ofSetColor(_cC,160); // executions
						ofRect(_x0, _y+_h*0.6, MAX(_vw,_xc)-_x0, _h*0.4);
						break;
					case plbmLoopForever  :
						ofRect( _x, _y+_h*0.2, _w, _h*0.8);
						ofSetColor(_cC,160); // executions
						if(_xc < _x0) {
							ofRect(MIN(_x0,_xc), _y+_h*0.6, _vw-MIN(_x0,_xc), _h*0.4);
						} else if(_xc > _vw) {
							ofRect(_x0, _y+_h*0.6, MAX(_vw,_xc)-_x0, _h*0.4);
						}
						break;
				}			
			} else {
				ofSetColor(_cC,160); // executions
				_vw = ofMap(executions_, _min, _max, _x, _x+_w, true);
				ofRect(_x0, _y+_h*0.2, _vw-_x0, _h*0.8);
			}

			ofSetColor(0,0,0,255); // cursors are blacks
			ofLine(_xc, _y,         _xc, _y+_h);
			ofLine(_xp, _y+_h*0.0, _xp, _y+_h*0.2);
			ofLine(_xe, _y+_h*0.2, _xe, _y+_h);
		ofPopStyle();
	}

	void drawStatus(float _x, float _y) {
		string _s = "";
		_s += "duration..:" + ofToString(duration_,     10) + "\n";
		_s += "velocity..:" + ofToString(velocity_,     10) + "\n";
		_s += "executions:" + ofToString(executions_,   10) + "\n";
		_s += "stepSize..:" + ofToString(stepSize_,     10) + "\n";
		_s += "curExec...:" + ofToString(curExec_,      10) + "\n";
		_s += "curPos....:" + ofToString(curPos_,       10) + "\n";
		_s += "curVel....:" + ofToString(curVel_,       10) + "\n";
		_s += "lastUpdate:" + ofToString(lastUpdate_,   10) + "\n";
		_s += "theClock..:" + ofToString(clock().get(), 10) + "\n";
		_s += "curClock..:" + ofToString(curClock_,     10) + "\n";
		_s += "active....:" + (active_     ? (string)"TRUE\n" : (string)"FALSE\n");
		_s += "isActive..:" + (isActive()  ? (string)"TRUE\n" : (string)"FALSE\n");
		_s += "isPaused..:" + (isPaused()  ? (string)"TRUE\n" : (string)"FALSE\n");
		_s += "isStopped.:" + (isStopped() ? (string)"TRUE\n" : (string)"FALSE\n");
		_s += "isPlaying.:" + (isPlaying() ? (string)"TRUE\n" : (string)"FALSE\n");
		if(loopMode_ != plmtNone) {
			if(loopMode_ == plmtNormal) _s +=        "loopMode..:plmtNormal\n"; else _s += "loopMode..:plmtPingPong\n";
			if(boundsMode_ == plbmLoopBelow)   _s += "boundsMode:<--- loop below =|\n";
			if(boundsMode_ == plbmLoopOver)    _s += "boundsMode:|= loop over ---->\n";
			if(boundsMode_ == plbmLoopForever) _s += "boundsMode:<--loop forever-->\n";
		} else {
			_s += "loopMode..:plmtNone\n";
		}
		ofDrawBitmapString(_s, _x, _y);
	}
};
/*
static const double _PI= 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
static const double _TWO_PI= 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696;

// Floating-point modulo
// The result (the remainder) has same sign as the divisor.
// Similar to matlab's mod(); Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
//template<typename T>
double Mod(double x, double y) {

	static_assert(!std::numeric_limits<double>::is_exact , "Mod: floating-point type expected");

    if (0. == y) return x;

    double m = x - y * floor(x/y);

    // handle boundary cases resulted from floating-point cut off:
    if (y > 0) {             // modulo range: [0..y)
        if (m>=y) return 0;  // Mod(-1e-16             , 360.    ): m= 360.

        if (m<0 ) {
            if (y+m == y)
                return 0  ; // just in case...
            else
                return y+m; // Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14 
        }
    } else {                // modulo range: (y..0]
        if (m<=y) return 0; // Mod(1e-16              , -360.   ): m= -360.
        if (m>0 ) {
            if (y+m == y)
                return 0;   // just in case...
            else
                return y+m; // Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14 
        }
    }
    return m;
}
*/