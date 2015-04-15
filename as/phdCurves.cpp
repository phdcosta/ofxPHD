#pragma once

#include "phdCurves.h"

//----------------------------------------------------------------------------------------------------------------------------
phdCurveSystem & phdMainCurveSystem() {
	return phdMainCurves;
}

//----------------------------------------------------------------------------------------------------------------------------
phdMultPointCurveCalculator::phdMultPointCurveCalculator(float _startValue, float _endValue, int _curveIndex) {

	pointEditMode = pcaAddReposAll;

	addPoint(0.0, _startValue, _curveIndex); // first point is always at 0.0 pos
	addPoint(1.0, _endValue,   _curveIndex); // last point is always at 1.0 pos
}

void phdMultPointCurveCalculator::setPointPos(unsigned int _index, double _pos) {

	int nPoints = points.size();

	if(nPoints == 0) { // no points yet? add one as last
		addPoint(1.0, 1.0, 0);
		return;
	}

	if(nPoints == 1) { // only one point changes it to position 1.0
		points[0].pos = 1.0;
		return;
	}

	_index -= 1; // _index 0 is reserved for startValue that is always at position 0.0

	if(_index < nPoints) { // index in range?

		if(_index == nPoints-1 && _pos > points[_index-1].pos) { // last point pos must be greater than previous

			if(pointEditMode == pcaAddReposAll) {
				points[_index].pos = _pos;
			} else if(pointEditMode == pcaIgnoreAdd) {
				points[_index].pos = 1.0;
			}

		} else if(_index == 0 && _pos > 0.0 && _pos < points[1].pos) { // first point

			points[_index].pos = 0.0;

		} else if(_pos > points[_index-1].pos && _pos < points[_index+1].pos) {

			points[_index].pos = _pos;

		}
	}
}

void phdMultPointCurveCalculator::setPointValue(unsigned int _index, double _value) {
	if(_index < points.size()) {
		points[_index].value = _value;
	}
}

void phdMultPointCurveCalculator::setSegmentCurveMode(unsigned int _index, int _curveIndex) {
	if(_index < points.size()-1) { // N curves == N segments == N-1 points
		points[_index+1].curve = phdMainCurveSystem().getCurveByIndex(_curveIndex); // changes curve proc on destination point of segment
		if(_index == 0) points[_index].curve = points[_index+1].curve;
	}
}

void phdMultPointCurveCalculator::setSegmentCurveMode(unsigned int _index, string _curveName) {
	if(_index < points.size()-1) { // N curves == N segments == N-1 points
		points[_index+1].curve = phdMainCurveSystem().getCurveByName(_curveName); // changes curve proc on destination point of segment
		if(_index == 0) points[_index].curve = points[_index+1].curve;
	}
}

int phdMultPointCurveCalculator::segmentIndexByPos(double _pos) {
	int nPoints = points.size();
	if(nPoints == 0) return -1;
	int _index = 0;
	while(_index < nPoints-1 && _pos >= points[_index].pos) { _index += 1; }
	return _index;
}

bool phdMultPointCurveCalculator::getSegment(unsigned int _index, phdSegment & _segment) {
	if(_index > 0 && _index < points.size()) {
		_segment.ta = points[_index-1].pos;
		_segment.va = points[_index-1].value;
		_segment.tb = points[_index+0].pos;
		_segment.vb = points[_index+0].value;
		_segment.curve = points[_index+0].curve;
		return true;
	}
	return false;
}

double phdMultPointCurveCalculator::getValueInPos(double _pos) {
	int _index = segmentIndexByPos(_pos);
	phdSegment _seg;
	if(getSegment(_index, _seg)) {
		return getValueInSegmentByPos(_seg, _pos);
	}
	return -1.0;
}

double phdMultPointCurveCalculator::getValueInSegmentByPos(phdSegment & _seg, double _position) {
	double _pos = (_position - _seg.ta) / (_seg.tb - _seg.ta);
	double _var = (_seg.vb - _seg.va);
	return _seg.va + (_var * _seg.curve->calculate(_pos));
}

void phdMultPointCurveCalculator::addPoint(double _pos, double _value, int _curveIndex) {

	int nPoints = points.size();

	if(nPoints == 0) { // no points? always add at first position --> 0.0

		points.push_back(phdPosPoint(0.0, _value, phdMainCurveSystem().getCurveByIndex(_curveIndex)));

	} else if(nPoints == 1) { // only the start point? always add at last position --> 1.0

		points.push_back(phdPosPoint(1.0, _value, phdMainCurveSystem().getCurveByIndex(_curveIndex)));

	} else {
		
		if(_pos >= 1.0) { // tryin to add at last position

			if(pointEditMode == pcaAddReposAll) {

				if(_pos > points[nPoints-1].pos) {
					points.push_back(phdPosPoint(_pos, _value, phdMainCurveSystem().getCurveByIndex(_curveIndex)));
					normalizePoints();
				} else { // _pos == 1.0
					points[nPoints-1].value = _value;
					points[nPoints-1].curve = phdMainCurveSystem().getCurveByIndex(_curveIndex);
				}

			} else if(pointEditMode == pcaIgnoreAdd) {

				points[nPoints-1].value = _value;
				points[nPoints-1].curve = phdMainCurveSystem().getCurveByIndex(_curveIndex);
			}

		} else if(_pos <= 0.0) { // trying to add at first position --> changes start value

			if(pointEditMode == pcaAddReposAll) {

				points.insert(points.begin(), phdPosPoint(_pos, _value, phdMainCurveSystem().getCurveByIndex(_curveIndex)));
				normalizePoints();

			} else if(pointEditMode == pcaIgnoreAdd) {

				points[0].value = _value;
				points[0].curve = phdMainCurveSystem().getCurveByIndex(_curveIndex);
			}

		} else {

			int _index = segmentIndexByPos(_pos); // find point index by position --> return -1 if pos > 1.0
			points.insert(points.begin() + _index, phdPosPoint(_pos, _value, phdMainCurveSystem().getCurveByIndex(_curveIndex)));
		}
	}
}

void phdMultPointCurveCalculator::delPoint(unsigned int _index) {
	if(_index < points.size() && points.size() > 1) { // leave at least two points
		points.erase(points.begin() + _index);
		normalizePoints();
	}
}

void phdMultPointCurveCalculator::normalizePoints() {

	int nPoints = points.size();

	if(nPoints == 0) return;
	if(nPoints == 1) { points[0].pos = 0.0; return; }
	if(nPoints == 2) { points[0].pos = 0.0; points[1].pos = 1.0; return; }

	double _low  = points[0].pos;
	double _high = points[nPoints-1].pos;

	for(int i = 0; i < nPoints; i++) {
		points[i].pos = ofMap(points[i].pos, _low, _high, 0.0, 1.0, true);
	}
}

void phdMultPointCurveCalculator::draw(float x, float y, float w, float h, float res, bool inv) {

	double _sp = 1.0 / (w / res); // line segment each RES px

	ofNoFill();
	ofRect(x, y, w, h);

	glBegin(GL_LINE_STRIP);
		for(int j = 0; j < points.size(); j++) {

			phdSegment _seg;
			if(getSegment(j, _seg)) {

				double _pos = _seg.ta;

				glVertex2d(x + INVIF(inv, w, (w * _seg.ta)), y + h - getValueInSegmentByPos(_seg, _seg.ta) * h);
				_pos += _sp;

				while(_pos + _sp < _seg.tb) {
					glVertex2d(x + INVIF(inv, w, (w * _pos)), y + h - getValueInSegmentByPos(_seg, _pos) * h);
					_pos += _sp;
				}

				glVertex2d(x + INVIF(inv, w, (w * _seg.tb)), y + h - getValueInSegmentByPos(_seg, _seg.tb) * h);
			}
		}
	glEnd();
		
	for(int j = 0; j < points.size(); j++) {
		phdSegment _seg;
		if(getSegment(j, _seg)) {
			if(j == 0) ofCircle(x + INVIF(inv, w, (w * _seg.ta)), y + h - getValueInSegmentByPos(_seg, _seg.ta) * h, 2);
			ofCircle(x + INVIF(inv, w, (w * _seg.tb)), y + h - getValueInSegmentByPos(_seg, _seg.tb) * h, 2);
		}
	}
}