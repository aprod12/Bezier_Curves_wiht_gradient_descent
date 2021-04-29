#pragma once

// Vector típus a szokåsos mŹveletekkel
#include "vector.hh"
#include <algorithm>

using Point = Vector;
using PointVector = std::vector<Point>;
using VectorVector = std::vector<Vector>;
using DoubleVector = std::vector<double>;
using DoubleMatrix = std::vector<DoubleVector>;
using PointMatrix = std::vector<PointVector>;

struct BezierCurve
{
	size_t n;                     // fokszåm => n + 1 kontrollpont
	PointVector cp;		            // kontrollpontok
	std::vector<int> coefficients;
	PointVector derivatedCp;

	Point getFunctionValue(double t);
	void init(int g);
	int factorial(int n);

	void regenerateControlPoints();
	Point randomPoint(double fMin, double  fMax);

	double bernstein(size_t i, size_t n, double u) const;
	Point evaluateOneByOne(double u) const;
	void bernsteinAll(size_t n, double u, DoubleVector & coeff) const;
	Point evaluate(double u) const;
	Point evaluateWithCachedCofficients(const DoubleVector & coeff) const;
	Point evaluateByDeCasteljau(double u) const;
	void derivativeControlPoints(size_t d, PointMatrix & dcp) const;
	void bernsteinAll(size_t n, double u, DoubleMatrix & coeff) const;
	void getFristDerivateControlPoints(PointVector& derivateCp);
	Point getFirstDerivatedValue(double u);
	size_t getIndexOfContorlPoint(Point p);
	Point derivativesByControlPoints(double u, size_t d, VectorVector & der) const;
	double arcLengthByFractions() const;
	double arcLengthByNumericalIntegral();
	double vectorLength(Vector v) const;
};