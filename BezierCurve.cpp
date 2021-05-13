#include "BezierCurve.h"

void BezierCurve::init(int g){
	coefficients = std::vector<int>(g+1);
	n = g;
	for (int i = 0; i < coefficients.size(); i++)
	{
		coefficients[i] = factorial(coefficients.size()) /
			(factorial(i) * factorial(coefficients.size() - i));
	}
}
Point BezierCurve::getFunctionValue(double t) {
	double sum = 0;
	for (int i = 0; i < coefficients.size(); i++) {
		sum += coefficients[i] * pow(t, i)
			* pow((1 - t), coefficients.size() - i);
	}
	return Point(t, sum, 0.0);
}

int BezierCurve::factorial(int n) {
	if (n = 1)
		return 1;
	else
		return n * factorial(n - 1);
}

void BezierCurve::regenerateControlPoints()
{
	if (cp.size() <= n) {
		derivatedCp.resize(n);
		secondDerivatedCp.resize(n - 1);
		cp.clear();
		for (size_t i = 0; i < n + 1; i++)
		{
			cp.push_back(randomPoint(-1.0, 1.0));
		}
	}
}
Point BezierCurve::randomPoint(double fMin, double  fMax) {
	double f = (double)rand() / RAND_MAX;
	double x = fMin + f * (fMax - fMin);

	double k = (double)rand() / RAND_MAX;
	double y = fMin + k * (fMax - fMin);

	double h = (double)rand() / RAND_MAX;
	double z = fMin + h * (fMax - fMin);

	return Point{ x, y, 0.0 };
}

double BezierCurve::bernstein(size_t i, size_t n, double u) const
{
	DoubleVector tmp(n + 1, 0.0);
	tmp[n - i] = 1.0;
	double u1 = 1.0 - u;
	for (size_t k = 1; k <= n; ++k)
		for (size_t j = n; j >= k; --j)
			tmp[j] = tmp[j - 1] * u + tmp[j] * u1;
	return tmp[n];
}

Point BezierCurve::evaluateOneByOne(double u) const
{
	Point p(0.0, 0.0, 0.0);
	for (size_t k = 0; k <= n; ++k)
		p += cp[k] * bernstein(k, n, u);
	return p;
}


// [3] Az Üsszes Bernstein polinom kiÊrtÊkelÊse egyszerre

void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector &coeff) const
{
	coeff.clear(); coeff.reserve(n + 1);
	coeff.push_back(1.0);
	double u1 = 1.0 - u;
	for (size_t j = 1; j <= n; ++j) {
		double saved = 0.0;
		for (size_t k = 0; k < j; ++k) {
			double tmp = coeff[k];
			coeff[k] = saved + tmp * u1;
			saved = tmp * u;
		}
		coeff.push_back(saved);
	}
}

Point BezierCurve::evaluate(double u) const
{
	DoubleVector coeff; 
	bernsteinAll(n, u, coeff);
	Point p(0.0, 0.0, 0.0);
	for (size_t k = 0; k <= n; ++k)
		p += cp[k] * coeff[k];
	return p;
}

Point BezierCurve::evaluateWithCachedCofficients(const DoubleVector &coeff) const
{
	Point p(0.0, 0.0, 0.0);
	for (size_t k = 0; k <= n; ++k)
		p += cp[k] * coeff[k];
	return p;
}


// [5] deCasteljau algoritmus

Point BezierCurve::evaluateByDeCasteljau(double u) const
{
	PointVector tmp = cp;
	double u1 = 1.0 - u;
	for (size_t k = 1; k <= n; ++k)
		for (size_t i = 0; i <= n - k; ++i)
			tmp[i] = tmp[i] * u1 + tmp[i + 1] * u;
	return tmp[0];
}


// [6] BÊzier derivåltak kontrollpontjainak kiszåmítåsa

// FeltÊtelezi, hogy d <= n
void BezierCurve::derivativeControlPoints(size_t d, PointMatrix &dcp) const
{
	dcp.clear(); dcp.resize(d + 1);
	dcp[0] = cp;
	for (size_t k = 1; k <= d; ++k) {
		size_t tmp = n - k + 1;
		dcp[k].reserve(tmp);
		for (size_t i = 0; i <= n - k; ++i)
			dcp[k].push_back((dcp[k - 1][i + 1] - dcp[k - 1][i]) * tmp);
	}
}

void BezierCurve::getFristDerivateControlPoints(PointVector& derivateCp) {
	derivateCp.clear();
	for (size_t i = 0; i < n; i++)
	{
		derivateCp.push_back((cp[i + 1] - cp[i]));
	}
}
void BezierCurve::getSecondDerivateControlPoints(PointVector& secondDerivateCp) {
	secondDerivateCp.clear();
	for (size_t i = 0; i < n - 1; i++)
	{
		secondDerivateCp.push_back(derivatedCp[i + 1] - derivatedCp[i]);
	}
}

Point BezierCurve::getSecondDerivatedValue(double u) {
	if (derivatedCp.size() == 0)
		getFristDerivateControlPoints(derivatedCp);
	getSecondDerivateControlPoints(secondDerivatedCp);
	DoubleVector coeff;
	bernsteinAll(n - 2, u, coeff);
	Point p(0.0, 0.0, 0.0);
	for (size_t k = 0; k <= n - 2; ++k)
		p += secondDerivatedCp[k] * coeff[k] * n;
	return p;
}

double BezierCurve::curvature(double u)
{
	double curve;
	curve = vectorLength(crossProduct(getFirstDerivatedValue(u), getSecondDerivatedValue(u))) /
		pow(vectorLength(getFirstDerivatedValue(u)), 3);
	return curve;
}

double BezierCurve::curveIntegralBaseFunction(double u, void * data)
{
	BezierCurve* b = (BezierCurve*)data;
	double result;
	result = pow(b->curvature(u), 2) * b->vectorLength(b->getFirstDerivatedValue(u));
	return result;
}

autodiff::dual BezierCurve::sumCurvature()
{
	autodiff::dual cur;
	for (size_t i = 2; i < 128; i++)
	{
		cur = gauss_legendre(i, curveIntegralBaseFunction, this, 0, 1);
	}
	return cur;
}
autodiff::dual BezierCurve::energyFunction(const autodiff::VectorXdual & x) {
	BezierCurve c;
	for (size_t i = 0; i < (x.size() / 3); i++)
	{
		c.cp.push_back(Point((double)x[3 * i], (double)x[3 * i + 1], 
			(double)x[3 * i + 2]));
	}
	c.n = c.cp.size() - 1;
	return c.sumCurvature();
}

VectorVector BezierCurve::getGradientVector()
{
	autodiff::VectorXdual x((cp.size() * 3));
	VectorVector vv;
	for (size_t i = 0; i < cp.size(); i++)
	{
		x[3 * i] = cp[i].x;
		x[3 * i + 1] = cp[i].y;
		x[3 * i + 2] = cp[i].z;
	}
	autodiff::dual u;
	Eigen::VectorXd g = autodiff::forward::gradient(energyFunction, autodiff::wrt(x), 
		autodiff::forward::at(x), u);
	for (size_t i = 0; i < cp.size(); i++)
	{
		Vector v = Vector(g[3 * i], g[3 * i + 1], g[3 * i + 2]);
		vv.push_back(v);
	}
	return vv;
}

void BezierCurve::gradientDescend()
{
	VectorVector vv(getGradientVector());
	for (size_t i = 1; i < cp.size() - 1; i++)
	{
		cp[i] -= vv[i];
	}
}

//TODO: ide lehet kell majd az n *-resz de meg nem biztos!!!
Point BezierCurve::getFirstDerivatedValue(double u) {
	getFristDerivateControlPoints(derivatedCp);
	DoubleVector coeff;
	bernsteinAll(n - 1, u, coeff);
	Point p(0.0, 0.0, 0.0);
	for (size_t k = 0; k <= n - 1; ++k)
		p +=  derivatedCp[k] * coeff[k] * n;
	return p;
}
// [7] Az Üsszes Bernstein polinom (minden kisebb fokszåmra is)

void BezierCurve::bernsteinAll(size_t n, double u, DoubleMatrix &coeff) const
{
	coeff.clear(); coeff.resize(n + 1);
	coeff[0].push_back(1.0);
	double u1 = 1.0 - u;
	for (size_t j = 1; j <= n; ++j) {
		coeff[j].reserve(j + 1);
		double saved = 0.0;
		for (size_t k = 0; k < j; ++k) {
			double tmp = coeff[j - 1][k];     // ... = coeff[k]  helyett
			coeff[j].push_back(saved + tmp * u1); // coeff[k] = ...  helyett
			saved = tmp * u;
		}
		coeff[j].push_back(saved);
	}
}

size_t BezierCurve::getIndexOfContorlPoint(Point point)
{
	for(size_t i = 0; i < cp.size(); i++)
	{
		if (cp[i].x == point.x && cp[i].y == point.y && cp[i].z == point.z)
			return i;
	}
}


// [8] Derivåltszåmítås

Point BezierCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
	size_t du = std::min(d, n);
	der.clear(); der.reserve(d + 1);
	DoubleMatrix coeff; bernsteinAll(n, u, coeff);
	PointMatrix dcp; derivativeControlPoints(du, dcp);
	for (size_t k = 0; k <= du; ++k) {
		der.emplace_back(0.0, 0.0, 0.0);
		for (size_t j = 0; j <= n - k; ++j)
			der[k] += dcp[k][j] * coeff[n - k][j];
	}
	for (size_t k = n + 1; k <= d; ++k)
		der.emplace_back(0.0, 0.0, 0.0);
	return der[0];
}

double BezierCurve::arcLengthByFractions() const
{
	float length = 0;
	for (double i = 0.0; i < 1.0; i += 0.01) {
		Point p1 = evaluateByDeCasteljau(i);
		Point p2 = evaluateByDeCasteljau(i + 0.01);
		Vector tmp = p1 - p2;
		float fractionLength = vectorLength(tmp);
		length += fractionLength;
	}
	return length;
}

double BezierCurve::derivatedValueLength(double x, void* data) {
	BezierCurve* tmp = (BezierCurve*)data;
	return tmp->vectorLength(tmp->getFirstDerivatedValue(x));
}

Point BezierCurve::crossProduct(Point p1, Point p2)
{
	Point result;
	result.x = p1.y * p2.z - p1.z * p2.y;
	result.y = p1.z * p2.x - p1.x * p2.z;
	result.z = p1.x * p2.y - p1.y * p2.x;
	return result;
}

double BezierCurve::arcLengthByNumericalIntegral()
{
	double length;
	for (size_t i = 2; i < 128; i++)
	{
		length = gauss_legendre(i, derivatedValueLength, this, 0, 1);
	}
	return length;
}

double BezierCurve::vectorLength(Vector v) const
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

