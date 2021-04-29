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
	for (size_t i = 0; i < n; i++)
	{
		derivateCp[i] =  cp[i + 1] - cp[i];
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
	for (double i = 0.0; i <= 1.0; i += 0.01) {
		Point p1 = evaluateByDeCasteljau(i);
		Point p2 = evaluateByDeCasteljau(i + 0.01);
		Vector tmp = p1 - p2;
		float fractionLength = vectorLength(tmp);
		length += fractionLength;
	}
	return length;
}

double BezierCurve::arcLengthByNumericalIntegral()
{
	float length;
	length = (vectorLength(getFirstDerivatedValue((1/2) * (1 / sqrt(3)) + (1 / 2)))
		+ vectorLength(getFirstDerivatedValue((1 / 2) * (-1 / sqrt(3)) + (1 / 2)))) / 2;
	return length;
}

double BezierCurve::vectorLength(Vector v) const
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}
