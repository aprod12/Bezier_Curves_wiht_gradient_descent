#include "BezierCurve.h"

BezierCurve::BezierCurve()
	:optimizationMode(CURVATURE), lineSearch(false),
	gradientDescendStepSize(0.01), iterationCount(0)
{
	lockedControlPointIndexes = { 0,4 };
	//TODO ez lehet ide nem kell
	//s_LockedControlPointIndexes = lockedControlPointIndexes;
}

BezierCurve::BezierCurve(bool forStatic)
{

}

BezierCurve::BezierCurve(size_t _n)
	: n(_n), optimizationMode(CURVATURE), lineSearch(false),
	gradientDescendStepSize(0.01), iterationCount(0)
{
	lockedControlPointIndexes = { 0,4 };
	//TODO ez lehet ide nem kell
	//s_LockedControlPointIndexes = lockedControlPointIndexes;
}

void BezierCurve::setQuadratureDegree(int n)
{
	switch (n)
	{
	case 20:
		numericalWeights = w20;
		numericalFunctionValues = x20;
		sizeOfNumericalArrays = 10;
		break;
	case 128:
		numericalWeights = w128;
		numericalFunctionValues = x128;
		sizeOfNumericalArrays = 64;
		break;
	case 256:
		numericalWeights = w256;
		numericalFunctionValues = x256;
		sizeOfNumericalArrays = 128;
		break;
	default:
		break;
	}
}

size_t BezierCurve::getDegree()
{
	return n;
}
void BezierCurve::setDegree(size_t _n)
{
	n = _n;
}
void BezierCurve::setControlPoints(PointVector&& _cp)
{
	cp = _cp;
}
PointVector& BezierCurve::getControlpoints() { return cp; }

void BezierCurve::setMode(int mode)
{
	optimizationMode = static_cast<OtimizationMode>(mode);
}

void BezierCurve::setLineSearch(bool ls)
{
	lineSearch = ls;
}

void BezierCurve::setGradientDescendStepSize(double newSize)
{
	gradientDescendStepSize = newSize;
}

double BezierCurve::getGradientDescendStepSize()
{
	return gradientDescendStepSize;
}

void BezierCurve::regenerateControlPoints()
{
	cp.push_back(Point(-0.98, -1.0));
	cp.push_back(Point(-0.5, 0.75));
	cp.push_back(Point(0, 1));
	cp.push_back(Point(0.5, 0.75));
	cp.push_back(Point(1.0, -1.0));
	lockedControlPointIndexes = { 0,4 };

}
void BezierCurve::increaseDegree()
{
	PointVector newCps(cp.size() + 1);
	newCps[0] = cp[0];
	newCps[cp.size()] = cp[cp.size() - 1];
	for (double i = 1; i < newCps.size() - 1; i++)
	{
		newCps[i] = Point((cp[i - 1] * (i / (n + 1))) + (cp[i] * (1 - (i / (n + 1)))));
	}
	lockedControlPointIndexes.erase(std::find(lockedControlPointIndexes.begin(),
		lockedControlPointIndexes.end(), n));
	cp = newCps;
	n += 1;
	lockedControlPointIndexes.push_back(n);
}

autodiff::VectorXdual BezierCurve::getNormalVector(double u)
{
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	autodiff::VectorXdual x(3 * (cp.size() - lockedControlPointIndexes.size()));
	autodiff::VectorXdual p(lockedControlPointIndexes.size() * 3);
	controlPointsToAutodiffArguments(x, p);

	autodiff::VectorXdual res(3);

	res = division(crossProduct(firstDerivate(x, p, u),
		crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u))),
		vectorLength(crossProduct(firstDerivate(x, p, u),
			crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u)))));
	return res;
}

autodiff::dual BezierCurve::getCurveSize(double u)
{
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	autodiff::VectorXdual x(3 * (cp.size() - lockedControlPointIndexes.size()));
	autodiff::VectorXdual p(lockedControlPointIndexes.size() * 3);
	controlPointsToAutodiffArguments(x, p);

	autodiff::dual res;

	res = vectorLength(crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u))) /
		pow(vectorLength(firstDerivate(x, p, u)), 3);
	return res * 0.5;
}

PointVector BezierCurve::fittingPoints(double scale)
{
	PointVector res;
	scaleNumberForFitting = scale;
	for (double i = 0.0; i < 1.01; i += scale)
	{
		double f = (double)rand() / RAND_MAX;
		double random = -0.2 + f * (0.2 + 0.2);
		Point p = evaluateByDeCasteljau(i);
		autodiff::VectorXdual normalVector = getNormalVector(i);
		Point p2(p.x + (normalVector[0].val * random),
			p.y + (normalVector[1].val * random),
			p.z + (normalVector[2].val * random));
		res.push_back(p2);
	}
	return res;
}

PointVector BezierCurve::genereateRandomReferencePoints(const double scale) const
{
	return PointVector();
}

auto BezierCurve::optimizationChoosing()
{
	auto fv = BezierCurve::curveEnergyFunction;
	switch (optimizationMode)
	{
	case CURVATURE:
		fv = curveEnergyFunction;
		break;
	case ARCHLENGTH:
		fv = arcLengthEnergyFunction;
		break;
	case TORSION:
		fv = torsionEnergyFunction;
		break;
	case FIRSTDERIVATE:
		fv = firstDerivateEnergyFunction;
		break;
	case SECONDDERIVATE:
		fv = secondDerivateEnergyFunction;
		break;
	case THIRDDERIVATE:
		fv = thirdDerivateEnergyFunction;
		break;
	case CURVEVARIANT:
		fv = curveVariantEnergyFunction;
		break;
	case FITTING:
		fv = fittingEnergyFunction;
		break;
	case FITTINGNORMA:
		fv = fittingNormaEnergyFunction;
		break;
	default:
		break;
	}
	return fv;
}

void BezierCurve::getInformations(double & curvatureSumary, double & fractionLength, double & archIntegral, double & error, double & torsion, double & firstDerivateEnergy, double & secondDerivateEnergy, double & thirdDerivateEnergy, double & curvatureVariant)
{
	autodiff::VectorXdual x(3 * (cp.size() - lockedControlPointIndexes.size()));
	autodiff::VectorXdual p(3 * lockedControlPointIndexes.size());
	controlPointsToAutodiffArguments(x, p);
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	curvatureSumary = curveEnergyFunction(x, p).val;
	fractionLength = arcLengthByFractions();
	archIntegral = arcLengthEnergyFunction(x, p).val;
	error = abs(fractionLength - archIntegral);
	torsion = torsionEnergyFunction(x, p).val;
	firstDerivateEnergy = firstDerivateEnergyFunction(x, p).val;
	secondDerivateEnergy = secondDerivateEnergyFunction(x, p).val;
	thirdDerivateEnergy = thirdDerivateEnergyFunction(x, p).val;
	curvatureVariant = curveVariantEnergyFunction(x, p).val;

}

double BezierCurve::getFittingEnergy()
{
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	autodiff::VectorXdual x(3 * (getControlpoints().size() - lockedControlPointIndexes.size()));
	autodiff::VectorXdual p(3 * lockedControlPointIndexes.size());
	controlPointsToAutodiffArguments(x, p);
	return fittingEnergyFunction(x, p).val;
}

void BezierCurve::makeFittingCurve()
{
	optimizationMode = FITTING;
	cp.clear();
	n = 4;
	cp.push_back(Point(-1.0, 1.0));
	cp.push_back(Point(-0.5, 1.0));
	cp.push_back(Point(0.0, 1.0));
	cp.push_back(Point(0.5, 1.0));
	cp.push_back(Point(1.0, 1.0));
	lockedControlPointIndexes.clear();
}

size_t n_choose_k(size_t n, size_t k) {
	// Base Cases
	if (k > n)
		return 0;
	if (k == 0 || k == n)
		return 1;

	// Recur
	return n_choose_k(n - 1, k - 1)
		+ n_choose_k(n - 1, k);
}

autodiff::dual BezierCurve::fittingEnergyFunction(const autodiff::VectorXdual & x, const autodiff::VectorXdual & p)
{
	BezierCurve c(false);
	autodiff::dual res = 0;
	for (double i = 0; i < 1.01; i += scaleNumberForFitting)
	{
		res += pow(c.vectorLength(c.subtraction(c.evaluateAutodiff(x, p, i), getReferencePoint(i))), 2);
	}
	return res;
}

autodiff::dual BezierCurve::fittingNormaEnergyFunction(const autodiff::VectorXdual & x, const autodiff::VectorXdual & p)
{
	BezierCurve c(false);
	autodiff::dual res = 0;
	for (double i = 0; i < 1.01; i += scaleNumberForFitting)
	{
		res += c.vectorLength(c.subtraction(c.evaluateAutodiff(x, p, i), getReferencePoint(i)));
	}
	return res;
}

autodiff::VectorXdual BezierCurve::getReferencePoint(double index)
{
	int ix = 0;
	autodiff::VectorXdual res(3);
	for (double i = 0; i < 1.01; i += scaleNumberForFitting)
	{
		if (i == index)
		{
			res[0] = referencePoints[ix].x;
			res[1] = referencePoints[ix].y;
			res[2] = referencePoints[ix].z;
			return res;
		}
		ix++;
	}
	return res;
}

autodiff::VectorXdual BezierCurve::evaluateAutodiff(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u)
{
	autodiff::VectorXdual res = autodiff::VectorXdual::Zero(3);
	int degree = ((x.size() + p.size()) / 3) - 1;
	autodiff::VectorXdual ctr1(3);
	size_t x_index = 0;
	size_t p_index = 0;
	autodiff::VectorXdual tmp(3 * (degree + 1));
	for (int i = 0; i < degree + 1; i++)
	{
		if (isStaticLocked(i))
		{
			tmp[3 * i] = p[3 * p_index];
			tmp[3 * i + 1] = p[3 * p_index + 1];
			tmp[3 * i + 2] = p[3 * p_index + 2];
			p_index++;
		}
		else
		{
			tmp[3 * i] = x[3 * x_index];
			tmp[3 * i + 1] = x[3 * x_index + 1];
			tmp[3 * i + 2] = x[3 * x_index + 2];
			x_index++;
		}

	}
	for (size_t i = 0; i < degree + 1; i++)
	{
		ctr1[0] = tmp[3 * i];
		ctr1[1] = tmp[3 * i + 1];
		ctr1[2] = tmp[3 * i + 2];
		res += n_choose_k(degree, i) * pow(u, i) * pow(1 - u, degree - i) * ctr1;
	}
	return res;
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


bool BezierCurve::isStaticLocked(int index)
{
	return (std::find(s_LockedControlPointIndexes.begin(), s_LockedControlPointIndexes.end(), index)
		!= s_LockedControlPointIndexes.end());
}

bool BezierCurve::isItLocked(int index)
{
	return (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), index)
		!= lockedControlPointIndexes.end());
}
void BezierCurve::lockControlPoint(int index)
{
	if (!isItLocked(index))
	{
		lockedControlPointIndexes.push_back(index);
	}
}
void BezierCurve::unlockControlPoint(int index)
{
	if (isItLocked(index))
	{
		lockedControlPointIndexes.erase(std::find(lockedControlPointIndexes.begin(),
			lockedControlPointIndexes.end(), index));
	}
}
#pragma region autodiff

autodiff::dual BezierCurve::combiningEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual & p)
{
	BezierCurve c(false);
	autodiff::dual res = 0;
	float weight = 0.0;
	res = (1 - weight) * curveEnergyFunction(x, p) + weight * arcLengthEnergyFunction(x, p);
	return res;

}

autodiff::VectorXdual BezierCurve::firstDerivate(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u) {
	autodiff::VectorXdual der = autodiff::VectorXdual::Zero(3);
	int degree = ((x.size() + p.size()) / 3) - 1;
	autodiff::VectorXdual ctr1(3), ctr2(3);
	size_t x_index = 0;
	size_t p_index = 0;
	autodiff::VectorXdual tmp(3 * (degree + 1));
	for (int i = 0; i < degree + 1; i++)
	{
		if (isStaticLocked(i))
		{
			tmp[3 * i] = p[3 * p_index];
			tmp[3 * i + 1] = p[3 * p_index + 1];
			tmp[3 * i + 2] = p[3 * p_index + 2];
			p_index++;
		}
		else
		{
			tmp[3 * i] = x[3 * x_index];
			tmp[3 * i + 1] = x[3 * x_index + 1];
			tmp[3 * i + 2] = x[3 * x_index + 2];
			x_index++;
		}

	}
	for (size_t i = 0; i < degree; i++)
	{
		ctr1[0] = tmp[3 * i];
		ctr1[1] = tmp[3 * i + 1];
		ctr1[2] = tmp[3 * i + 2];
		ctr2[0] = tmp[3 * (i + 1)];
		ctr2[1] = tmp[3 * (i + 1) + 1];
		ctr2[2] = tmp[3 * (i + 1) + 2];
		der += degree * n_choose_k(degree - 1, i) * pow(1 - u, degree - (i)-1) *
			pow(u, i)*(ctr2 - ctr1);
	}
	return der;
}
autodiff::VectorXdual BezierCurve::thirdDerivate(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u) {
	autodiff::VectorXdual der = autodiff::VectorXdual::Zero(3);
	DoubleVector coeff;
	int degree = ((x.size() + p.size()) / 3) - 1;
	size_t i = 0;
	autodiff::VectorXdual ctr1(3), ctr2(3), ctr3(3), ctr4(3), ctr5(3);
	size_t x_index = 0;
	size_t p_index = 0;
	autodiff::VectorXdual tmp(3 * (degree + 1));
	for (int i = 0; i < degree + 1; i++)
	{
		if (isStaticLocked(i))
		{
			tmp[3 * i] = p[3 * p_index];
			tmp[3 * i + 1] = p[3 * p_index + 1];
			tmp[3 * i + 2] = p[3 * p_index + 2];
			p_index++;
		}
		else
		{
			tmp[3 * i] = x[3 * x_index];
			tmp[3 * i + 1] = x[3 * x_index + 1];
			tmp[3 * i + 2] = x[3 * x_index + 2];
			x_index++;
		}

	}
	for (size_t i = 0; i < degree - 2; i++)
	{
		ctr1[0] = tmp[3 * i];
		ctr1[1] = tmp[3 * i + 1];
		ctr1[2] = tmp[3 * i + 2];
		ctr2[0] = tmp[3 * (i + 1)];
		ctr2[1] = tmp[3 * (i + 1) + 1];
		ctr2[2] = tmp[3 * (i + 1) + 2];
		ctr3[0] = tmp[3 * (i + 2)];
		ctr3[1] = tmp[3 * (i + 2) + 1];
		ctr3[2] = tmp[3 * (i + 2) + 2];
		ctr4[0] = tmp[3 * (i + 3)];
		ctr4[1] = tmp[3 * (i + 3) + 1];
		ctr4[2] = tmp[3 * (i + 3) + 2];
		der += degree * (degree - 1) * (degree - 2) * n_choose_k(degree - 3, i) * pow(u, (i))
			* pow(1 - u, degree - (i)-3) *
			((((ctr4 - ctr3)) - (ctr3 - ctr2)) - ((ctr3 - ctr2) - (ctr2 - ctr1)));
	}
	return der;
}

autodiff::VectorXdual BezierCurve::secondDerivate(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u) {
	autodiff::VectorXdual der(3);
	DoubleVector coeff;
	int degree = ((x.size() + p.size()) / 3) - 1;
	bernsteinAll(degree - 2, u, coeff);
	size_t i = 0;
	autodiff::VectorXdual ctr1(3), ctr2(3), ctr3(3);
	size_t x_index = 0;
	size_t p_index = 0;
	autodiff::VectorXdual tmp(3 * (degree + 1));
	for (int i = 0; i < degree + 1; i++)
	{
		if (isStaticLocked(i))
		{
			tmp[3 * i] = p[3 * p_index];
			tmp[3 * i + 1] = p[3 * p_index + 1];
			tmp[3 * i + 2] = p[3 * p_index + 2];
			p_index++;
		}
		else
		{
			tmp[3 * i] = x[3 * x_index];
			tmp[3 * i + 1] = x[3 * x_index + 1];
			tmp[3 * i + 2] = x[3 * x_index + 2];
			x_index++;
		}

	}
	for (size_t i = 0; i < degree - 1; i++)
	{
		ctr1[0] = tmp[3 * i];
		ctr1[1] = tmp[3 * i + 1];
		ctr1[2] = tmp[3 * i + 2];
		ctr2[0] = tmp[3 * (i + 1)];
		ctr2[1] = tmp[3 * (i + 1) + 1];
		ctr2[2] = tmp[3 * (i + 1) + 2];
		ctr3[0] = tmp[3 * (i + 2)];
		ctr3[1] = tmp[3 * (i + 2) + 1];
		ctr3[2] = tmp[3 * (i + 2) + 2];
		der += degree * (degree - 1) * n_choose_k(degree - 2, i) * pow(u, (i)) * pow(1 - u, degree - (i)-2) *
			((ctr3 - ctr2) - (ctr2 - ctr1));
	}
	return der;
}
autodiff::dual BezierCurve::curvatureAutodiff(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u) {
	autodiff::dual curve;
	curve = vectorLength(crossProduct(firstDerivate(x, p, u),
		secondDerivate(x, p, u))) / pow(vectorLength(firstDerivate(x, p, u)), 3);
	return curve;
}
autodiff::dual BezierCurve::torsionAutodiff(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u) {
	autodiff::dual torsion;
	torsion = (dotProduct(crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u)), thirdDerivate(x, p, u))) /
		(pow(vectorLength(crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u))), 2));
	return torsion;
}
autodiff::dual BezierCurve::curveVariantAutodiff(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p, double u)
{
	autodiff::dual curveVariant;
	curveVariant = vectorLength(subtraction(scalarProduct((1 / pow(vectorLength(firstDerivate(x, p, u)), 3)),
		crossProduct(firstDerivate(x, p, u), thirdDerivate(x, p, u))),
		scalarProduct(3 * ((dotProduct(firstDerivate(x, p, u), secondDerivate(x, p, u))) /
		(pow(vectorLength(firstDerivate(x, p, u)), 4))),
			(crossProduct(firstDerivate(x, p, u), secondDerivate(x, p, u))))));
	return curveVariant;
}

autodiff::dual BezierCurve::curveVariantEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p)
{
	BezierCurve c(false);
	autodiff::dual curveVariant = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		curveVariant += numericalWeights[i] *
			((pow(c.curveVariantAutodiff(x, p, numericalFunctionValues[i]), 2)) / 
				c.vectorLength(c.firstDerivate(x, p, numericalFunctionValues[i])));
	}
	return curveVariant;
}
autodiff::dual BezierCurve::curveEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual sumCurve = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		sumCurve += numericalWeights[i] *
			(pow(c.curvatureAutodiff(x, p, numericalFunctionValues[i]), 2) *
				c.vectorLength(c.firstDerivate(x, p, numericalFunctionValues[i])));
	}
	return sumCurve;
}
autodiff::dual BezierCurve::torsionEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual torsion = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		torsion += numericalWeights[i] *
			(pow(c.torsionAutodiff(x, p, numericalFunctionValues[i]), 2)) 
			* c.vectorLength(c.firstDerivate(x,p, numericalFunctionValues[i]));
	}
	return torsion;
}
autodiff::dual BezierCurve::firstDerivateEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual energy = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		energy += numericalWeights[i] *
			(pow(c.vectorLength(c.firstDerivate(x, p, numericalFunctionValues[i])), 2));
	}
	return energy;
}
autodiff::dual BezierCurve::secondDerivateEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual energy = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		energy += numericalWeights[i] *
			(pow(c.vectorLength(c.secondDerivate(x, p, numericalFunctionValues[i])), 2));
	}
	return energy;
}
autodiff::dual BezierCurve::thirdDerivateEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual energy = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		energy += numericalWeights[i] *
			(pow(c.vectorLength(c.thirdDerivate(x, p, numericalFunctionValues[i])), 2));
	}
	return energy;
}


autodiff::dual BezierCurve::arcLengthEnergyFunction(const autodiff::VectorXdual & x,
	const autodiff::VectorXdual& p) {
	BezierCurve c(false);
	autodiff::dual arcLength = 0;
	for (size_t i = 0; i < sizeOfNumericalArrays; i++)
	{
		arcLength += numericalWeights[i]
			* c.vectorLength(c.firstDerivate(x, p, numericalFunctionValues[i]));
	}
	return arcLength;
}
#pragma endregion
VectorVector BezierCurve::getGradientVector()
{
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	autodiff::VectorXdual x((cp.size() * 3) - (lockedControlPointIndexes.size() * 3));
	autodiff::VectorXdual p(lockedControlPointIndexes.size() * 3);
	VectorVector vv;
	size_t index = 0;
	size_t p_index = 0;
	for (size_t i = 0; i < cp.size(); i++)
	{
		if (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), i)
			!= lockedControlPointIndexes.end())
		{
			p[3 * p_index] = cp[i].x;
			p[3 * p_index + 1] = cp[i].y;
			p[3 * p_index + 2] = cp[i].z;
			p_index++;
		}
		else
		{
			x[3 * index] = cp[i].x;
			x[3 * index + 1] = cp[i].y;
			x[3 * index + 2] = cp[i].z;
			index++;
		}
	}
	autodiff::dual u;
	Eigen::VectorXd g;
	auto fv = optimizationChoosing();

	g = autodiff::forward::gradient(fv,
		autodiff::wrt(x), autodiff::forward::at(x, p), u);
	if (g.size() != 0)
	{
		for (size_t i = 0; i < cp.size() - lockedControlPointIndexes.size(); i++)
		{
			Vector v = Vector(g[3 * i], g[3 * i + 1], g[3 * i + 2]);
			vv.push_back(v);
		}
	}

	return vv;
}
void BezierCurve::gradientDescend()
{
	VectorVector vv(getGradientVector());
	s_LockedControlPointIndexes = lockedControlPointIndexes;
	double previosEnergy = 0, currentEnergy = 100;
	if (vv.size() == 0)
	{
		return;
	}
	if (!lineSearch)
	{
		//iterationCount = 0;
		//while (abs(previosEnergy - currentEnergy) > pow(10,-5))
		//{
			autodiff::VectorXdual x1(3 * (getControlpoints().size() - lockedControlPointIndexes.size()));
			autodiff::VectorXdual p1(3 * lockedControlPointIndexes.size());
			auto optimizationFunctions = optimizationChoosing();
			controlPointsToAutodiffArguments(x1, p1);
			previosEnergy = optimizationFunctions(x1, p1).val;
			int index = 0;
			for (size_t i = 0; i < cp.size(); i++)
			{
				if (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), i)
					== lockedControlPointIndexes.end())
				{
					cp[i] -= vv[index++] * gradientDescendStepSize;
				}
			}
			controlPointsToAutodiffArguments(x1, p1);
			currentEnergy = optimizationFunctions(x1, p1).val;
			iterationCount++;
		//}

	}
	else
	{
		BezierCurve movedCurve;
		while (true)
		{
			int index = 0;
			movedCurve.n = n;
			movedCurve.cp = cp;
			movedCurve.lockedControlPointIndexes = this->lockedControlPointIndexes;
			for (size_t i = 0; i < movedCurve.cp.size(); i++)
			{
				if (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), i)
					== lockedControlPointIndexes.end())
				{
					movedCurve.cp[i] -= vv[index++] * gradientDescendStepSize;
				}
			}
			autodiff::VectorXdual x1(3 * (movedCurve.getControlpoints().size() - movedCurve.lockedControlPointIndexes.size()));
			autodiff::VectorXdual p1(3 * movedCurve.lockedControlPointIndexes.size());
			autodiff::VectorXdual x2(3 * (this->getControlpoints().size() - this->lockedControlPointIndexes.size()));
			autodiff::VectorXdual p2(3 * this->lockedControlPointIndexes.size());
			movedCurve.controlPointsToAutodiffArguments(x1, p1);
			this->controlPointsToAutodiffArguments(x2, p2);
			auto optimizationFunctions = optimizationChoosing();
			if (optimizationFunctions(x1, p1).val < optimizationFunctions(x2, p2).val)
			{
				index = 0;
				for (size_t i = 0; i < cp.size(); i++)
				{
					if (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), i)
						== lockedControlPointIndexes.end())
					{
						cp[i] -= vv[index++] * gradientDescendStepSize;
					}
				}
				gradientDescendStepSize *= 1.25;
				iterationCount++;
				break;
			}
			else
			{
				gradientDescendStepSize /= 2;
			}
		}

	}

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
	for (size_t i = 0; i < cp.size(); i++)
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

Point BezierCurve::crossProduct(Point p1, Point p2)
{
	Point result;
	result.x = p1.y * p2.z - p1.z * p2.y;
	result.y = p1.z * p2.x - p1.x * p2.z;
	result.z = p1.x * p2.y - p1.y * p2.x;
	return result;
}
autodiff::VectorXdual BezierCurve::crossProduct(const autodiff::VectorXdual& p1,
	const autodiff::VectorXdual& p2)
{
	autodiff::VectorXdual result(3);
	result[0] = p1[1] * p2[2] - p1[2] * p2[1];
	result[1] = p1[2] * p2[0] - p1[0] * p2[2];
	result[2] = p1[0] * p2[1] - p1[1] * p2[0];
	return result;
}

double BezierCurve::vectorLength(Vector v) const
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}
autodiff::dual BezierCurve::vectorLength(const autodiff::VectorXdual& v) const
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
autodiff::dual BezierCurve::dotProduct(const autodiff::VectorXdual& v1,
	const autodiff::VectorXdual& v2) const
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
autodiff::VectorXdual BezierCurve::scalarProduct(const autodiff::dual& scalar,
	const autodiff::VectorXdual& vec) const
{
	autodiff::VectorXdual res(3);
	for (size_t i = 0; i < 3; i++)
	{
		res[i] = scalar * vec[i];
	}
	return res;
}

autodiff::VectorXdual BezierCurve::subtraction(const autodiff::VectorXdual & vec1, const autodiff::VectorXdual & vec2) const
{
	autodiff::VectorXdual res(3);
	for (size_t i = 0; i < 3; i++)
	{
		res[i] = vec1[i] - vec2[i];
	}
	return res;
}

void BezierCurve::controlPointsToAutodiffArguments(autodiff::VectorXdual & x, autodiff::VectorXdual & p)
{
	int p_index = 0;
	int index = 0;
	for (size_t i = 0; i < cp.size(); i++)
	{
		if (std::find(lockedControlPointIndexes.begin(), lockedControlPointIndexes.end(), i)
			!= lockedControlPointIndexes.end())
		{
			p[3 * p_index] = cp[i].x;
			p[3 * p_index + 1] = cp[i].y;
			p[3 * p_index + 2] = cp[i].z;
			p_index++;
		}
		else
		{
			x[3 * index] = cp[i].x;
			x[3 * index + 1] = cp[i].y;
			x[3 * index + 2] = cp[i].z;
			index++;
		}
	}
}

autodiff::VectorXdual BezierCurve::division(const autodiff::VectorXdual & v, const autodiff::dual & u)
{
	autodiff::VectorXdual res(3);
	for (size_t i = 0; i < 3; i++)
	{
		res[i] = v[i] / u;
	}
	return res;
}