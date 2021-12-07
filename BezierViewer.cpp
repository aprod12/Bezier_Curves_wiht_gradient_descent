#include "BezierViewer.h"
#include "BezierWindow.h"
#include <random>


BezierViewer::BezierViewer(QWidget * parent, BezierCurve* c) : showControlPoints(false),
curve(c), QGLViewer(parent), selected_vertex(-1), window(parent), showControlPoligon(false),
showCurvature(false), sumCurve(0), arcLengthByFractions(0), arcLenghtIntegral(0), m_Error(0),
showCurve(true), showFittingCurve(false)
{
	setSelectRegionWidth(10);
	setSelectRegionHeight(10);
}
BezierViewer::~BezierViewer()
{
	delete fittingCurve;
	delete curve;
}
void BezierViewer::setupDerivate() {
	derivate.setDegree(curve->getDegree() - 1);
	derivate.setControlPoints(PointVector(derivate.getDegree() + 1));
	curve->getFristDerivateControlPoints(derivate.getControlpoints());
}

void BezierViewer::draw() {
	glDisable(GL_LIGHTING);
	if (showControlPoints) {
		drawControlPoints(*curve, Vector(1.0, 0.0, 1.0));
	}
	if (showControlPoligon)
	{
		drawControlPoligon(Vector(0, 0.75, 1.0));
	}
	if (showCurve)
	{
		glLineWidth(3.0);
		glColor3d(0.96, 0.91, 0.28);
		glBegin(GL_LINE_STRIP);
		for (double i = 0.0; i <= 1.01; i += 0.01) {
			Point p = curve->evaluateByDeCasteljau(i);
			glVertex3f(p.x, p.y, p.z);
		}
		std::string tmp = "Iteration: " + std::to_string(curve->iterationCount) + 
		"| Fitting Curve Iteration: " + std::to_string(fittingCurve->iterationCount);
		displayMessage(QString::fromUtf8(tmp.c_str()));
		glEnd();
	}
	if (showCurvature)
	{
		for (double i = 0.0; i <= 1.01; i += 0.01)
		{
			glLineWidth(1.5);
			glColor3d(0.96, 0.1, 0.0);
			glBegin(GL_LINE_STRIP);
			Point p = curve->evaluateByDeCasteljau(i);
			autodiff::VectorXdual normal = curve->getNormalVector(i);
			double length = curve->getCurveSize(i).val;
			Point p2(p.x + normal[0].val * length,
				p.y + normal[1].val * length,
				p.z + normal[2].val * length);
			glVertex3f(p.x, p.y, p.z);
			glVertex3f(p2.x, p2.y, p2.z);
			glEnd();
		}
	}
	if (showFittingCurve)
	{
		glLineWidth(3.0);
		glColor3d(0.0, 1.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (double i = 0.0; i <= 1.01; i += 0.01)
		{
			Point p = fittingCurve->evaluateByDeCasteljau(i);
			glVertex3f(p.x, p.y, p.z);
		}
		glEnd();
		drawControlPoints(*fittingCurve, Vector(0, 0, 1));
		drawFittingPoints(*fittingCurve, Vector(0.99, 0.53, 0.07));

	}
	if (axes.shown)
		drawAxes();
}
void BezierViewer::keyPressEvent(QKeyEvent * e)
{
	if (e->modifiers() == Qt::NoModifier)
		switch (e->key()) {
		case Qt::Key_C:
			showControlPoints = !showControlPoints;
			update();
			break;
		case Qt::Key_R: {
			double curvatureSumary, fractionLength, archIntegral, error, torsion, firstDerivateEnergy,
				secondDerivateEnergy, thirdDerivateEnergy,
				curvatureVariant, fittingEnergy;
			curve->getInformations(curvatureSumary, fractionLength, archIntegral, error, torsion, firstDerivateEnergy,
				secondDerivateEnergy, thirdDerivateEnergy,
				curvatureVariant);
			fittingEnergy = curve->getFittingEnergy();
			static_cast<BezierWindow*>(window)->updateListOfDatas(
				curvatureSumary, fractionLength, archIntegral, error, torsion, firstDerivateEnergy,
				secondDerivateEnergy, thirdDerivateEnergy,
				curvatureVariant, fittingEnergy);

			sumCurve = curvatureSumary;
			arcLengthByFractions = fractionLength;
			arcLenghtIntegral = archIntegral;
			m_Error = error;
			update();
			break;
		}
		case Qt::Key_0:
			curve->gradientDescend();
			update();
			break;
		case Qt::Key_P:
			showControlPoligon = !showControlPoligon;
			update();
			break;
		case Qt::Key_T:
			showCurvature = !showCurvature;
			update();
		case Qt::Key_3:
			fittingCurve->gradientDescend();
			update();
			break;
		case Qt::Key_1:
			showCurve = !showCurve;
			update();
			break;
		case Qt::Key_2:
			showFittingCurve = !showFittingCurve;
			update();
			break;

		}
	QGLViewer::keyPressEvent(e);
}
void BezierViewer::drawControlPoints(BezierCurve c, Vector color) const {
	glLineWidth(1.0);
	glPointSize(8.0);
	glBegin(GL_POINTS);
	int index = 0;
	for (auto &p : c.getControlpoints())
	{
		if (c.isItLocked(index))
		{
			glColor3d(1.0, 0, 0);
		}
		else
		{
			glColor3d(color.x, color.y, color.z);
		}

		glVertex3f(p.x, p.y, p.z);
		index++;
	}
	glEnd();
}

void BezierViewer::drawFittingPoints(BezierCurve c, Vector color) const {
	glLineWidth(1.0);
	glPointSize(8.0);
	glBegin(GL_POINTS);
	int index = 0;
	for (auto &p : c.referencePoints)
	{
		if (c.isItLocked(index))
		{
			glColor3d(1.0, 0, 0);
		}
		else
		{
			glColor3d(color.x, color.y, color.z);
		}

		glVertex3f(p.x, p.y, p.z);
		index++;
	}
	glEnd();
}
size_t BezierViewer::getDegree()
{
	return curve->getDegree();
}
void BezierViewer::setDegree(size_t newDegree)
{
	curve->setDegree(newDegree);
}
void BezierViewer::changeCurve(BezierCurve * p_curve)
{
	delete curve;
	curve = p_curve;
	update();
}
QString BezierViewer::helpString() const
{
	QString text("<h2>Bezier Framework</h2>"
		"<p>This is a minimal framework for Bezirt curve manipulation.</p>"
		"<p>The following hotkeys are available:</p>"
		"<ul>"
		"<li>&nbsp;C: Toggle control point visualization</li>"
		"<li>&nbsp;R: Refresh curve information</li>"
		"<li>&nbsp;0: Gradient descent step</li>"
		"</ul>");
	return text;
}

void BezierViewer::init() {
	curve->regenerateControlPoints();
	fittingCurve = new BezierCurve((size_t)4);
	fittingCurve->makeFittingCurve();
	PointVector tmp = curve->fittingPoints(0.02);
	fittingCurve->referencePoints = tmp;
	curve->referencePoints = tmp;
}

void BezierViewer::drawAxes() const {
	const Vec &p = axes.position;
	glColor3d(1.0, 0.0, 0.0);
	drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
	glColor3d(0.0, 1.0, 0.0);
	drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
	glColor3d(0.0, 0.0, 1.0);
	drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
	glEnd();
}
void BezierViewer::drawWithNames() {
	if (axes.shown)
		return drawAxesWithNames();
	if (!showControlPoints)
		return;
	for (size_t i = 0; i < curve->getControlpoints().size(); ++i) {
		Vec const &p = Vec(curve->getControlpoints()[i].x,
			curve->getControlpoints()[i].y,
			curve->getControlpoints()[i].z);
		glPushName(i);
		glRasterPos3fv(p);
		glPopName();
	}
}

void BezierViewer::drawAxesWithNames() const {
	const Vec &p = axes.position;
	glPushName(0);
	drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
	glPopName();
	glPushName(1);
	drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
	glPopName();
	glPushName(2);
	drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
	glPopName();
}

void BezierViewer::postSelection(const QPoint &p) {
	int sel = selectedName();
	if (sel == -1) {
		axes.shown = false;
		return;
	}

	if (axes.shown) {
		axes.selected_axis = sel;
		bool found;
		axes.grabbed_pos = camera()->pointUnderPixel(p, found);
		axes.original_pos = axes.position;
		if (!found)
			axes.shown = false;
		return;
	}
	selected_vertex = sel;
	axes.position = Vec(curve->getControlpoints()[sel].x,
		curve->getControlpoints()[sel].y,
		curve->getControlpoints()[sel].z);
	double depth = camera()->projectedCoordinatesOf(axes.position)[2];
	Vec q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
	Vec q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
	axes.size = (q1 - q2).norm() / 10.0;
	axes.shown = true;
	axes.selected_axis = -1;
}

void BezierViewer::mouseMoveEvent(QMouseEvent *e) {
	if (!axes.shown || axes.selected_axis < 0 ||
		!(e->modifiers() & Qt::ShiftModifier) ||
		!(e->buttons() & Qt::LeftButton))
		return QGLViewer::mouseMoveEvent(e);

	Vec from, dir, axis(axes.selected_axis == 0,
		axes.selected_axis == 1,
		axes.selected_axis == 2);
	camera()->convertClickToLine(e->pos(), from, dir);
	auto p = intersectLines(axes.grabbed_pos, axis, from, dir);
	float d = (p - axes.grabbed_pos) * axis;
	Point cont(axes.position.x, axes.position.y, axes.position.z);
	axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
	curve->getControlpoints()[curve->getIndexOfContorlPoint(cont)] =
		Point(axes.position.x, axes.position.y, axes.position.z);
	update();
}
Vec BezierViewer::intersectLines(const Vec &ap, const Vec &ad,
	const Vec &bp, const Vec &bd) {
	double a = ad * ad, b = ad * bd, c = bd * bd;
	double d = ad * (ap - bp), e = bd * (ap - bp);
	if (a * c - b * b < 1.0e-7)
		return ap;
	double s = (b * e - c * d) / (a * c - b * b);
	return ap + s * ad;
}

void BezierViewer::drawControlPoligon(Vector color)
{
	glColor3d(color.x, color.y, color.z);
	glBegin(GL_LINE_STRIP);
	for (auto p : curve->getControlpoints()) {
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}
