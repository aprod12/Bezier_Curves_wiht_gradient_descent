#include "BezierViewer.h"
#include <random>


BezierViewer::BezierViewer(QWidget * parent): showControlPoints(false), 
QGLViewer(parent), selected_vertex(-1),
sumCurve(0), arcLengthByFractions(0), arcLenghtIntegral(0), error(0)
{
	setSelectRegionWidth(10);
	setSelectRegionHeight(10);
}
BezierViewer::~BezierViewer()
{
}
void BezierViewer::setupDerivate() {
	derivate.n = curve.n - 1;
	derivate.cp = PointVector(derivate.n + 1);
	curve.getFristDerivateControlPoints(derivate.cp);
}

void BezierViewer::draw() {
	glDisable(GL_LIGHTING);
	if (showControlPoints) {
		drawControlPoints(curve, Vector(1.0, 0.0, 1.0));
	}
	glLineWidth(3.0);
	glColor3d(0.96, 0.91, 0.28);
	glBegin(GL_LINE_STRIP);
	curve.regenerateControlPoints();
	for (double i = 0.0; i <= 1.01; i += 0.01) {
		Point p = curve.evaluateByDeCasteljau(i);
		glVertex3f(p.x, p.y, p.z);
	}
	auto a = curve.getGradientVector();
	double ossz = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		ossz += pow(a[i].norm(), 2);
	}
	ossz = sqrt(ossz);
	std::string tmp = "Arc length by fractions: " + std::to_string(arcLengthByFractions) + " | " +
		"Arc length by integral: "  + std::to_string(arcLenghtIntegral) + " | Error: " + 
		std::to_string(error) + " | Curve sum: " + std::to_string(sumCurve) 
		+" | Gradient norm: " + std::to_string(ossz);
	displayMessage(QString::fromUtf8(tmp.c_str()));
	glEnd();
	glLineWidth(3.0);
	glColor3d(0.0, 0.91, 0.28);
	glBegin(GL_LINE_STRIP);
	for (double i = 0.0; i <= 1.01; i += 0.01) {
		Point p = curve.getFirstDerivatedValue(i);
		//glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
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
				autodiff::VectorXdual x((curve.cp.size() * 3) - 6);
				autodiff::VectorXdual p(6);
				VectorVector vv;
				size_t index = 0;
				for (size_t i = 1; i < curve.cp.size() - 1; i++)
				{
					x[3 * index] = curve.cp[i].x;
					x[3 * index + 1] = curve.cp[i].y;
					x[3 * index + 2] = curve.cp[i].z;
					index++;
				}
				p[0] = curve.cp[0].x;
				p[1] = curve.cp[0].y;
				p[2] = curve.cp[0].z;
				p[3] = curve.cp[curve.cp.size() - 1].x;
				p[4] = curve.cp[curve.cp.size() - 1].y;
				p[5] = curve.cp[curve.cp.size() - 1].z;
				sumCurve = curve.energyFunction(x,p).val;
				arcLengthByFractions = curve.arcLengthByFractions();
				arcLenghtIntegral = curve.arcLengthByNumericalIntegral();
				error = abs(arcLengthByFractions - arcLenghtIntegral);
				update();
				break;
			}
			case Qt::Key_0:
				curve.gradientDescend();
				update();
				break;

		}
		QGLViewer::keyPressEvent(e);
}
void BezierViewer::drawControlPoints(BezierCurve c, Vector color) const {
	glLineWidth(1.0);
	glPointSize(8.0);
	glColor3d(color.x, color.y, color.z);
	glBegin(GL_POINTS);
	for (auto &p : c.cp)
		glVertex3f(p.x, p.y, p.z);
	glEnd();
}
size_t BezierViewer::getDegree()
{
	return curve.n;
}
void BezierViewer::setDegree(size_t newDegree)
{
	curve.n = newDegree;
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
	curve.n = 3;
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
	for (size_t i = 0; i < curve.cp.size(); ++i) {
			Vec const &p = Vec(curve.cp[i].x, curve.cp[i].y, curve.cp[i].z);
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
	axes.position = Vec(curve.cp[sel].x, curve.cp[sel].y, curve.cp[sel].z);
	double depth = camera()->projectedCoordinatesOf(axes.position)[2];
	Vec q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
	Vec q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
	axes.size = (q1 - q2).norm() / 10.0;
	axes.shown = true;
	axes.selected_axis = -1;
}

void BezierViewer::mouseMoveEvent(QMouseEvent *e) {
	if (!axes.shown || axes.selected_axis < 0 ||
		!(e->modifiers() & Qt::ShiftModifier)||
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
	curve.cp[curve.getIndexOfContorlPoint(cont)] =
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