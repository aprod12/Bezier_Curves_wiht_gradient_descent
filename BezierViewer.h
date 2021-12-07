#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include "IDrawable.h"
#include "BezierCurve.h"
#include <string.h>

using qglviewer::Vec;

class BezierViewer : public QGLViewer{
	Q_OBJECT

public: 
	BezierViewer(QWidget* parent, BezierCurve* c);
	virtual ~BezierViewer();
	size_t getDegree();
	void setDegree(size_t newDegree);
	void changeCurve(BezierCurve* p_curve);

protected:
	virtual void init();
	virtual void draw();
	virtual void keyPressEvent(QKeyEvent * e) override;
	virtual void drawWithNames() override;
	virtual QString helpString() const override;
	virtual void postSelection(const QPoint &p) override;
	virtual void mouseMoveEvent(QMouseEvent *e) override;
	void drawControlPoints(BezierCurve c, Vector color) const;
	void drawFittingPoints(BezierCurve c, Vector color) const;
	void drawAxes() const;	
	void setupDerivate();
private:
	void drawAxesWithNames() const;
	static Vec intersectLines(const Vec &ap, const Vec &ad,
		const Vec &bp, const Vec &bd);
	//TODO IDrawable legyen
	BezierCurve* curve;
	BezierCurve* fittingCurve;
	BezierCurve derivate;
	QWidget* window;
	bool showControlPoints;	
	bool showControlPoligon;
	bool showCurvature;
	bool showCurve;
	bool showFittingCurve;
	struct ModificationAxes {
		bool shown;
		float size;
		int selected_axis;
		Vec position, grabbed_pos, original_pos;
	} axes;
	int selected_vertex;
	double sumCurve, arcLengthByFractions, arcLenghtIntegral, m_Error;
	void drawControlPoligon(Vector color);
};

