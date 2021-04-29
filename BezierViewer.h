#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include "BezierCurve.h"
#include <string.h>

using qglviewer::Vec;

class BezierViewer : public QGLViewer{
	Q_OBJECT

public: 
	BezierViewer(QWidget* parent);
	virtual ~BezierViewer();
	size_t getDegree();
	void setDegree(size_t newDegree);

protected:
	virtual void init();
	virtual void draw();
	virtual void keyPressEvent(QKeyEvent * e) override;
	virtual void drawWithNames() override;
	virtual QString helpString() const override;
	virtual void postSelection(const QPoint &p) override;
	virtual void mouseMoveEvent(QMouseEvent *e) override;
	void drawControlPoints(BezierCurve c, Vector color) const;
	void drawAxes() const;	
	void setupDerivate();
private:
	void drawAxesWithNames() const;
	static Vec intersectLines(const Vec &ap, const Vec &ad,
		const Vec &bp, const Vec &bd);
	BezierCurve curve; 
	BezierCurve derivate;
	bool showControlPoints;	
	struct ModificationAxes {
		bool shown;
		float size;
		int selected_axis;
		Vec position, grabbed_pos, original_pos;
	} axes;
	int selected_vertex;
};

