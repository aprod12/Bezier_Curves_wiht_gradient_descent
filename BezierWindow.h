#pragma once
#include <QtWidgets/QMainWindow>
#include "DockWidget.h"
#include "BezierViewer.h"

class BezierWindow : public QMainWindow {
	Q_OBJECT

public:
	explicit BezierWindow(QApplication *parent);
	~BezierWindow();
	void createDockWindows(BezierCurve* c);
	void updateListOfDatas(double curvatureSumary, double fractionLength,
		double archIntegral, double error, double torsion, double firstDerivateEnergy,
		double secondDerivateEnergy, double thirdDerivateEnergy,
		double curvatureVariant, double fittingEnergy);
	void updateViewer();
private slots:
	void setDegreeInWindow();
	void reset();
private:
	QApplication* parent;
	BezierViewer* viewer;
	DockWidget* widget;

	
};
