#pragma once



#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QtWidgets>
#include <QWidget>
#include "BezierCurve.h"

class DockWidget : public QWidget {
	Q_OBJECT
private:

	QLabel* controlPointLabel;
	QPushButton* lockButton;
	QComboBox* listOfControlPoints;
	QPushButton* quadratureDegreeButton;
	QComboBox* quadratureDegreeList;
	QLabel* degreeLabel;
	QLabel* quadratureDegreeLabel;
	QPlainTextEdit* degreeCount;
	QPushButton* increaseButton;
	QComboBox* listOfEnergies;
	QPushButton* energyConfirm;
	QListWidget* listOfDatas;
	QPushButton* gradientDescendButton;
	QDoubleSpinBox* gradientDescendStepSize;
	QCheckBox* lineSearchCheckBox;
	QLabel* lineSearchLabel;
	BezierCurve* curve;
	QWidget* parentWidget;

private slots:
	void setDegree();
	void setOptimization();
	void setLineSearch();
	void setStepSize();
	void gradientDescendStep();
	void setQuadratureDegree();
	void changeToCorrectButtonText();
	void lock();

public:
	void changeCurve(BezierCurve* p_curve);
	DockWidget(BezierCurve* c, QWidget *parent = nullptr);
	void updateListOfDatas(double curvatureSumary, double fractionLength,
		double archIntegral, double error, double torsion, double firstDerivateEnergy,
		double secondDerivateEnergy, double thirdDerivateEnergy,
		double curvatureVariant, double fittingEnergy);
	
};