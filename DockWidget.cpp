#include "DockWidget.h"
#include "BezierWindow.h"

void DockWidget::changeCurve(BezierCurve * p_curve)
{
	curve = p_curve;
	listOfControlPoints->clear();
	degreeCount->setPlainText(QString::number(curve->getDegree()));
	for (size_t i = 0; i < curve->getDegree() + 1; i++)
	{
		listOfControlPoints->addItem(QString::number(i));
	}
}

DockWidget::DockWidget(BezierCurve* c, QWidget *parent)
	: degreeLabel(nullptr), degreeCount(nullptr), increaseButton(nullptr), parentWidget(parent),
	  listOfEnergies(nullptr), energyConfirm(nullptr), listOfDatas(nullptr), gradientDescendButton(nullptr),
	  gradientDescendStepSize(nullptr), lineSearchCheckBox(nullptr), lineSearchLabel(nullptr), curve(c), QWidget(parent)
{


	auto* hb1 = new QHBoxLayout;
	auto* hbBottom = new QHBoxLayout;
	auto* hb2 = new QHBoxLayout;
	auto* hb3 = new QHBoxLayout;
	auto* hb4 = new QHBoxLayout;
	auto* hb5 = new QHBoxLayout;

	degreeLabel = new QLabel("Deegre:", this);
	degreeCount = new QPlainTextEdit(this);
	increaseButton = new QPushButton(tr("Increese"));

	listOfEnergies = new QComboBox(this);
	energyConfirm = new QPushButton(tr("Set"));
	listOfDatas = new QListWidget(this);
	gradientDescendButton = new QPushButton(tr("Step"));
	gradientDescendStepSize = new QDoubleSpinBox();
	lineSearchCheckBox = new QCheckBox();
	lineSearchLabel = new QLabel("Line Search");
	quadratureDegreeButton= new QPushButton("Set");
	quadratureDegreeList = new QComboBox();
	quadratureDegreeLabel = new QLabel("Quadrature Degree:");
	listOfControlPoints = new QComboBox(this);
	lockButton = new QPushButton(tr("Unlock"));
	controlPointLabel = new QLabel(tr("ControlPoint Index:"));



	degreeCount->setReadOnly(true);
	degreeCount->setFixedHeight(25);
	degreeCount->setPlainText(QString::number(curve->getDegree()));
	connect(increaseButton, SIGNAL(clicked()), this, SLOT(setDegree()));

	hb1->addWidget(degreeLabel);
	hb1->addWidget(degreeCount);
	hb1->addWidget(increaseButton);

	connect(energyConfirm, SIGNAL(clicked()), this, SLOT(setOptimization()));
	listOfEnergies->addItem("Curvature");
	listOfEnergies->addItem("ArcLength");
	listOfEnergies->addItem("Torsion");
	listOfEnergies->addItem("First derivate square");
	listOfEnergies->addItem("Second derivate square");
	listOfEnergies->addItem("Third derivate square");
	listOfEnergies->addItem("Curve variant");
	hb2->addWidget(listOfEnergies);
	hb2->addWidget(energyConfirm);

	gradientDescendStepSize->setMinimumWidth(200);
	gradientDescendStepSize->setSingleStep(0.0001);
	gradientDescendStepSize->setDecimals(5);
	gradientDescendStepSize->setValue(curve->getGradientDescendStepSize());

	connect(gradientDescendButton, SIGNAL(clicked()), this, SLOT(gradientDescendStep()));
	hb3->addWidget(gradientDescendButton);

	connect(gradientDescendStepSize, SIGNAL(editingFinished()), this, SLOT(setStepSize()));
	hb3->addWidget(gradientDescendStepSize);
	hb3->addWidget(lineSearchLabel);

	connect(lineSearchCheckBox, SIGNAL(clicked()), this, SLOT(setLineSearch()));
	hb3->addWidget(lineSearchCheckBox);

	quadratureDegreeList->addItem(QString::number(20));
	quadratureDegreeList->addItem(QString::number(128));
	quadratureDegreeList->addItem(QString::number(256));
	connect(quadratureDegreeButton, SIGNAL(clicked()), this, SLOT(setQuadratureDegree()));
	hb4->addWidget(quadratureDegreeLabel);
	hb4->addWidget(quadratureDegreeList);
	hb4->addWidget(quadratureDegreeButton);

	for (size_t i = 0; i < curve->getDegree() + 1; i++)
	{
		listOfControlPoints->addItem(QString::number(i));
	}

	connect(listOfControlPoints, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [=](int index) { changeToCorrectButtonText(); });
	connect(lockButton, SIGNAL(pressed()), this, SLOT(lock()));
	hb5->addWidget(controlPointLabel);
	hb5->addWidget(listOfControlPoints);
	hb5->addWidget(lockButton);
	hbBottom->addWidget(listOfDatas);


	auto *vb = new QVBoxLayout;

	vb->addLayout(hb1);
	vb->addLayout(hb2);
	vb->addLayout(hb3);
	vb->addLayout(hb4);
	vb->addLayout(hb5);
	vb->addLayout(hbBottom);

	setLayout(vb);
}
void DockWidget::updateListOfDatas(double curvatureSumary, double fractionLength,
	double archIntegral, double error, double torsion, double firstDerivateEnergy,
	double secondDerivateEnergy, double thirdDerivateEnergy,
	double curvatureVariant, double fittingEnergy)
{
	listOfDatas->clear();

	QString sumCurveString = "Curvature summary: ";
	sumCurveString.append(QString::number(curvatureSumary));
	QString arcLengthByFractionsString = "Arc length by fractions: ";
	arcLengthByFractionsString.append(QString::number(fractionLength));
	QString arcLenghtIntegralString = "Arc length integral: ";
	arcLenghtIntegralString.append(QString::number(archIntegral));
	QString errorString = "Error between the arc lengths: ";
	errorString.append(QString::number(error));
	QString torsionString = "Torsion: ";
	torsionString.append(QString::number(torsion));
	QString firstDerivateString = "Energy calculated by first derivate: ";
	firstDerivateString.append(QString::number(firstDerivateEnergy));
	QString secondDerivateString = "Energy calculated by second derivate: ";
	secondDerivateString.append(QString::number(secondDerivateEnergy));
	QString thirdDerivateString = "Energy calculated by third derivate: ";
	thirdDerivateString.append(QString::number(thirdDerivateEnergy));
	QString curvatureVariantString = "Curvature variance: ";
	curvatureVariantString.append(QString::number(curvatureVariant));
	QString fittingString = "The fitting curve's error compared to the reference points: ";
	fittingString.append(QString::number(fittingEnergy));
	
	listOfDatas->addItem(sumCurveString);
	listOfDatas->addItem(arcLengthByFractionsString);
	listOfDatas->addItem(arcLenghtIntegralString);
	listOfDatas->addItem(errorString);
	listOfDatas->addItem(torsionString);
	listOfDatas->addItem(firstDerivateString);
	listOfDatas->addItem(secondDerivateString);
	listOfDatas->addItem(thirdDerivateString);
	listOfDatas->addItem(curvatureVariantString);
	listOfDatas->addItem(fittingString);

}
void DockWidget::setOptimization()
{
	curve->setMode(listOfEnergies->currentIndex());
}

void DockWidget::setDegree()
{
	curve->increaseDegree();
	degreeCount->setPlainText(QString::number(curve->getDegree()));
	listOfControlPoints->addItem(QString::number(curve->getDegree() + 1));
	static_cast<BezierWindow*>(parentWidget)->updateViewer();

}

void DockWidget::setLineSearch()
{
	curve->setLineSearch(lineSearchCheckBox->isChecked());
}

void DockWidget::setStepSize() 
{
	curve->setGradientDescendStepSize(gradientDescendStepSize->value());
}

void DockWidget::gradientDescendStep()
{
	curve->gradientDescend();
	static_cast<BezierWindow*>(parentWidget)->updateViewer();
}

void DockWidget::setQuadratureDegree()
{
	QString listElement = quadratureDegreeList->itemText(quadratureDegreeList->currentIndex());
	curve->setQuadratureDegree(listElement.toInt());
}

void DockWidget::changeToCorrectButtonText()
{
	if (curve->isItLocked(listOfControlPoints->currentIndex()))
	{
		lockButton->setText("Unlock");
	}
	else
	{
		lockButton->setText("Lock");
	}
}

void DockWidget::lock()
{
	if (curve->isItLocked(listOfControlPoints->currentIndex()))
	{
		curve->unlockControlPoint(listOfControlPoints->currentIndex());
		lockButton->setText("Lock");
	}
	else
	{
		curve->lockControlPoint(listOfControlPoints->currentIndex());
		lockButton->setText("Unlock");
	}
	static_cast<BezierWindow*>(parentWidget)->updateViewer();
}
