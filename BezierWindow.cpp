#include "BezierWindow.h"

#include <QtWidgets>

BezierWindow::BezierWindow(QApplication * _parent) :parent(_parent),
QMainWindow() {


	BezierCurve* c = new BezierCurve((size_t)4);
	setWindowTitle(tr("Bezier Framework"));
	setStatusBar(new QStatusBar);

	auto quitAction = new QAction(tr("&Quit"), this);
	quitAction->setShortcut(tr("Ctrl+Q"));
	quitAction->setStatusTip(tr("Quit the program"));
	connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

	auto degreeAction = new QAction(tr("Set &degree"), this);
	degreeAction->setStatusTip(tr("Set the degree of the curve"));
	connect(degreeAction, SIGNAL(triggered()), this, SLOT(setDegreeInWindow()));

	auto resetAction = new QAction(tr("&Reset"), this);
	quitAction->setShortcut(tr("Ctrl+R"));
	quitAction->setStatusTip(tr("Resets the curve"));
	connect(resetAction, SIGNAL(triggered()), this, SLOT(reset()));

	auto fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(quitAction);
	fileMenu->addAction(resetAction);

	auto visMenu = menuBar()->addMenu(tr("&Visualization"));
	visMenu->addAction(degreeAction);


	viewer = new BezierViewer(this, c);
	setCentralWidget(viewer);
	createDockWindows(c);
}
void BezierWindow::setDegreeInWindow() {

	auto dlg = std::make_unique<QDialog>(this);
	auto *hb1 = new QHBoxLayout,
		*hb2 = new QHBoxLayout;
	auto *vb = new QVBoxLayout;
	auto *text = new QLabel(tr("Degree:"));
	auto *sb = new QSpinBox;
	auto *cancel = new QPushButton(tr("Cancel"));
	auto *ok = new QPushButton(tr("Ok"));

	sb->setSingleStep(1);
	sb->setValue(viewer->getDegree());
	connect(cancel, SIGNAL(pressed()), dlg.get(), SLOT(reject()));
	connect(ok, SIGNAL(pressed()), dlg.get(), SLOT(accept()));
	ok->setDefault(true);

	hb1->addWidget(text);
	hb1->addWidget(sb);
	hb2->addWidget(cancel);
	hb2->addWidget(ok);
	vb->addLayout(hb1);
	vb->addLayout(hb2);

	dlg->setWindowTitle(tr("Set degree"));
	dlg->setLayout(vb);

	if (dlg->exec() == QDialog::Accepted) {
		viewer->setDegree(sb->value());

	}
}

void BezierWindow::createDockWindows(BezierCurve* c)
{
	QDockWidget *dock = new QDockWidget(tr("Actions"), this);
	dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	widget = new DockWidget(c, this);
	dock->setWidget(widget);
	addDockWidget(Qt::RightDockWidgetArea, dock);
}

void BezierWindow::updateListOfDatas(double curvatureSumary, double fractionLength,
	double archIntegral, double error, double torsion, double firstDerivateEnergy,
	double secondDerivateEnergy, double thirdDerivateEnergy,
	double curvatureVariant, double fittingEnergy)
{
	widget->updateListOfDatas(curvatureSumary, fractionLength, archIntegral, error, torsion, firstDerivateEnergy,
		secondDerivateEnergy, thirdDerivateEnergy,
		curvatureVariant, fittingEnergy);
}

BezierWindow::~BezierWindow() {
}

void BezierWindow::updateViewer()
{
	viewer->update();
}

void BezierWindow::reset()
{
	BezierCurve* c = new BezierCurve((size_t)4);
	c->regenerateControlPoints();
	widget->changeCurve(c);
	viewer->changeCurve(c);
}
