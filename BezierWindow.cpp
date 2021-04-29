#include "BezierWindow.h"

#include <QtWidgets>

BezierWindow::BezierWindow(QApplication * _parent) :parent(_parent), 
QMainWindow() {

	setWindowTitle(tr("Bezier Framework"));
	setStatusBar(new QStatusBar);

	auto quitAction = new QAction(tr("&Quit"), this);
	quitAction->setShortcut(tr("Ctrl+Q"));
	quitAction->setStatusTip(tr("Quit the program"));
	connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

	auto degreeAction = new QAction(tr("Set &degree"), this);
	degreeAction->setStatusTip(tr("Set the degree of the curve"));
	connect(degreeAction, SIGNAL(triggered()), this, SLOT(setDegreeInWindow()));

	auto fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(quitAction);

	auto visMenu = menuBar()->addMenu(tr("&Visualization"));
	visMenu->addAction(degreeAction);

	viewer = new BezierViewer(this);
	setCentralWidget(viewer);
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
		viewer->update();
	}
}

BezierWindow::~BezierWindow() {
}