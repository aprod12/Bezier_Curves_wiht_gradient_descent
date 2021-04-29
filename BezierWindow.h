#pragma once
#include <QtWidgets/QMainWindow>

#include "BezierViewer.h"

class BezierWindow : public QMainWindow {
	Q_OBJECT

public:
	explicit BezierWindow(QApplication *parent);
	~BezierWindow();
private slots:
	void setDegreeInWindow();

private:
	QApplication* parent;
	BezierViewer* viewer;
	
};
