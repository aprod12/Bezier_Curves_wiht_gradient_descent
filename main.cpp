#include <QtWidgets/QApplication>

#include "BezierWindow.h"
#include "MyWindow.h"

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  BezierWindow window(&app);
  int x = window.width()* 2.0;
  int y = window.height()* 2.0;
  window.setFixedSize(x, y);
  window.show();
  return app.exec();
}
