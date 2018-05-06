#ifndef AUTORALLYRTWINDOW_H
#define AUTORALLYRTWINDOW_H

#include <QMainWindow>
#include <qcustomplot/qcustomplot.h>

namespace Ui {
class AutoRallyRTWindow;
}

class AutoRallyRTWindow : public QMainWindow
{
    Q_OBJECT

public:
    using Scalar = double;
    explicit AutoRallyRTWindow(Scalar dt, QWidget *parent = 0);
    ~AutoRallyRTWindow();

public slots:
    void update_graph(const QVector<Scalar> &state, const QVector<Scalar> &control);

private:
    Scalar dt_;
    Ui::AutoRallyRTWindow *ui;
    QCPCurve *position;
};

#endif // AUTORALLYRTWINDOW_H
