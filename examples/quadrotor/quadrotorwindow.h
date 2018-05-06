#ifndef QUADROTORWINDOW_H
#define QUADROTORWINDOW_H

#include <QMainWindow>
#include <vector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    using Scalar = double;
    explicit MainWindow(Scalar dt, QWidget *parent = 0);
    ~MainWindow();

public slots:
    void update(const QVector<Scalar> &x, const QVector<Scalar> &u, const QVector<Scalar> &xf, Scalar true_cost);

private:
    Scalar dt_;
    Ui::MainWindow *ui;
};

#endif // QUADROTORWINDOW_H
