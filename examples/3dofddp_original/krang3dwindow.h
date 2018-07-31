#ifndef KRANG3DWINDOW_H
#define KRANG3DWINDOW_H

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
    explicit MainWindow(Scalar pos_xf, Scalar pos_yf, Scalar anglef, Scalar dt, QWidget *parent = 0);
    ~MainWindow();

public slots:
    void update(const QVector<Scalar> &x, const QVector<Scalar> &u, const QVector<Scalar> &x_ref,
                const QVector<Scalar> &u_ref);

private:
    Scalar anglef_;
    Scalar pos_xf_;
    Scalar pos_yf_;
    Scalar dt_;
    Ui::MainWindow *ui;
};

#endif // KRANG3DWINDOW_H
