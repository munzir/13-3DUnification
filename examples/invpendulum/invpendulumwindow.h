#ifndef INVPENDULUMWINDOW_H
#define INVPENDULUMWINDOW_H

#include <QMainWindow>
#include <vector>

namespace Ui {
class InvPendulumWindow;
}

class InvPendulumWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit InvPendulumWindow(const std::vector<double> &angle, const std::vector<double> &anglef,
                               const std::vector<double> &velocity, const std::vector<double> &velocityf,
                               const std::vector<double> &time, const std::vector<double> &cost, QWidget *parent = 0);
    ~InvPendulumWindow();

private:
    Ui::InvPendulumWindow *ui;
};

#endif // INVPENDULUMWINDOW_H
