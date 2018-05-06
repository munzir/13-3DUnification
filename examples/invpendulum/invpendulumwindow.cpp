#include "invpendulumwindow.h"
#include "ui_invpendulumwindow.h"

InvPendulumWindow::InvPendulumWindow(const std::vector<double> &angle, const std::vector<double> &anglef,
                                     const std::vector<double> &velocity, const std::vector<double> &velocityf,
                                     const std::vector<double> &time, const std::vector<double> &cost, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::InvPendulumWindow)
{
    ui->setupUi(this);

    // Angle
    ui->angle_plot->addGraph();
    ui->angle_plot->addGraph();
    ui->angle_plot->legend->setVisible(true);
    ui->angle_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->angle_plot->graph(0)->setName("Target");
    ui->angle_plot->graph(1)->setName("Trajectory");
    ui->angle_plot->xAxis->setLabel("Time (s)");
    ui->angle_plot->yAxis->setLabel("Angle (rad)");
    ui->angle_plot->graph(0)->setPen(QPen(Qt::red));
    ui->angle_plot->graph(0)->setData(QVector<double>::fromStdVector(time), QVector<double>::fromStdVector(anglef));
    ui->angle_plot->graph(1)->setData(QVector<double>::fromStdVector(time), QVector<double>::fromStdVector(angle));
    ui->angle_plot->rescaleAxes(true);
    ui->angle_plot->yAxis->setRange(ui->angle_plot->yAxis->range().upper + 1, ui->angle_plot->yAxis->range().lower - 1);

    // Position
    ui->velocity_plot->addGraph();
    ui->velocity_plot->addGraph();
    ui->velocity_plot->legend->setVisible(true);
    ui->velocity_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->velocity_plot->graph(0)->setName("Target");
    ui->velocity_plot->graph(1)->setName("Trajectory");
    ui->velocity_plot->xAxis->setLabel("Time (s)");
    ui->velocity_plot->yAxis->setLabel("Angular Velocity (rad/s)");
    ui->velocity_plot->graph(0)->setPen(QPen(Qt::red));
    ui->velocity_plot->graph(0)->setData(QVector<double>::fromStdVector(time), QVector<double>::fromStdVector(velocityf));
    ui->velocity_plot->graph(1)->setData(QVector<double>::fromStdVector(time), QVector<double>::fromStdVector(velocity));
    ui->velocity_plot->rescaleAxes(true);
    ui->velocity_plot->yAxis->setRange(ui->velocity_plot->yAxis->range().upper + 1, ui->velocity_plot->yAxis->range().lower - 1);

    // Cost
    ui->cost_plot->addGraph();
    ui->cost_plot->xAxis->setLabel("Iteration #");
    ui->cost_plot->yAxis->setLabel("Cost");
    std::vector<double> iter(cost.size(), 1);
    double n = 1.0;
    std::generate(iter.begin(), iter.end(), [&n]{ return n++; });
    ui->cost_plot->graph(0)->setData(QVector<double>::fromStdVector(iter), QVector<double>::fromStdVector(cost));
    ui->cost_plot->rescaleAxes(true);
    ui->cost_plot->xAxis->setAutoTickStep(false);
    ui->cost_plot->xAxis->setTickStep(2.0);
}

InvPendulumWindow::~InvPendulumWindow()
{
    delete ui;
}
