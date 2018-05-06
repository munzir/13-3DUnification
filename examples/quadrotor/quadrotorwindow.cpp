#include "quadrotorwindow.h"
#include "ui_quadrotorwindow.h"

MainWindow::MainWindow(double dt, QWidget *parent) :
    QMainWindow(parent),
    dt_(dt),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // x Position
    ui->x_plot->addGraph();
    ui->x_plot->addGraph();
    ui->x_plot->legend->setVisible(true);
    ui->x_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->x_plot->graph(0)->setName("Target");
    ui->x_plot->graph(1)->setName("Trajectory");
    ui->x_plot->xAxis->setLabel("Time (s)");
    ui->x_plot->yAxis->setLabel("x Position (m)");
    ui->x_plot->graph(0)->setPen(QPen(Qt::red));
    ui->x_plot->rescaleAxes(true);

    // y Position
    ui->y_plot->addGraph();
    ui->y_plot->addGraph();
    ui->y_plot->legend->setVisible(true);
    ui->y_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->y_plot->graph(0)->setName("Target");
    ui->y_plot->graph(1)->setName("Trajectory");
    ui->y_plot->xAxis->setLabel("Time (s)");
    ui->y_plot->yAxis->setLabel("y Position (m)");
    ui->y_plot->graph(0)->setPen(QPen(Qt::red));
    ui->y_plot->rescaleAxes(true);

    // z Position
    ui->z_plot->addGraph();
    ui->z_plot->addGraph();
    ui->z_plot->legend->setVisible(true);
    ui->z_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->z_plot->graph(0)->setName("Target");
    ui->z_plot->graph(1)->setName("Trajectory");
    ui->z_plot->xAxis->setLabel("Time (s)");
    ui->z_plot->yAxis->setLabel("z Position (m)");
    ui->z_plot->graph(0)->setPen(QPen(Qt::red));
    ui->z_plot->rescaleAxes(true);

    // Cost
    ui->cost_plot->addGraph();
    ui->cost_plot->xAxis->setLabel("Time Step");
    ui->cost_plot->yAxis->setLabel("True Cost");
    ui->cost_plot->rescaleAxes(true);
}

void MainWindow::update(const QVector<Scalar> &x, const QVector<Scalar> &u, const QVector<Scalar> &xf, Scalar cost)
{
    static Scalar time = 0.0;
    time += dt_;
    ui->x_plot->graph(0)->addData(time, xf[0]);
    ui->x_plot->graph(1)->addData(time, x[0]);
    ui->y_plot->graph(0)->addData(time, xf[1]);
    ui->y_plot->graph(1)->addData(time, x[1]);
    ui->z_plot->graph(0)->addData(time, xf[2]);
    ui->z_plot->graph(1)->addData(time, x[2]);
    ui->cost_plot->graph(0)->addData(time, cost);
    ui->x_plot->rescaleAxes(true);
    ui->x_plot->replot();
    ui->y_plot->rescaleAxes(true);
    ui->y_plot->replot();
    ui->z_plot->rescaleAxes(true);
    ui->z_plot->replot();
    ui->cost_plot->rescaleAxes(true);
    ui->cost_plot->replot();
}

MainWindow::~MainWindow()
{
    delete ui;
}
