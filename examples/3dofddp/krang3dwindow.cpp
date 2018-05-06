#include "krang3dwindow.h"
#include "ui_krang3dwindow.h"

MainWindow::MainWindow(Scalar pos_xf, Scalar pos_yf, Scalar anglef, Scalar dt, QWidget *parent) :
    QMainWindow(parent),
    anglef_(anglef), pos_xf_(pos_xf), pos_yf_(pos_yf), dt_(dt),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Angle
    ui->angle_plot->addGraph();
    ui->angle_plot->addGraph();
    ui->angle_plot->addGraph();
    ui->angle_plot->legend->setVisible(true);
    ui->angle_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->angle_plot->graph(0)->setName("Target");
    ui->angle_plot->graph(1)->setName("Trajectory");
    ui->angle_plot->graph(2)->setName("Reference Trajectory");
    ui->angle_plot->xAxis->setLabel("Time (s)");
    ui->angle_plot->yAxis->setLabel("Angle (rad)");
    ui->angle_plot->graph(0)->setPen(QPen(Qt::red));
    ui->angle_plot->graph(2)->setPen(QPen(Qt::green));
    ui->angle_plot->rescaleAxes(true);
    ui->angle_plot->yAxis->setRange(-2.0 * M_PI, 2.0 * M_PI);

    // Position
    ui->position_plot->addGraph();
    ui->position_plot->addGraph();
    ui->position_plot->addGraph();
    ui->position_plot->legend->setVisible(true);
    ui->position_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignLeft);
    ui->position_plot->graph(0)->setName("Target");
    ui->position_plot->graph(1)->setName("Trajectory");
    ui->position_plot->graph(2)->setName("Reference Trajectory");
    ui->position_plot->xAxis->setLabel("Position X");
    ui->position_plot->yAxis->setLabel("Position Y");
    QPen targetPen;
    targetPen.setColor(Qt::red);
    targetPen.setWidth(5);
    QPen trajPen;
    trajPen.setColor(Qt::blue);
    trajPen.setWidth(5);
    QPen refPen;
    refPen.setColor(Qt::green);
    refPen.setWidth(3);
    ui->position_plot->graph(0)->setPen(targetPen);
    ui->position_plot->graph(1)->setPen(trajPen);
    ui->position_plot->graph(2)->setPen(refPen);
//    ui->position_plot->rescaleAxes(true);
    ui->position_plot->yAxis->setRange(-5.0, 5.0);
    ui->position_plot->xAxis->setRange(-5.0, 5.0);

    // Control
    ui->control_plot->addGraph();
    ui->control_plot->addGraph();
    ui->control_plot->graph(0)->setName("Control");
    ui->control_plot->graph(1)->setName("Reference Control");
    ui->control_plot->graph(1)->setPen(QPen(Qt::green));
    ui->control_plot->xAxis->setLabel("Time (s)");
    ui->control_plot->yAxis->setLabel("Control");
    ui->control_plot->rescaleAxes(true);
}

void MainWindow::update(const QVector<Scalar> &x, const QVector<Scalar> &u, const QVector<Scalar> &x_ref, const QVector<Scalar> &u_ref)
{
    static Scalar time = 0.0;
    time += dt_;
    double last_traj_x;
    double last_ref_x;

    ui->position_plot->graph(0)->addData(pos_xf_, pos_yf_);
//    ui->position_plot->graph(1)->data()->remove(last_traj_x);
    ui->position_plot->graph(1)->addData(x[6], x[7]);
//    ui->position_plot->graph(2)->data()->remove(last_ref_x);
    ui->position_plot->graph(2)->addData(x_ref[6], x_ref[7]);
    ui->angle_plot->graph(0)->addData(time, anglef_);
    ui->angle_plot->graph(1)->addData(time, x[0]);
    ui->angle_plot->graph(2)->addData(time, x_ref[0]);
    ui->control_plot->graph(0)->addData(time, u[0]);
    ui->control_plot->graph(1)->addData(time, u_ref[0]);
//    ui->position_plot->rescaleAxes(true);
    ui->position_plot->replot();
    ui->angle_plot->rescaleAxes(true);
    ui->angle_plot->replot();
    ui->control_plot->rescaleAxes(true);
    ui->control_plot->replot();

    last_traj_x = x[6];
    last_ref_x = x_ref[6];
}

MainWindow::~MainWindow()
{
    delete ui;
}
