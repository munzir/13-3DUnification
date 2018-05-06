#include "autorallyrtwindow.h"
#include "ui_autorallyrtwindow.h"

AutoRallyRTWindow::AutoRallyRTWindow(Scalar dt, QWidget *parent):
    dt_(dt),
    QMainWindow(parent),
    ui(new Ui::AutoRallyRTWindow)
{
    ui->setupUi(this);
    position = new QCPCurve(ui->pos_plot->xAxis, ui->pos_plot->yAxis);
    ui->pos_plot->xAxis->setLabel("x Position");
    ui->pos_plot->xAxis->setRangeLower(-25.0);
    ui->pos_plot->xAxis->setRangeUpper(5.0);
    ui->pos_plot->yAxis->setLabel("y Position");
    ui->pos_plot->yAxis->setRangeLower(-12.0);
    ui->pos_plot->yAxis->setRangeUpper(12.0);
//    ui->pos_plot->addPlottable(position);

    // Plot velocity
    ui->vel_plot->addGraph();
    ui->vel_plot->graph(0)->setName("Velocity");
    ui->vel_plot->xAxis->setLabel("Time (s)");
    ui->vel_plot->yAxis->setLabel("Vx (m/s)");
    ui->vel_plot->rescaleAxes(true);

    // Plot control
    ui->control_plot->addGraph();
    ui->control_plot->addGraph();
    ui->control_plot->graph(0)->setName("Steering");
    ui->control_plot->graph(1)->setName("Throttle");
    ui->control_plot->graph(1)->setPen(QPen(Qt::red));
    ui->control_plot->legend->setVisible(true);
    ui->control_plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop|Qt::AlignRight);
    ui->control_plot->xAxis->setLabel("Time (s)");
    ui->control_plot->yAxis->setLabel("Control Magnitude");
    ui->control_plot->yAxis->setRangeLower(-2.0);
    ui->control_plot->yAxis->setRangeUpper(2.0);
    ui->control_plot->xAxis->rescale(true);
}

void AutoRallyRTWindow::update_graph(const QVector<Scalar> &state, const QVector<Scalar> &control)
{
    static Scalar time = 0.0;
    time += dt_;
    position->addData(state[0], state[1]);
    ui->pos_plot->replot();
    ui->vel_plot->graph(0)->addData(time, state[4]);
    ui->vel_plot->rescaleAxes(true);
    ui->vel_plot->replot();
    ui->control_plot->graph(0)->addData(time, control[0]);
    ui->control_plot->graph(1)->addData(time, control[1]);
    ui->control_plot->xAxis->rescale(true);
    ui->control_plot->replot();
}

AutoRallyRTWindow::~AutoRallyRTWindow()
{
    delete ui;
    if(position) delete position;
}
