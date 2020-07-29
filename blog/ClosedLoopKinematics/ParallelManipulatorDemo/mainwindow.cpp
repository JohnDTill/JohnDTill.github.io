#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <math.h>
#include <QGraphicsLineItem>
#include <QGraphicsView>
#include <QMatrix4x4>
#include <QVector2D>
#include <QVector4D>

static constexpr float l1 = 1;
static constexpr float l2 = 2;
static constexpr float X1 = 0;
static constexpr float Y1 = 0;
static constexpr float X3 = -1;
static constexpr float Y3 = 0;

static float Xe = -0.5;
static float Ye = 2;
static QVector4D guess(3.141f/2, 3.141f/2, 0, 0);

static QVector4D objFunc(const QVector4D& guess){
    float alpha = guess[0];
    float beta = guess[1];
    float theta = guess[2];
    float phi = guess[3];

    return QVector4D(
                    Xe - (X1 + l1*cos(alpha) + l2*cos(alpha+theta)),
                    Ye - (Y1 + l1*sin(alpha) + l2*sin(alpha+theta)),
                    Xe - (X3 + l1*cos(beta) + l2*cos(beta+phi)),
                    Ye - (Y3 + l1*sin(beta) + l2*sin(beta+phi))
                );
}

static QMatrix4x4 getJacobian(QVector4D& guess, const QVector4D& err){
    //Approximate the Jacobian of the objective function w.r.t. the guess via finite differences

    static constexpr float delta = 1e-4f;
    QMatrix4x4 J;

    for(int i = 0; i < 4; i++){
        guess[i] += delta;
        const QVector4D finite_difference = (objFunc(guess) - err) / delta; //First order forward difference
        J.setColumn(i, finite_difference);
        guess[i] -= delta;
    }

    return J;
}

static QVector4D solveLevenbergMarquardt(QVector4D y0, bool& success){
    //Convex optimization routine- this is similar to the logic fsolve uses.

    static constexpr float sos_tol = 1e-12f;
    static constexpr int max_iter = 500;
    static float damp = 1e-2f;
    static float alm_adptv_coeff = 0.5;
    success = true;

    QVector4D r = objFunc(y0);
    float E = r.lengthSquared();
    if( E <= sos_tol ) return y0;

    QVector4D y  = y0;
    QMatrix4x4 I;
    I.setToIdentity();
    int counter = 0;

    while( E > sos_tol ){
        if(counter++ > max_iter){
            success = false;
            return y;
        }

        QMatrix4x4 J = getJacobian(y,r);
        QMatrix4x4 lhs = J.transposed()*J + damp*I;
        QVector4D rhs = -J.transposed()*r;

        damp *= alm_adptv_coeff;
        //lhs is symmetric, positive semi-definite, positive definite for damp > 0
        QVector4D pc = lhs.inverted() * rhs;
        r = objFunc(y + pc);
        float Ec = r.lengthSquared();

        if( Ec < E ){
            y += pc;
            E = Ec;
        }else{
            while( Ec > E ){
                if(counter++ > max_iter){
                    success = false;
                    return y;
                }

                float old_damp = damp;
                damp /= alm_adptv_coeff;

                lhs += (damp - old_damp)*I;
                pc = lhs.inverted() * rhs;
                r = objFunc(y + pc);
                Ec = r.lengthSquared();
            }

            y += pc;
            E = Ec;
        }
    }

    return y;
}

static constexpr float max_reach = l1+l2;
static float min_reach = abs(l2-l1);

static constexpr float scale = 100;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow){
    ui->setupUi(this);

    //Set up the chart view
    chart.setSceneRect(-2.5*scale, -3*scale, 4*scale, 4*scale);
    QGraphicsView* view = new QGraphicsView(&chart);
    view->setFocusProxy(this);
    ui->verticalLayout->addWidget(view);

    QPen link_pen(Qt::blue);
    link_pen.setWidthF(1.5);
    QPen joint_pen(Qt::red);
    joint_pen.setWidthF(2.0);

    for(int i = 0; i < 4; i++){
        links[i] = new QGraphicsLineItem();
        links[i]->setPen(link_pen);
        chart.addItem(links[i]);
    }

    for(int i = 0; i < 5; i++){
        joints[i] = new QGraphicsEllipseItem();
        joints[i]->setPen(joint_pen);
        joints[i]->setRect(-3,-3,6,6);
        chart.addItem(joints[i]);
    }

    update();
}

MainWindow::~MainWindow(){
    delete ui;
}

void MainWindow::update(){
    //Check for workspace violations- can do this ahead of time since robot has simple design
    QVector2D pe(Xe,Ye);
    QVector2D p1(X1,Y1);
    QVector2D p3(X3,Y3);

    const float distance_from_p1 = (pe-p1).length();
    const float distance_from_p3 = (pe-p3).length();

    if( distance_from_p1 > max_reach ||
        distance_from_p3 > max_reach ||
        distance_from_p1 < min_reach ||
        distance_from_p3 < min_reach ){

        ui->label->setText("<em>Workspace Violation</em><br>Goal: (" + QString::number(Xe) + ", " + QString::number(Ye) + ")");
        return;
    }else{
        ui->label->setText("<b>Five-bar Robot</b><br>Goal: (" + QString::number(Xe) + ", " + QString::number(Ye) + ")");
    }

    //Solve the IK and update the chart
    Xe = pe.x();
    Ye = pe.y();

    bool success;
    guess = solveLevenbergMarquardt(guess, success);
    if(!success){
        ui->label->setText("<em>SOLVER FAILURE</em><br>Goal: (" + QString::number(Xe) + ", " + QString::number(Ye) + ")");
        return;
    }

    float alpha = guess[0];
    float beta = guess[1];

    float X2 = X1 + l1*cos(alpha);
    float Y2 = Y1 + l1*sin(alpha);
    float X4 = X3 + l1*cos(beta);
    float Y4 = Y3 + l1*sin(beta);

    links[0]->setLine(X1*scale, -Y1*scale, X2*scale, -Y2*scale);
    links[1]->setLine(X2*scale, -Y2*scale, Xe*scale, -Ye*scale);
    links[2]->setLine(Xe*scale, -Ye*scale, X4*scale, -Y4*scale);
    links[3]->setLine(X4*scale, -Y4*scale, X3*scale, -Y3*scale);

    joints[0]->setPos(X1*scale,-Y1*scale);
    joints[1]->setPos(X2*scale,-Y2*scale);
    joints[2]->setPos(Xe*scale,-Ye*scale);
    joints[3]->setPos(X3*scale,-Y3*scale);
    joints[4]->setPos(X4*scale,-Y4*scale);
}

void MainWindow::keyPressEvent(QKeyEvent* event){
    //Respond to arrow and WASD keys

    constexpr float incr = 3e-2f;

    switch(event->key()){
        case Qt::Key_A: Xe -= incr; break;
        case Qt::Key_Left: Xe -= incr; break;
        case Qt::Key_D: Xe += incr; break;
        case Qt::Key_Right: Xe += incr; break;
        case Qt::Key_S: Ye -= incr; break;
        case Qt::Key_Down: Ye -= incr; break;
        case Qt::Key_W: Ye += incr; break;
        case Qt::Key_Up: Ye += incr; break;
    }

    update();
}
