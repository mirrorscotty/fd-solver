#include <QtGui>

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_math.h>

#include <stdio.h>
#include <string.h>

#include "solver.h"
#include "fd-solver.h"
#include "datafile.h"
#include "can.h"

Solver::Solver(QWidget *parent)
{
    setupUi(this); // Initialize the user interface

    // Connect all the relevant signals
    connect(buttonSolve, SIGNAL( clicked() ), this, SLOT( solveProblems() ));
    connect(actionQuit, SIGNAL( triggered() ), this, SLOT( quitApplication() ));
    connect(actionSave_Simulation, SIGNAL( triggered() ), this, SLOT( saveSimulation() ));
    connect(actionOpen, SIGNAL( triggered() ), this, SLOT( loadSimulation() ));
    connect(actionAbout, SIGNAL( triggered() ), this, SLOT( about() ));
    connect(comboLeftBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeLBC(int) ));
    connect(comboRightBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeRBC(int) ));

    domain = NULL;
    datalist = NULL;

    leftBCHideAll();
    rightBCHideAll();

    radioTime->setChecked(true);
        Temp = new QwtPlotCurve( "Temperature" );
        Temp->setPen( QPen( Qt::red ) );
        Prod = new QwtPlotCurve( "Product Concentration" );
        Prod->setPen( QPen( Qt::blue ) );
        Bact = new QwtPlotCurve( "Bacteria Concentration" );
        Bact->setPen( QPen( Qt::green ) );
}

Solver::~Solver()
{
    return;
}

void Solver::solve()
{
    return;
}

void Solver::quitApplication()
{
    qApp->exit(0);
}

void Solver::changeLBC(int index)
{
    leftBCHideAll();
    switch(index) {
        case 1:
            leftBCShowConv();
            break;
        case 2:
            leftBCShowConvTime();
            break;
    }
    return;
}

void Solver::changeRBC(int index)
{
    rightBCHideAll();
    switch(index) {
        case 1:
            rightBCShowConv();
            break;
        case 2:
            rightBCShowConvTime();
            break;
    }
    return;
}

void Solver::leftBCHideAll()
{
    lbch->hide();
    lbctext->hide();
    lbctexthot->hide();
    lbctextcold->hide();
    lbctheat->hide();
    spinHLeft->hide();
    spinTExtLeft->hide();
    spinTExtHotLeft->hide();
    spinTExtColdLeft->hide();
    spinTHeatLeft->hide();
    return;
}

void Solver::rightBCHideAll()
{
    rbch->hide();
    rbctext->hide();
    rbctexthot->hide();
    rbctextcold->hide();
    rbctheat->hide();
    spinHRight->hide();
    spinTExtRight->hide();
    spinTExtHotRight->hide();
    spinTExtColdRight->hide();
    spinTHeatRight->hide();

    return;
}

void Solver::leftBCShowConv()
{
    lbch->show();
    lbctext->show();
    spinHLeft->show();
    spinTExtLeft->show();
 
    return;
}

void Solver::rightBCShowConv()
{
    rbch->show();
    rbctext->show();
    spinHRight->show();
    spinTExtRight->show();
    return;
}

void Solver::rightBCShowConvTime()
{
    rbch->show();
    rbctexthot->show();
    rbctextcold->show();
    rbctheat->show();
    spinHRight->show();
    spinTExtHotRight->show();
    spinTExtColdRight->show();
    spinTHeatRight->show();

    return;
}

void Solver::leftBCShowConvTime()
{
    lbch->show();
    lbctexthot->show();
    lbctextcold->show();
    lbctheat->show();
    spinHLeft->show();
    spinTExtHotLeft->show();
    spinTExtColdLeft->show();
    spinTHeatLeft->show();
    return;
}

#define GETVARIABLE( VARNAME, BOXNAME ) { tmp->name = #VARNAME; tmp->value = (BOXNAME)->value(); tmp->next = new_var(); tmp = tmp->next; }
void Solver::storeVariables()
{
    struct var *tmp;
    datalist = new_var();
    tmp = datalist;

    GETVARIABLE(Mpro, spinMPro)
    GETVARIABLE(Mfat, spinMFat)
    GETVARIABLE(Mcar, spinMCar)
    GETVARIABLE(Mfib, spinMFib)
    GETVARIABLE(Mash, spinMAsh)
    GETVARIABLE(Mwat, spinMWat)
    GETVARIABLE(Mice, spinMIce)
    GETVARIABLE(AA, spinA1)
    GETVARIABLE(EaA, spinEa1)
    GETVARIABLE(AB, spinA2)
    GETVARIABLE(EaB, spinEa2)
    GETVARIABLE(To, spinTInit)
    GETVARIABLE(Text_hot, spinTExtHotRight)
    GETVARIABLE(Text_cold, spinTExtColdRight)
    GETVARIABLE(NNodes, spinNodes)
    GETVARIABLE(Deltat, spinDt)
    GETVARIABLE(t_heat, spinTHeatRight)
    
    tmp->name = "Deltax";
    tmp->value = spinWidth->value()/spinNodes->value();
    tmp->next = new_var();
    tmp = tmp->next;

    tmp->name = "NTimeSteps";
    tmp->value = (int) (spintEnd->value()/spinDt->value());

    return;
}

#define RESTOREVAR( VARNAME, BOXNAME ) { if(strcmp(tmp->name, #VARNAME) == 0){(BOXNAME)->setValue(tmp->value);} }
void Solver::loadVars()
{
    struct var *tmp;
    tmp = datalist;
    while(tmp) {
        RESTOREVAR(Mpro, spinMPro)
        RESTOREVAR(Mfat, spinMFat)
        RESTOREVAR(Mcar, spinMCar)
        RESTOREVAR(Mfib, spinMFib)
        RESTOREVAR(Mash, spinMAsh)
        RESTOREVAR(Mwat, spinMWat)
        RESTOREVAR(Mice, spinMIce)
        RESTOREVAR(AA, spinA1)
        RESTOREVAR(EaA, spinEa1)
        RESTOREVAR(AB, spinA2)
        RESTOREVAR(EaB, spinEa2)
        RESTOREVAR(To, spinTInit)
        RESTOREVAR(Text_hot, spinTExtHotRight)
        RESTOREVAR(Text_cold, spinTExtColdRight)
        RESTOREVAR(NNodes, spinNodes)
        RESTOREVAR(Deltat, spinDt)
        RESTOREVAR(t_heat, spinTHeatRight)

        tmp = tmp->next;
    }
    return;
}

void Solver::setupDomain()
{
    struct Node1D *node;
    int NNodes = spinNodes->value();
    double Deltax = spinWidth->value()/NNodes;
    double Deltat = spinDt->value();
    int NTimeSteps = (int) (spintEnd->value()/Deltat);

    domain = CreateDomain1D("x", NNodes, Deltax, Deltat, 3, NTimeSteps);

    DomainApplyInitialCondition(domain, 0, spinTInit->value());
    DomainApplyInitialCondition(domain, 1, spinC1Init->value());
    DomainApplyInitialCondition(domain, 2, spinC2Init->value());

    node = domain->Nodes[0];
    while(node) {
        NodeSetUpdateFunction(node, 0, 'T', &UpdateSubdomain);
        NodeSetUpdateFunction(node, 1, 'c', &UpdateSubdomainRxn1);
        NodeSetUpdateFunction(node, 2, 'd', &UpdateSubdomainRxn2);
        node = node->Next;
    }

    NodeSetUpdateFunction(domain->Nodes[0], 0, 'T', &UpdateInsulatedBoundary);
    NodeSetUpdateFunction(domain->Nodes[NNodes-1], 0, 'T', &UpdateConvectiveBoundary);

    printf("Step Size: %g\n", Deltat);
    printf("End Time: %g\n", spintEnd->value());
    printf("Time Steps: %d\n", NTimeSteps);

    return;
}

void Solver::solveProblems()
{
    struct var *tmp;
    initialize_variables();
    //init("can_data.dat");
    datalist = NULL;
    storeVariables();
    tmp = datalist;
    while(tmp) {
        store_data(tmp);
        tmp = tmp->next;
    }

    setupDomain();

    // Solve
    while(UpdateDomain(domain) != 1);

    plotResultsTime(domain->Nodes[0]);

    DestroyDomain1D(domain);
}

/* Example, since the Qwt docs are terrible

void Solver::plotResultsTime(struct Node1D* n)
{
    // Insert new curves
    QwtPlotCurve *cSin = new QwtPlotCurve( "y = sin(x)" );
    cSin->setPen( QPen( Qt::red ) ); // Draw it in red
    cSin->attach(qwtPlot); // Stick it in the plot

    // Create test data
    int npoints = 100;
    double *x, *y;
    x = (double*) calloc(npoints, sizeof(double));
    y = (double*) calloc(npoints, sizeof(double));

    int i;
    for(i=0; i<npoints; i++) {
        x[i] = i*3.14159/npoints;
        y[i] = sin(x[i]);
    }

    cSin->setRawSamples(x, y, npoints);

    // Show the plots
    qwtPlot->replot();
}
*/

void Solver::plotResultsTime(struct Node1D *n)
{
    int i;
    int npts = n->NumValues;
    double *t, *T, *c, *d;
    
    /* Create an array with the time values. */
    t = (double*) calloc(npts, sizeof(double));
    for(i=0; i<npts; i++) {
        t[i] = i*n->dt;
    }

    Temp->detach();
    Prod->detach();
    Bact->detach();

    if(checkTemp->isChecked()) {
        T = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            T[i] = n->Value[0][i];
        }

        Temp->attach(qwtPlot);

        Temp->setSamples(t, T, npts);
        free(T);
    }

    if(checkProduct->isChecked()) {
        c = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            c[i] = n->Value[1][i];
        }

        Prod->attach(qwtPlot);

        Prod->setSamples(t, c, npts);
        free(c);
    }

    if(checkBact->isChecked()) {
        d = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            d[i] = n->Value[2][i];
        }

        Bact->attach(qwtPlot);

        Bact->setSamples(t, d, npts);
        free(d);
    }

    qwtPlot->replot();
}

void Solver::loadSimulation()
{
    struct var *tmp;
    QString path;
    QByteArray ba;
    char **buffer;
    int i;
    
    buffer = NULL;
    tmp = NULL;

    // Get the filename to save to.
    path = QFileDialog::getOpenFileName(
            this,
            "Open",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();

    buffer = read_datafile(ba.data());

    for(i=0; i<200; i++) {
        buffer[i] = remove_comments(buffer[i]);
        tmp = read_line(buffer[i]);
        if(strcmp(tmp->name, "NULL") != 0)
            datalist = push_var(datalist, tmp);
    }

    delete_buffer(buffer);

    printf("something: %s: %g\n", datalist->next->name, datalist->next->value);

    loadVars();
    return;
}

void Solver::saveSimulation()
{
    struct var *tmp;
    QString path;
    QByteArray ba;
    FILE *fp;
    storeVariables();
    tmp = datalist;


    // Get the filename to save to.
    path = QFileDialog::getSaveFileName(
            this,
            "Save As",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();

    fp = fopen( (char*) ba.data(), "w+" );
    if(!fp) {
        fprintf(stderr, "Failed to open file for writing.");
        return;
    }

    while(tmp) {
        fprintf(fp, "%s = %g\n", tmp->name, tmp->value);
        tmp = tmp->next;
    }

    fclose(fp);

    return;
}

void Solver::about()
{
    QMessageBox::about(this, "About", "Simple finite difference solver.\n");
}

