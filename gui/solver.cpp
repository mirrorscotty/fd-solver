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

extern double R;

Solver::Solver(QWidget *parent)
{
    setupUi(this); // Initialize the user interface

    // Connect all the relevant signals
    connect(buttonSolve, SIGNAL( clicked() ), this, SLOT( solveProblems() ));

    connect(actionQuit, SIGNAL( triggered() ), this, SLOT( quitApplication() ));
    connect(actionSave_Simulation, SIGNAL( triggered() ), this, SLOT( saveSimulation() ));
    connect(actionOpen, SIGNAL( triggered() ), this, SLOT( loadSimulation() ));
    connect(actionAbout, SIGNAL( triggered() ), this, SLOT( about() ));
    connect(actionSave_CSV, SIGNAL( triggered() ), this, SLOT( saveCSV() ));

    connect(comboLeftBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeLBC(int) ));
    connect(comboRightBC, SIGNAL( currentIndexChanged(int) ), this, SLOT( changeRBC(int) ));

    connect(spinNodes, SIGNAL( valueChanged(int) ), this, SLOT( setMaxNode(int) ));
    connect(spinDt, SIGNAL( valueChanged(double) ), this, SLOT(setMaxTIndex(double) ));
    connect(spintEnd, SIGNAL( valueChanged(double) ), this, SLOT(setMaxTIndex(double) ));

    // Set the domain and datalist variables to NULL. They're both initialized
    // later in the program.
    domain = NULL;
    datalist = NULL;

    // Since the default boundary condition in the program is "insulation,"
    // none of the options need to be shown initially.
    leftBCHideAll();
    rightBCHideAll();

    // Check one of the radio buttons for the graph
    radioTime->setChecked(true);

    // Setup curves for each of the variables to plot. Also, make them
    // different colors.
    Temp = new QwtPlotCurve( "Temperature" );
    Temp->setPen( QPen( Qt::red ) );
    Prod = new QwtPlotCurve( "Product Concentration" );
    Prod->setPen( QPen( Qt::blue ) );
    Bact = new QwtPlotCurve( "Bacteria Concentration" );
    Bact->setPen( QPen( Qt::green ) );
    alpha = new QwtPlotCurve("Thermal Diffusivity");
    alpha->setPen( QPen( Qt::magenta ) );

    // Initialize the GUI
    comboLeftBC->setCurrentIndex(0);
    comboLeftBC->setEnabled(false);
    comboRightBC->setCurrentIndex(2);
    comboRightBC->setEnabled(false);
}

Solver::~Solver()
{
    DestroyDomain1D(domain);
    delete Temp;
    delete Prod;
    delete Bact;
    delete alpha;

    return;
}

// Close the window and exit the program.
void Solver::quitApplication()
{
    qApp->exit(0);
}

// Set the maximum value for the spin box that allows the user to select the
// node to plot data for.
void Solver::setMaxNode(int i) {
    spinPlotNode->setMaximum(i-1);
    return;
}

void Solver::setMaxTIndex(double i) {
    spinPlotTIndex->setMaximum((int) spintEnd->value()/spinDt->value()-1);
    return;
}

// When the user changes the boundary condition settings, show/hide options as
// appropriate.
void Solver::changeLBC(int index)
{
    // Hide everything to begin with.
    leftBCHideAll();
    switch(index) {
        case 1:
            // Show just options for convection.
            leftBCShowConv();
            break;
        case 2:
            // Show options for convection and heating/cooling.
            leftBCShowConvTime();
            break;
    }
    return;
}

// Same as above.
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

// Function to hide all the widget associated with the left boundary condition.
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

// Same for the right BC.
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

// Show options associated with convection.
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

// Show options for convection where the external temperature can take on two
// values.
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

// Macro to simplify the code for saving data from the gui to a linked list.
#define GETVARIABLE( VARNAME, BOXNAME ) { tmp->name = #VARNAME; tmp->value = (BOXNAME)->value(); tmp->next = new_var(); tmp = tmp->next; }
// Store all (most of) the data that the user has inputted in the gui in the
// datalist variable.
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
    GETVARIABLE(HConv, spinHRight)
    GETVARIABLE(CAinit, spinC1Init)
    GETVARIABLE(CBinit, spinC2Init)
    
    tmp->name = "Deltax";
    tmp->value = spinWidth->value()/spinNodes->value();
    tmp->next = new_var();
    tmp = tmp->next;

    tmp->name = "NTimeSteps";
    tmp->value = (int) (spintEnd->value()/spinDt->value());

    return;
}

// Simplify code for moving values from a linked list into the gui.
#define RESTOREVAR( VARNAME, BOXNAME ) { if(strcmp(tmp->name, #VARNAME) == 0){(BOXNAME)->setValue(tmp->value);} }
// Take all the data from the datalist variable and stick it into the
// appropriate fields in the gui.
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
        RESTOREVAR(HConv, spinHRight)
        RESTOREVAR(CAinit, spinC1Init)
        RESTOREVAR(CBinit, spinC1Init)

        if(strcmp(tmp->name, "R") == 0) {
            puts("Success!");
            R = tmp->value;
        }

        tmp = tmp->next;
    }
    return;
}

// Take the information about the domain from the appropriate fields in the gui
// and initialize the "domain" variable. Three dependant variables are solved
// for: Temperature, product concentration, and bacteria concentration.
void Solver::setupDomain()
{
    struct Node1D *node;

    // Load the required parameters
    int NNodes = spinNodes->value();
    double Deltax = spinWidth->value()/NNodes;
    double Deltat = spinDt->value();
    int NTimeSteps = (int) (spintEnd->value()/Deltat);

    // Set up the domain variable.
    domain = CreateDomain1D("x", NNodes, Deltax, Deltat, 3, NTimeSteps);

    // Apply the appropriate initial conditions.
    DomainApplyInitialCondition(domain, 0, spinTInit->value());
    DomainApplyInitialCondition(domain, 1, spinC1Init->value());
    DomainApplyInitialCondition(domain, 2, spinC2Init->value());

    // Define the update functions for each variable at all the nodes.
    node = domain->Nodes[0];
    while(node) {
        NodeSetUpdateFunction(node, 0, 'T', &UpdateSubdomainCyl);
        NodeSetUpdateFunction(node, 1, 'c', &UpdateSubdomainRxn1);
        NodeSetUpdateFunction(node, 2, 'd', &UpdateSubdomainRxn2);
        node = node->Next;
    }

    // Apply the boundary conditions. Currently, these are hardcoded, but they
    // should be configurable from the gui.
    NodeSetUpdateFunction(domain->Nodes[0], 0, 'T', &UpdateRadialSymmetryBoundary);
    NodeSetUpdateFunction(domain->Nodes[NNodes-1], 0, 'T', &UpdateConvectiveBoundaryCyl);

    return;
}

void Solver::solveProblems()
{
    struct var *tmp;
    int nodenum = 0; // Node to use to plot stuff
    int timeindex = 0;
    
    // Set all the global variables to 0 so that things don't break hs orribly
    //  if everything wasn't defined.
    initialize_variables();

    // Set the datalist to NULL since we don't want to deal with old values.
    // This is a horrible idea, and the old data should be cleaned up instead.
    datalist = NULL;

    storeVariables();

    // Save the data to the global variables used by the material property
    // funcitons. (Also bad.)
    tmp = datalist;
    while(tmp) {
        store_data(tmp);
        tmp = tmp->next;
    }

    // Clean up stuff that was there before.
    DestroyDomain1D(domain);
    // Initialize the domain
    setupDomain();

    // Solve
    while(UpdateDomain(domain) != 1)
        progressBar->setValue( (int) domain->TimeIndex/domain->NumTimes * 100);

    // Plot the data for node 0. This should be configureable, but isn't.
    if(!radioTime->isChecked()) {
        nodenum = spinPlotNode->value();
        plotResultsTime(domain->Nodes[nodenum]);
    } else {
        timeindex = spinPlotTIndex->value();
        plotResultsSpace(timeindex);
    }
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

// Function to plot results at a single node as a function of time.
void Solver::plotResultsTime(struct Node1D *n)
{
    int i;
    int npts = n->NumValues;
    double *t, *T, *c, *d, *a, tmp;
    char *title;
    title = (char*) calloc(20, sizeof(char));
    
    /* Create an array with the time values. */
    t = (double*) calloc(npts, sizeof(double));
    for(i=0; i<npts; i++) {
        t[i] = i*n->dt;
    }

    // Detach the curves from the plot so that they don't show up. If we want
    // them to show up, we'll reattach them later.
    Temp->detach();
    Prod->detach();
    Bact->detach();
    alpha->detach();

    // Plot temperature data if the box for it is checked.
    if(checkTemp->isChecked()) {
        T = (double*) calloc(npts, sizeof(double));
        // Get the values from the solution.
        for(i=0; i<npts; i++) {
            T[i] = n->Value[0][i];
        }

        // Save the data to the curve. This makes a copy of the data, so we can
        // delete the two arrays later.
        Temp->setSamples(t, T, npts);

        // Put the curve on the graph.
        Temp->attach(qwtPlot);

        // Delete the temperature data we retrieved.
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
    
    if(checkk->isChecked()) {
        a = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            tmp = n->Value[0][i];
            a[i] = k(tmp)/(rho(tmp)*Cp(tmp));
        }

        alpha->attach(qwtPlot);

        alpha->setSamples(t, a, npts);
        free(a);
    }

    sprintf(title, "x = %g", n->dx*n->NodeNum);
    qwtPlot->setTitle(title);
    free(title);

    // Set the axis title. 2 corresponds to the lower x axis.
    qwtPlot->setAxisTitle(2, "Time (sec)");

    // Update the graph to show what we just did.
    qwtPlot->replot();

    // Get rid of the "t" array. We don't want it anymore.
    free(t);
}

// Plot the results at a fixed time as a function of the spatial coordinate.
void Solver::plotResultsSpace(int t)
{
    int i;
    int npts = domain->NumNodes;
    double *x, *T, *c, *d, *a, tmp;
    char *title;
    title = (char*) calloc(20, sizeof(char));
    
    /* Create an array with the x values. */
    x = (double*) calloc(npts, sizeof(double));
    for(i=0; i<npts; i++) {
        x[i] = i*domain->Nodes[i]->dx;
    }

    // Detach the curves from the plot so that they don't show up. If we want
    // them to show up, we'll reattach them later.
    Temp->detach();
    Prod->detach();
    Bact->detach();
    alpha->detach();

    // Plot temperature data if the box for it is checked.
    if(checkTemp->isChecked()) {
        T = (double*) calloc(npts, sizeof(double));
        // Get the values from the solution.
        for(i=0; i<npts; i++) {
            T[i] = domain->Nodes[i]->Value[0][t];
        }

        // Save the data to the curve. This makes a copy of the data, so we can
        // delete the two arrays later.
        Temp->setSamples(x, T, npts);

        // Put the curve on the graph.
        Temp->attach(qwtPlot);

        // Delete the temperature data we retrieved.
        free(T);
    }

    if(checkProduct->isChecked()) {
        c = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            c[i] = domain->Nodes[i]->Value[1][t];
        }

        Prod->attach(qwtPlot);

        Prod->setSamples(x, c, npts);
        free(c);
    }

    if(checkBact->isChecked()) {
        d = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            d[i] = domain->Nodes[i]->Value[2][t];
        }

        Bact->attach(qwtPlot);

        Bact->setSamples(x, d, npts);
        free(d);
    }
    
    if(checkk->isChecked()) {
        a = (double*) calloc(npts, sizeof(double));
        for(i=0; i<npts; i++) {
            tmp = domain->Nodes[i]->Value[0][i];

            a[i] = k(tmp)/(rho(tmp)*Cp(tmp));
        }

        alpha->attach(qwtPlot);

        alpha->setSamples(x, a, npts);
        free(a);
    }

    // Set the title.
    sprintf(title, "t = %g", t*domain->dt);
    qwtPlot->setTitle(title);
    free(title);

    // Set the axis title. 2 corresponds to the lower x axis.
    qwtPlot->setAxisTitle(2, "Radius (m)");

    // Update the graph to show what we just did.
    qwtPlot->replot();

    // Get rid of the "x" array. We don't want it anymore.
    free(x);
}
// Load all the simulation data from a text file.
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
    
    // Read the data file into a buffer.
    buffer = read_datafile(ba.data());

    // Parse the file line by line. (Hopefully there's not more than 200 lines.)
    for(i=0; i<200; i++) {
        // Strip any comments that might be lurking in the file.
        buffer[i] = remove_comments(buffer[i]);
        // Parse the line and save it to a variable.
        tmp = read_line(buffer[i]);
        // If the variable actually contained something, save it. Otherwise just
        // keep going and pretend like it never existed.
        // TODO: Fix the memory leak.
        if(strcmp(tmp->name, "NULL") != 0)
            datalist = push_var(datalist, tmp);
    }

    // Clean up the buffer.
    delete_buffer(buffer);

    // Load the variables into the interface.
    loadVars();
    return;
}

// Save everything the user was working on to a file.
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

    // Open up the requested file.
    fp = fopen( (char*) ba.data(), "w+" );
    // Check for errors in opening the file.
    // TODO: Pop up an error message instead of printing stuff to the console.
    if(!fp) {
        fprintf(stderr, "Failed to open file for writing.");
        return;
    }

    //while(tmp) {
    //    fprintf(fp, "%s = %g\n", tmp->name, tmp->value);
    //    tmp = tmp->next;
    //}

    // Close the file.
    fclose(fp);

    return;
}

void Solver::saveCSV()
{
    QString path;
    QByteArray ba;

    // Check to see if the user has run the simulation before trying to save
    // the output. Stops the program from seg faulting.
    if(!domain) {
        QMessageBox::warning(this, "Error", "Please run the simulation before saving the results.");
        return;
    }

    // Get the filename to save to.
    path = QFileDialog::getSaveFileName(
            this,
            "Save As",
            QString::null,
            QString::null);

    // Convert the filename to a char*
    ba = path.toLocal8Bit();

    if(radioTime->isChecked()) {
        CSVOutFixedTime(domain, spinPlotTIndex->value(), ba.data());
    } else {
        CSVOutFixedNode(domain->Nodes[spinPlotNode->value()], ba.data());
    }
    return;
}
        


// Bring up an about dialog.
void Solver::about()
{
    QMessageBox::about(this, "About", "Simple finite difference solver.\n");
}

