/*
 * Copyright (C) 2006 Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <math.h>
#include <model_cell.h>

/*
 * Model Functions
 */


extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new Cell();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "I",
        "A",
        DefaultGUIModel::INPUT,
    },
    {
        "Vm",
        "V",
        DefaultGUIModel::OUTPUT,
    },
    {
        "V0",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "Cm",
        "pF",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "R_L_MC",
        "GOhm - Leak of model cell",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "Iapp_offset",
        "pA - Current added to the output.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "rate",
        "Hz - The rate of integration.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER,
    },
};

static size_t num_vars = sizeof(vars)/sizeof(DefaultGUIModel::variable_t);

/*
 * Macros for making the code below a little bit cleaner.
 */

#define Iapp (input(0)*1e12) // Convert from A to pA  
#define Vm (y[0]) // mV
#define dVm (dydt[0])

Cell::Cell(void)
    : DefaultGUIModel("Model cell",::vars,::num_vars) {
    createGUI(vars, num_vars);
    /*
     * Initialize Parameters
     */
    V0 = 0.0;
    Cm = 33.0;
    R_L_MC = 0.5181;
    Iapp_offset = 0.0;
    rate = 40000;

    /*
     * Initialize Variables
     */
    Vm = V0;
    period = RT::System::getInstance()->getPeriod()*1e-9;
    steps = static_cast<int>(ceil(period*rate));

    /*
     * Initialize GUI
     */
    setParameter("V0",V0);
    setParameter("Cm",Cm);
    setParameter("R_L_MC",R_L_MC);
    setParameter("Iapp_offset",Iapp_offset);
    setParameter("rate",rate);
    refresh();
}

Cell::~Cell(void) {}

/*
 * Simple Euler solver.
 */

void Cell::solve(double dt, double *y) {
    double dydt[MODEL_DIM];

    derivs(y,dydt);

    for(size_t i = 0;i < MODEL_DIM;++i)
        y[i] += dt*dydt[i];
}

void Cell::derivs(double *y,double *dydt) {
    dVm = (Iapp + Iapp_offset - Vm/R_L_MC)/Cm;
}

void Cell::execute(void) {

    /*
     * Because the real-time thread may run much slower than we want to
     *   integrate we need to run multiple interations of the solver.
     */

    for(int i = 0;i < steps;++i)
        solve(period/steps,y);

    // Membrane capacitance Cm is added to enable the use of nS/pF
    output(0) = Vm*1e-3; //convert from mV to V
}

void Cell::update(DefaultGUIModel::update_flags_t flag) {
    switch(flag) {
        case MODIFY:
            V0 = getParameter("V0").toDouble();
            Cm = getParameter("Cm").toDouble();
            R_L_MC = getParameter("R_L_MC").toDouble();
            Iapp_offset = getParameter("Iapp_offset").toDouble();
            rate = getParameter("rate").toDouble();
            steps = static_cast<int>(ceil(period*rate));
    
            Vm = V0;
        break;
        case PERIOD:
            period = RT::System::getInstance()->getPeriod()*1e-9; // ms
            steps = static_cast<int>(ceil(period*rate));
        break;
        case PAUSE:
            output(0) = V0*1e-3; 
        break;
        case UNPAUSE:
        break;
        default:
        break;
    }
}
