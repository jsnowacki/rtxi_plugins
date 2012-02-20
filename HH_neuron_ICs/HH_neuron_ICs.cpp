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
#include <HH_neuron_ICs.h>

/*
 * Model Functions
 */

static inline double alpha_m(double V) {
    double x = -(V+40.0);
    double y = 10.0;

    if(fabs(x/y) < 1e-6)
        return 0.1*y*(1-x/y/2.0);
    else
        return 0.1*x/(exp(x/y)-1.0);
}

static inline double beta_m(double V) {
    return 4.0*exp(-(V+65.0)/18.0);
}

static inline double m_inf(double V) {
    return alpha_m(V)/(alpha_m(V)+beta_m(V));
}

static inline double tau_m(double V) {
    return 1.0/(alpha_m(V)+beta_m(V));
}

static inline double alpha_h(double V) {
    return 0.07*exp(-(V+65.0)/20.0);
}

static inline double beta_h(double V) {
    return 1.0/(1.0+exp(-(V+35.0)/10.0));
}

static inline double h_inf(double V) {
    return alpha_h(V)/(alpha_h(V)+beta_h(V));
}

static inline double tau_h(double V) {
    return 1.0/(alpha_h(V)+beta_h(V));
}

static inline double alpha_n(double V) {
    double x = -(V+55.0);
    double y = 10.0;

    if(fabs(x/y) < 1e-6)
        return 0.01*y*(1-x/y/2.0);
    else
        return 0.01*x/(exp(x/y)-1.0);
}

static inline double beta_n(double V) {
    return 0.125*exp(-(V+65.0)/80.0);
}

static inline double n_inf(double V) {
    return alpha_n(V)/(alpha_n(V)+beta_n(V));
}

static inline double tau_n(double V) {
    return 1.0/(alpha_n(V)+beta_n(V));
}

extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new Conductance();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "I",
        "A",
        DefaultGUIModel::OUTPUT,
    },
    {
        "Vm",
        "V",
        DefaultGUIModel::INPUT,
    },
    {
        "Cm",
        "pF",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "G_Na_max",
        "nS/pF",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "E_Na",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "G_K_max",
        "nS/pF",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "E_K",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "G_L",
        "nS/pF",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "E_L",
        "mV",
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
    {
        "m",
        "Sodium Activation",
        DefaultGUIModel::STATE,
    },
    {
        "h",
        "Sodium Inactivation",
        DefaultGUIModel::STATE,
    },
    {
        "n",
        "Potassium Activation",
        DefaultGUIModel::STATE,
    },
};

static size_t num_vars = sizeof(vars)/sizeof(DefaultGUIModel::variable_t);

/*
 * Macros for making the code below a little bit cleaner.
 */

#define V (input(0)/1e3) // mV
#define m (y[0])
#define h (y[1])
#define n (y[2])
#define dm (dydt[0])
#define dh (dydt[1])
#define dn (dydt[2])
#define G_Na (G_Na_max*m*m*m*h)
#define G_K  (G_K_max*n*n*n*n)

Conductance::Conductance(void)
    : DefaultGUIModel("Conductance",::vars,::num_vars) {
    createGUI(vars, num_vars);
    /*
     * Initialize Parameters
     */
    Cm = 1.0;
    G_Na_max = 1.0;
    E_Na = 50.0;
    G_K_max = 0.4;
    E_K = -79.0;
    G_L = 0.005;
    E_L = -54.4;
    R_L_MC = 0.5181;
    Iapp_offset = 0.0;
    rate = 40000;

    /*
     * Initialize Variables
     */
    m = m_inf(V);
    h = h_inf(V);
    n = n_inf(V);
    period = RT::System::getInstance()->getPeriod()*1e-9;
    steps = static_cast<int>(ceil(period*rate));

    /*
     * Initialize States
     */
    setState("m",m);
    setState("h",h);
    setState("n",n);

    /*
     * Initialize GUI
     */
    setParameter("Cm",Cm);
    setParameter("G_Na_max",G_Na_max);
    setParameter("E_Na",E_Na);
    setParameter("G_K_max",G_K_max);
    setParameter("E_K",E_K);
    setParameter("G_L",G_L);
    setParameter("E_L",E_L);
    setParameter("R_L_MC",R_L_MC);
    setParameter("Iapp_offset",Iapp_offset);
    setParameter("rate",rate);
    refresh();
}

Conductance::~Conductance(void) {}

/*
 * Simple Euler solver.
 */

void Conductance::solve(double dt, double *y) {
    double dydt[MODEL_DIM];

    derivs(y,dydt);

    for(size_t i = 0;i < MODEL_DIM;++i)
        y[i] += dt*dydt[i];
}

void Conductance::derivs(double *y,double *dydt) {
    dm = (m_inf(V)-m)/tau_m(V);
    dh = (h_inf(V)-h)/tau_h(V);
    dn = (n_inf(V)-n)/tau_n(V);
}

#define I_ionic (G_Na*(V-E_Na)+G_K*(V-E_K)+G_L*(V-E_L))
void Conductance::execute(void) {

    /*
     * Because the real-time thread may run much slower than we want to
     *   integrate we need to run multiple interations of the solver.
     */

    for(int i = 0;i < steps;++i)
        solve(period/steps,y);

    // Membrane capacitance Cm is added to enable the use of nS/pF; default 1 pF   
    output(0) = (Iapp_offset - I_ionic*Cm + V/R_L_MC)*1e-12; //convert from pA to A
}

void Conductance::update(DefaultGUIModel::update_flags_t flag) {
    switch(flag) {
        case MODIFY:
            Cm = getParameter("Cm").toDouble();
            G_Na_max = getParameter("G_Na_max").toDouble();
            E_Na = getParameter("E_Na").toDouble();
            G_K_max = getParameter("G_K_max").toDouble();
            E_K = getParameter("E_K").toDouble();
            G_L = getParameter("G_L").toDouble();
            E_L = getParameter("E_L").toDouble();
            R_L_MC = getParameter("R_L_MC").toDouble();
            Iapp_offset = getParameter("Iapp_offset").toDouble();
            rate = getParameter("rate").toDouble();
            steps = static_cast<int>(ceil(period*rate));

            m = m_inf(V);
            h = h_inf(V);
            n = n_inf(V);
        break;
        case PERIOD:
            period = RT::System::getInstance()->getPeriod()*1e-9; // ms
            steps = static_cast<int>(ceil(period*rate));
        break;
        case PAUSE:
            output(0) = 0.0; 
        break;
        case UNPAUSE:
        break;
        default:
        break;
    }
}
