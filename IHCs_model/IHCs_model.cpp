/*
 * Copyright (C) 2006 University of Bristol, UK
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
#include <IHCs_model.h>

/*
 * Model Functions
 */
#define vml -26.7
#define kml 11.5
static inline double mlInf(double v) {
    return 1/(1+exp(-(v-vml)/kml));
}

#define vNDR -16.0
#define kNDR 10.0
static inline double nDRInf(double v) {
    return 1/(1+exp(-(v-vNDR)/kNDR));
}

static inline double tauNDR(double v) {
    return 0.0022+0.0029*exp(-v/14.3);
}

#define vSDR1 -60.5
#define kSDR1 6.8
#define vSDR2 -17.8
#define kSDR2 7.1
static inline double sDRInf(double v) {
    return 0.214+0.355/(1+exp((v-vSDR1)/kSDR1))+0.448/(1+exp((v-vSDR2)/kSDR2));
}

static inline double sl(double c) {
    return 1/(1+(ca/ksl));
}

#define kMMP 0.08
static inline double jEff(double c) {
    return nuMP*c^2/(c^2+kMMP^2);
}

#define pi M_PI
#define acell pi*dcell^2
#define vcell pi/6000.0*dcell^3
#define cm    acell/1e5
#define alpha 1e5/(2*9.65*acell)
#define beta  acell/(1000*vcell)
#define betaer acell/(100*vcell)

/*
 *  Plugin body
 */

extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new Cell();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "Iapp",
        "A",
        DefaultGUIModel::INPUT,
    },
    {
        "Vm",
        "V",
        DefaultGUIModel::OUTPUT,
    },
    {
        "pER",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "caER",
        "[Ca2+]",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gKCaMax",
        "mS/cm^2",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gCaL",
        "mS/cm^2",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "nuER",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "nuMP",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "kmKCa",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "ksl",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gKDR",
        "mS/cm^2",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "tauSDR",
        "ms",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gLeak",
        "mS/cm^2",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "eLeak",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "f",
        "Flux",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "eCa",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "eK",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "Iapp_offset",
        "uA/cm^2 - Current added to the input.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "rate",
        "Hz - The rate of integration.",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER,
    },
    {
        "c",
        "Ca2+ concentration",
        DefaultGUIModel::STATE,
    },
    {
        "nDR",
        "Potassium activation",
        DefaultGUIModel::STATE,
    },
    {
        "sDR",
        "Potassium activation",
        DefaultGUIModel::STATE,
    },
};

static size_t num_vars = sizeof(vars)/sizeof(DefaultGUIModel::variable_t);

/*
 * Macros for making the code below a little bit cleaner.
 */

#define V (y[0])
#define c (y[1])
#define nDR (y[2])
#define sDR (y[3])
#define dV (dydt[0])
#define dc (dydt[1])
#define dnDR (dydt[2])
#define dsDR (dydt[3])
#define Iapp (input(0)+iAppOffset)

Neuron::Neuron(void)
    : DefaultGUIModel("Cell",::vars,::num_vars) {
    createGUI(vars, num_vars);
    /*
     * Initialize Parameters
     */
    pER = 0.0003;
    caER = 500.0;
    gKCaMax 3.0;
    gCaL = 2.4;
    nuER = 1.2;
    nuMP = 3.6;
    kmKCa = 1.25;
    ksl = 0.6;
    gKDR = 2.85;
    tauSDR = 0.55;
    gLeak = 0.12;
    eLeak = -20.0;
    f = 0.01;
    eCa = 60.0;
    eK = -60.0;
    dCell = 15.0;
    //
    iAppOffset = 0.0;
    rate = 40000;

    /*
     * Initialize Variables
     */
    V = V0;
    m = m_inf(V0);
    h = h_inf(V0);
    n = n_inf(V0);
    period = RT::System::getInstance()->getPeriod()*1e-6;
    steps = static_cast<int>(ceil(period/rate/1000.0));

    /*
     * Initialize States
     */
    setState("m",m);
    setState("h",h);
    setState("n",n);

    /*
     * Initialize GUI
     */
    setParameter("V0",V0);
    setParameter("Cm",Cm);
    setParameter("G_Na_max",G_Na_max);
    setParameter("E_Na",E_Na);
    setParameter("G_K_max",G_K_max);
    setParameter("E_K",E_K);
    setParameter("G_L",G_L);
    setParameter("E_L",E_L);
    setParameter("Iapp_offset",Iapp_offset);
    setParameter("rate",rate);

    refresh();
}

Neuron::~Neuron(void) {}

/*
 * Simple Euler solver.
 */

void Neuron::solve(double dt, double *y) {
    double dydt[4];

    derivs(y,dydt);

    for(size_t i = 0;i < 4;++i)
        y[i] += dt*dydt[i];
}

void Neuron::derivs(double *y,double *dydt) {
    dV = (Iapp - G_Na*(V-E_Na) - G_K*(V-E_K) - G_L*(V-E_L)) / Cm;
    dm = (m_inf(V)-m)/tau_m(V);
    dh = (h_inf(V)-h)/tau_h(V);
    dn = (n_inf(V)-n)/tau_n(V);
}

void Neuron::execute(void) {

    /*
     * Because the real-time thread may run much slower than we want to
     *   integrate we need to run multiple interations of the solver.
     */

    for(int i = 0;i < steps;++i)
        solve(period/steps,y);

    output(0) = V/1000; //convert to mV
}

void Neuron::update(DefaultGUIModel::update_flags_t flag) {
    if(flag == MODIFY) {
        V0 = getParameter("V0").toDouble();
        Cm = getParameter("Cm").toDouble();
        G_Na_max = getParameter("G_Na_max").toDouble();
        E_Na = getParameter("E_Na").toDouble();
        G_K_max = getParameter("G_K_max").toDouble();
        E_K = getParameter("E_K").toDouble();
        G_L = getParameter("G_L").toDouble();
        E_L = getParameter("E_L").toDouble();
        Iapp_offset = getParameter("Iapp_offset").toDouble();
        rate = getParameter("rate").toDouble();
        steps = static_cast<int>(ceil(period*rate/1000.0));

        V = V0;
        m = m_inf(V0);
        h = h_inf(V0);
        n = n_inf(V0);
    } else if(flag == PERIOD) {
        period = RT::System::getInstance()->getPeriod()*1e-6; // ms
        steps = static_cast<int>(ceil(period*rate/1000.0));
    }
}
