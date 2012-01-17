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

static inline double sl(double c, double ksl) {
    return 1/(1+(c/ksl));
}

#define kMMP 0.08
static inline double jEff(double c, double nuMP) {
    return nuMP*c*c/(c*c+kMMP*kMMP);
}

#define pi M_PI
#define dcell 15.0
#define acell pi*dcell*dcell
#define vcell pi/6000.0*dcell*dcell*dcell
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
        "v",
        "V",
        DefaultGUIModel::OUTPUT,
    },
    {
        "v0",
        "mV",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "c0",
        "[Ca2+]",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "pER",
        "Unit",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "cER",
        "[Ca2+]",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gKCa",
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
        "iAppOffset",
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

#define v (y[0])
#define c (y[1])
#define nDR (y[2])
#define sDR (y[3])
#define dv (dydt[0])
#define dc (dydt[1])
#define dnDR (dydt[2])
#define dsDR (dydt[3])
#define iApp (input(0)*1e6+iAppOffset)

Cell::Cell(void)
    : DefaultGUIModel("IHCs model",::vars,::num_vars) {
    createGUI(vars, num_vars);
    /*
     * Initialize Parameters
     */
    v0 = -52.25;
    c0 = 0.335;
    pER = 0.0003;
    cER = 500.0;
    gKCa = 3.0;
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
    //
    iAppOffset = 0.0;
    rate = 1e5;    

    /*
     * Initialize Variables
     */
    v = v0;
    c = c0;
    nDR = nDRInf(v0);
    sDR = sDRInf(v0);
    period = RT::System::getInstance()->getPeriod()*1e-9;
    steps = static_cast<int>(ceil(period*rate));

    /*         
     * Initialize States
     */
    setState("c",c);
    setState("nDR",nDR);
    setState("sDR",sDR);

    /*
     * Initialize GUI
     */
    setParameter("v0", v0);
    setParameter("c0", c0);
    setParameter("pER", pER);
    setParameter("cER", cER);
    setParameter("gKCa", gKCa);
    setParameter("gCaL", gCaL);
    setParameter("nuER", nuER);
    setParameter("nuMP", nuMP);
    setParameter("kmKCa", kmKCa);
    setParameter("ksl", ksl);
    setParameter("gKDR", gKDR);
    setParameter("tauSDR", tauSDR);
    setParameter("gLeak", gLeak);
    setParameter("eLeak", eLeak);
    setParameter("f", f);
    setParameter("eCa", eCa);
    setParameter("eK", eK);
    //
    setParameter("iAppOffset", iAppOffset);
    setParameter("rate", rate);

    refresh();
}

Cell::~Cell(void) {}

/*
 * Simple Euler solver.
 */

void Cell::solve(double dt, double *y) {
    double dydt[4];

    derivs(y,dydt);

    for(size_t i = 0;i < 4;++i)    
        y[i] += dt*dydt[i];
}

void Cell::derivs(double *y,double *dydt) {
    // Ionic currents and conductances
    double iCaL, iKCa, iKDR, iLeak;
    iCaL  = gCaL*sl(c,ksl)*(mlInf(v))*(mlInf(v))*(v-eCa);
    iKCa  = gKCa*(c*c*c*c)/((c*c*c*c)+(kmKCa*kmKCa*kmKCa*kmKCa))*(v-eK);
    iKDR  = gKDR*nDR*sDR*(v-eK);
    iLeak = gLeak*(v-eLeak);
    // RHS
    dv = 1/cm*(iApp-iCaL-iKDR-iKCa-iLeak);
    dc = f*beta*(-alpha*iCaL-jEff(c,nuMP) - (nuER/(f*beta))*c + (pER/(f*beta))*(cER-c));
    dnDR = (nDRInf(v)-nDR)/tauNDR(v);
    dsDR = (sDRInf(v)-sDR)/tauSDR;
}

void Cell::execute(void) {

    /*
     * Because the real-time thread may run much slower than we want to
     *   integrate we need to run multiple interations of the solver.
     */

    for(int i = 0;i < steps;++i)
        solve(period/steps,y);

    output(0) = v/1000; //convert to mV
}

void Cell::update(DefaultGUIModel::update_flags_t flag) {
    if(flag == MODIFY) {
        v0 = getParameter("v0").toDouble();
        c0 = getParameter("c0").toDouble();
        pER = getParameter("pER").toDouble();
        cER = getParameter("cER").toDouble();
        gKCa = getParameter("gKCa").toDouble();
        gCaL = getParameter("gCaL").toDouble();
        nuER = getParameter("nuER").toDouble();
        nuMP = getParameter("nuMP").toDouble();
        kmKCa = getParameter("kmKCa").toDouble();
        ksl = getParameter("ksl").toDouble();
        gKDR = getParameter("gKDR").toDouble();
        tauSDR = getParameter("tauSDR").toDouble();
        gLeak = getParameter("gLeak").toDouble();
        eLeak = getParameter("eLeak").toDouble();
        f = getParameter("f").toDouble();
        eCa = getParameter("eCa").toDouble();
        eK = getParameter("eK").toDouble();
        //
        iAppOffset = getParameter("iAppOffset").toDouble();
        rate = getParameter("rate").toDouble();
        steps = static_cast<int>(ceil(period*rate));

        v = v0;
        c = c0;
        nDR = nDRInf(v0);
        sDR = sDRInf(v0);
    } else if(flag == PERIOD) {
        period = RT::System::getInstance()->getPeriod()*1e-9; // s
        steps = static_cast<int>(ceil(period*rate));
    }
}
