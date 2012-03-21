/*
 * Copyright (C) 2012 University of Bristol, UK
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

// All includes are in the below header
#include <IHCs_model_ICs_cont.h>

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

// When defining some expression always use brackets.
#define pi M_PI
#define dcell 15.0
#define acell (pi*dcell*dcell)
#define vcell (pi/6000.0*dcell*dcell*dcell)
#define cm    (acell/1e5)
#define alpha (1e5/(2*9.65*acell))
#define beta  (acell/(1000*vcell))
#define betaer (acell/(100*vcell))
/*
 *  Plugin body
 */

extern "C" Plugin::Object *createRTXIPlugin(void) {
    return new Conductance();
}

static DefaultGUIModel::variable_t vars[] = {
    {
        "i",
        "A",
        DefaultGUIModel::OUTPUT,
    },
    {
        "v",
        "V",
        DefaultGUIModel::INPUT,
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
    // This is still a parameter, we just not allow to change it from GUI.
    // Also, it can be monitored using Oscilloscope as a state.
    // Form a plugin point of view, it is still a variable.
    {
        "gKCa",
        "nS",
        DefaultGUIModel::STATE,
    },
    {
        "gCaL",
        "nS",
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
        "nS",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "tauSDR",
        "ms",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "gLeak",
        "nS",    
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
        "rLeakMC",
        "GOhms",
        DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,
    },
    {
        "iAppOffset",
        "pA - Current added to the input.",
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

#define v (input(0)*1e3)
#define c (y[0])
#define nDR (y[1])
#define sDR (y[2])
#define dc (dydt[0])
#define dnDR (dydt[1])
#define dsDR (dydt[2])

// Our non-RT pthread function, i.e., what it does
// Function can be called anything, we call it change_parameter
static void *change_parameter(void *arg)
{
    double *par = (double *)arg;
    long long i = 0;

    while(1) // Infinite loop
    {
        *par = 3.0*sin(pi*i/10); // Some sin value
        i++;
        usleep(1e5); // sleep 100 ms
    }

    return 0;
}

Conductance::Conductance(void)
    : DefaultGUIModel("IHCs model ICs",::vars,::num_vars) {
    createGUI(vars, num_vars);
    /*
     * Initialize Parameters
     */
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
    rLeakMC = 0.5181;
    //
    iAppOffset = 0.0;
    rate = 1e4;    

    /*
     * Initialize Variables
     */
    c = c0;
    nDR = nDRInf(v);
    sDR = sDRInf(v);
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
    setParameter("c0", c0);
    setParameter("pER", pER);
    setParameter("cER", cER);
    setState("gKCa", gKCa); // Now, a state; updating just GUI
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
    setParameter("rLeakMC", rLeakMC);
    //
    setParameter("iAppOffset", iAppOffset);
    setParameter("rate", rate);

    refresh();

    // Start our thread
    int retval = pthread_create(&thread,NULL,&::change_parameter,&gKCa);
    if(retval)
        ERROR_MSG("RT::OS::createTask : pthread_create failed\n");
}

Conductance::~Conductance(void) 
{
    pthread_join(thread, 0);
}

/*
 * Simple Euler solver.
 */

void Conductance::solve(double dt, double *y) {
    double dydt[MODEL_DIM];

    derivs(y,dydt);

    for(size_t i = 0;i < MODEL_DIM;++i)    
        y[i] += dt*dydt[i];
}

// Ionic currents and conductances
#define iCaL  (gCaL*sl(c,ksl)*(mlInf(v))*(mlInf(v))*(v-eCa))
#define iKCa  (gKCa*(c*c*c*c)/((c*c*c*c)+(kmKCa*kmKCa*kmKCa*kmKCa))*(v-eK))
#define iKDR  (gKDR*nDR*sDR*(v-eK))
#define iLeak (gLeak*(v-eLeak))
void Conductance::derivs(double *y,double *dydt) {
    // RHS
    dc = f*beta*(-alpha*iCaL-jEff(c,nuMP) - (nuER/(f*beta))*c + (pER/(f*beta))*(cER-c));
    dnDR = (nDRInf(v)-nDR)/tauNDR(v);
    dsDR = (sDRInf(v)-sDR)/tauSDR;
}

void Conductance::execute(void) {

    /*
     * Because the real-time thread may run much slower than we want to
     *   integrate we need to run multiple interations of the solver.
     */

    for(int i = 0;i < steps;++i)
        solve(period/steps,y);

    // Convert pA to A
    output(0) = (iAppOffset-(iCaL+iKDR+iKCa+iLeak) + 1/rLeakMC*v)*1e-12;  
}

void Conductance::update(DefaultGUIModel::update_flags_t flag) {
    switch(flag) {
        case MODIFY:
            c0 = getParameter("c0").toDouble();
            pER = getParameter("pER").toDouble();
            cER = getParameter("cER").toDouble();
            //gKCa = getParameter("gKCa").toDouble();
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
            rLeakMC = getParameter("rLeakMC").toDouble();
            //
            iAppOffset = getParameter("iAppOffset").toDouble();
            rate = getParameter("rate").toDouble();
            steps = static_cast<int>(ceil(period*rate));

            c = c0;
            nDR = nDRInf(v);
            sDR = sDRInf(v);
        break;
        case PERIOD:
            period = RT::System::getInstance()->getPeriod()*1e-9; // s
            steps = static_cast<int>(ceil(period*rate));
        break;
        case PAUSE:
            output(0) = 0.0;
        break;
        default:
        // Something unexpected happened
        break;
    }
}
