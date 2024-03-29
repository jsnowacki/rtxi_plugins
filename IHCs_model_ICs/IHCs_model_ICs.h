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

#include <default_gui_model.h>

#define MODEL_DIM 3

class Conductance : public DefaultGUIModel
{

public:

    Conductance(void);
    virtual ~Conductance(void);

    void execute(void);

protected:

    void update(DefaultGUIModel::update_flags_t);

private:

    void derivs(double *,double *);
    void solve(double,double *);

    double y[MODEL_DIM];
    double period;
    int steps;

    double c0;
    double pER;
    double cER;
    double gKCa;
    double gCaL;
    double nuER;
    double nuMP;
    double kmKCa;
    double ksl;
    double gKDR;
    double tauSDR;
    double gLeak;
    double eLeak;
    double f;
    double eCa;
    double eK;
    double rLeakMC;

    double iAppOffset;
    double rate;

};
