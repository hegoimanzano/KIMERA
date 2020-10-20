/*----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 *                                                              V2.0 22-03-2019
 *
 *                            'K I M E R A'
 *
 *                              C + + 11
 *
 *                       By Pablo Martín García
 *
 *
 * ----------------------------------------------------------------------
 *  Copyright (C) 2019, Pablo Martín García
 *
 *  Contact Addresses: Pablo Martín García          pablo.martinga@gmail.com
 *
 * KIMERA  is free software;  you  can  redistribute  it and/or
 * modify  it  under  the  terms  of  the  GNU General Public License  as
 * published by the Free Software Foundation
 *
 *---------------------------------------------------------------------------*/


#include "Simulation.h"



using namespace std;

int main(int argc, char *argv[])
{
Simulation sim1=Simulation();

if (argc==2)
{
    sim1.complete_simulation(argv[1]);
}
else
{
    sim1.complete_simulation("kimera.input");
}

return 0;

}


