#ifndef SIMULATION_H
#define SIMULATION_H

#define IDUNSTUCKSTACK 2

#include "Atom.h"
#include "Box.h"
#include "Event.h"
#include "Event_stack.h"
#include "Event_stack_unstuck.h"
#include "Util.h"
#include "Control.h"
#include "Models.h"
#include "RandomGenerator.h"
#include "Event_definition.h"
#include "Tracker.h"
#include "Cell.h"
#include "Data_printer.h"
#include "Data_reader.h"
#include "Position.h"
#include "Record.h"
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <chrono>
#include <thread>



using namespace std;

class Simulation
{
    public:
        Simulation();
        virtual ~Simulation();





        void update_surface_and_events(Atom * dissolved_atom,Box *box, Event_stack * eventstack, Control *control);
        long long int dissolve_atom_by_position_event( long long int position, Box * box, Event_stack * eventstack, Control * control);

        long long int KMC_step(Event_stack * mainstack, Control *control, Box * box);
        long long int KMC_step_resolution(Event_stack * mainstack, Control *control, Box * box);


        long long int dissolve_atom_by_position_event_complex( long long int position, Box * box, Event_stack * eventstack, Control * control, Tracker * tracker);
        void update_surface_and_events_complex(Atom * dissolved_atom,Box *box, Event_stack * eventstack, Control *control);
        long long int KMC_step_complex(Event_stack * mainstack, Control *control, Box * box, Tracker * tracker);
        long long int KMC_step_complex(Event_stack * mainstack, Control * control, Box * box, Tracker * tracker,long long int seed);

        long long int complete_simulation(string name_input_file);

    protected:

    private:

};

#endif // SIMULATION_H
