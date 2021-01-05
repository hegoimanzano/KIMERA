#ifndef EVENT_STACK_H
#define EVENT_STACK_H

#include "Atom.h"
#include "Box.h"
#include "Event.h"
#include "Control.h"
#include "Models.h"
#include "Util.h"
#include "Record.h"
#include "Event_definition.h"
#include <vector>
#include <stdexcept>

using namespace std;


class Event_stack
{
    public:
        Event_stack(long long int id_);
        virtual ~Event_stack();

        vector<Event*> get_stack_events();
        Event* get_event_from_stack( long long int pos_);
        long long int get_stack_size();

        void add_event(Event *event);
        void add_event(Atom *atom, Control *control);


        double get_eventrate( long long int position);

        long long int rm_event_by_id(long long int eventid);
        long long int rm_event_by_position( long long int position);

        long long int increment_by1_totaleventsid();

        long long int get_totaleventsid();
        void set_totaleventsid(long long int totaleventsid_);

        long long int get_stack_id();
        double get_characteristicrate(long long int position);

        long long int update_stackrate_propensity(Control * control);

        double get_stackpropensity();
        void set_stackpropensity();

        long long int get_pos_by_atomid(long long int atomid);

        long long int set_initial_events(Box * box, Control * control);

        Atom * get_involved_first_atom( long long int pos);


        long long int event_binary_search(double searchvalue); // return the possition of the event in the vector
        long long int event_lineal_search(double searchvalue); // return the possition of the event in the vector

        long long int set_initial_events_complex(Box * box, Control * control);

        long long int update_stackrate_propensity_complex();

        long long int add_event_complex(Atom *atom, Control *control);

    protected:

    private:

    long long int totaleventsid;
    long long int id;
    double characteristicrate;
    double stackpropensity;
    vector<Event*> stack_events;
};

#endif // EVENT_STACK_H
