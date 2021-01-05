#ifndef EVENT_STACK_UNSTUCK_H
#define EVENT_STACK_UNSTUCK_H

#include "Event_stack.h"

class Event_stack_unstuck: public Event_stack
{
    public:
        Event_stack_unstuck(long long int id_);
        virtual ~Event_stack_unstuck();

        long long int get_correspondingposition(long long int unstuckpos);
        void set_correspondingposition(long long int pos);
        void add_correspondingposition(long long int pos);


        long long int set_unstuck_events(Event_stack * mainstack, Control * control);
        long long int set_unstuck_events_complex(Event_stack * mainstack, Control * control);

        vector<long long int> get_correspondingpostion();



    protected:

    private:

        vector<long long int> correspondingposition;

};

#endif // EVENT_STACK_UNSTUCK_H
