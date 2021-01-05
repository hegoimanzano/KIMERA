#include "Event_stack_unstuck.h"

Event_stack_unstuck::Event_stack_unstuck(long long int id_):Event_stack(id_)
{

}


Event_stack_unstuck::~Event_stack_unstuck()
{

    for (long long int i=0; i<get_stack_size();i++)
    {
        Event *todelete= get_stack_events().at(i);
        delete todelete;  //deberia liberar la memoria guardada para todos los eventos del stack
    }

}


vector<long long int> Event_stack_unstuck::get_correspondingpostion()
{
    return correspondingposition;
}

long long int Event_stack_unstuck::get_correspondingposition(long long int unstuckpos)
{
    return correspondingposition.at(unstuckpos);
}

void Event_stack_unstuck::add_correspondingposition(long long int pos)
{
    correspondingposition.push_back(pos);
}

long long int Event_stack_unstuck::set_unstuck_events(Event_stack * mainstack, Control * control) //seria un stack2 espontanio que se eliminaria rapidamente
{
    double compare=control->get_unstuck_compare_time();
    double rate_;
    long long int eventid;
    long long int limit = mainstack->get_stack_size();

    //ya si eso..
    for (long long int i=0;i<limit;i++)
    {
        rate_=mainstack->get_eventrate(i);
        if ((1.0/rate_)>compare)
        {
            eventid=increment_by1_totaleventsid();
            Event * evento= new Event(eventid,DISSOLUTION_TYPE,rate_,mainstack->get_involved_first_atom(i));  //!!!WATCH OUT
            add_event(evento);
            add_correspondingposition(i);
        }
    }
    //set_totaleventsid(limit+1);
    return limit; //number of initial events;
}

long long int Event_stack_unstuck::set_unstuck_events_complex(Event_stack * mainstack, Control * control) //it is a temporal stack fastly removed
{
    double compare=control->get_unstuck_compare_time();
    double rate_;
    long long int eventid=0;
    long long int limit = mainstack->get_stack_size();
    Event * event1;
    #pragma omp parallel for
    for (long long int i=0;i<limit;i++)
    {
        rate_=mainstack->get_eventrate(i);
        if ((1.0/rate_)>compare)
        {
            #pragma omp critical
            {
                event1=mainstack->get_event_from_stack(i);

                Event * evento= new Event(eventid,DISSOLUTION_TYPE,rate_,mainstack->get_involved_first_atom(i));  //!!!WATCH OUT
                eventid=increment_by1_totaleventsid();
                evento->set_involved_atom_type(event1->get_involved_atom_type());
                evento->set_distance_neighbours(event1->get_distance_neighbours());
                evento->set_type_neighbours(event1->get_type_neighbours());
                evento->set_list_linked_neighbour(event1->get_list_linked_neighbour());
                evento->set_Ed_neighbours(event1->get_Ed_neighbours());
                evento->set_Ep_neighbours(event1->get_Ep_neighbours());
                evento->set_Ed_neighbours_direct(event1->get_Ed_neighbours_direct());
                evento->set_Ep_neighbours_direct(event1->get_Ep_neighbours_direct());
                evento->set_deltaG(event1->get_deltaG());
                evento->set_ffd(event1->get_ffd());
                evento->set_ffp(event1->get_ffp()); //habia error

                add_event(evento);
                add_correspondingposition(i);
            }
        }
    }
    //set_totaleventsid(limit+1); // esto no lo pillo creo que es sin mas.. para poner un numero
    return limit; //number of initial events;
}
