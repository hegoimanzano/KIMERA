#include "Event_stack.h"

Event_stack::Event_stack(long long int id_)
{
    id=id_;
    totaleventsid=0;
}

Event_stack::~Event_stack() //No lo borra bien, no se porque.. se me queda 1 sin borar....
{
/**
    long long int stack_size=stack_events.size();
    for (long long int i=0; i<stack_size-1;i++)
    {
        Event *todelete= stack_events.at(i);
        delete todelete;  //deberia liberar la memoria guardada para todos los eventos del stack
    }
*/
}


vector<Event*> Event_stack::get_stack_events()
{
    return stack_events;
}

long long int Event_stack::get_stack_size()
{
    return stack_events.size();
}

void Event_stack::add_event(Event *event)
{
    stack_events.push_back(event);
}

Event* Event_stack::get_event_from_stack( long long int pos_)
{
    return stack_events.at(pos_);
}

void Event_stack::add_event(Atom *atom, Control *control)//mejor no usar
{
        long long int neighbour= atom->get_size_neighbour();
        double rate_=Models::model42_dissolution(control->get_fd(),neighbour,control->get_Ed());
        long long int eventid=increment_by1_totaleventsid();
        Event * evento= new Event(eventid,DISSOLUTION_TYPE,rate_,atom);  //!!!WATCH OUT
        add_event(evento);
}

long long int Event_stack::rm_event_by_id(long long int eventid)
{
     long long int value;
    volatile bool flag=false;

    #pragma omp parallel for shared(flag)
    for (long long int i=0; i<get_stack_size();i++)
    {
        if(flag) continue;
        if (eventid==stack_events.at(i)->get_id())
        {
            value=i;
            flag=true;
        }
    }

    if (flag)
    {
        Event *todelete= stack_events.at(value);
        stack_events.erase(stack_events.begin()+value);
        delete todelete;  //free space saved for the event

        return 0;
    }
    return 1; //fail
}

long long int Event_stack::rm_event_by_position( long long int position)
{
    Event *todelete= stack_events.at(position);
    stack_events.erase(stack_events.begin()+position);
    delete todelete;  //deberia liberar la memoria guardada para el evento
    return 0;
}

long long int Event_stack::increment_by1_totaleventsid()
{
    totaleventsid++;
    return totaleventsid;
}

long long int Event_stack::get_totaleventsid()
{
    return totaleventsid;
}

void Event_stack::set_totaleventsid(long long int totaleventsid_)
{
    totaleventsid=totaleventsid_;
}

long long int Event_stack::get_stack_id()
{
    return id;
}

double Event_stack::get_eventrate( long long int position)
{
    return get_event_from_stack(position)->get_rate();
}

double Event_stack::get_stackpropensity()
{
    return stackpropensity;
}

long long int Event_stack::update_stackrate_propensity(Control * control)
{
    long long int neighbour=-1;
    double rate_=0.0;
    double propensity=0.0;

    long long int limit=get_stack_size();
        //Ya se hara si eso
        for (long long int i=0;i<limit;i++)
        {
            neighbour = stack_events.at(i)->get_involved_atoms().at(0)->get_size_neighbour();
            rate_= Models::model42_dissolution(control->get_fd(),neighbour,control->get_Ed());
            stack_events.at(i)->set_rate(rate_);
            propensity+=rate_;
            stack_events.at(i)->set_propensity(propensity);
        }

    stackpropensity=propensity;

    return 0;
}

void Event_stack::set_stackpropensity()
{
    long long int limit=get_stack_size();

    double propensity=0.0;

        for (long long int i=0;i<limit;i++)
        {
            propensity+=stack_events.at(i)->get_rate();
            stack_events.at(i)->set_propensity(propensity);
        }
    stackpropensity=propensity;
}

long long int Event_stack::get_pos_by_atomid(long long int atomid)
{
     long long int value;
    volatile bool flag=false;

    long long int limit=get_stack_size();

    #pragma omp parallel for shared(flag)
    for (long long int i=0;i<limit;i++)
    {
        if (flag) continue;
        if(atomid==stack_events.at(i)->get_involved_atoms().at(0)->get_id())
        {
            value=i;
        }
    }
    if (flag) return value;
    return -1; //failed
}

long long int Event_stack::set_initial_events(Box * box, Control * control) //tiene que setearse después de que hayamos quitado los vacios
{
    double propensity=0.0;
    Atom * involvedatom;
    double rate_;
    long long int neighbour;
    long long int eventid=0;
    long long int limit = box->get_surface_size();
    //ya ser hara si eso
    for (long long int i=0;i<limit;i++)
    {
        involvedatom=box->get_surface_atoms().at(i);
        neighbour= involvedatom->get_size_neighbour();
        rate_=Models::model42_dissolution(control->get_fd(),neighbour,control->get_Ed());
        //cout<< "  rate_  " << rate_ <<endl;
        propensity+=rate_;
        Event * evento= new Event(eventid,DISSOLUTION_TYPE,rate_,involvedatom);  //!!!WATCH OUT
        evento->set_propensity(propensity);
        add_event(evento);
        eventid=increment_by1_totaleventsid();
    }
    stackpropensity=propensity;
    totaleventsid=limit+1;
    return totaleventsid; //number of initial events;
}


Atom * Event_stack::get_involved_first_atom( long long int pos)
{
    return stack_events.at(pos)->get_involved_atoms().at(0);
}


long long int Event_stack::event_binary_search(double searchvalue)
{
    double compara=0.0;
    long long int topearriba=get_stack_size()-1;
    if (get_stack_size()==0) {cout<< "There are no events, simulation must be aborted. "<<endl; throw exception();}

    long long int topeabajo=0;
    long long int enelbucle=0;

    long long int l=topearriba/2;

    bool encontrado=true;

    while (encontrado)
    {
        enelbucle++;
        if (enelbucle>get_stack_size())
        {
            cout<< "Watch out! in binary search loop for more than the size of the event stack " << endl;
            cout<< "value of l: "<< l << endl;
            cout<< "value of uplimit: "<< topearriba << endl;
            cout<< "value of botlimit: "<< topeabajo << endl;
            cout<< "value of propensity: "<< stack_events.at(topearriba)->get_propensity() << endl;
            cout<< "value of selected event propensity : "<< stack_events.at(l)->get_propensity() << endl;
            if(enelbucle>4*get_stack_size())
            {
                cout<< "Watch out! in binary search loop for more than four times the size of the event stack , aborting simulation" << endl;
                exit(1);
            }
        }
        if (l==0) {compara=0.0;}
        else compara=stack_events.at(l-1)->get_propensity();

        if (searchvalue>compara) //mayor que el anterior
        {
            if (searchvalue<=stack_events.at(l)->get_propensity())
            {
                encontrado=!encontrado;
            }
            else
            {
                topeabajo=l;
                l=l+ceil((topearriba-topeabajo)/2.0);
            }
        }

        if (searchvalue<compara )
        {
            topearriba=l;
            l=l-ceil((topearriba-topeabajo)/2.0);
        }
    }
    return l;
}

long long int Event_stack::event_lineal_search(double searchvalue)
{
     long long int value;
    volatile bool flag=false;

    if (get_stack_size()==0) {cout<< "There are no events, simulation must be aborted. "<<endl; throw exception();}

    if (searchvalue<stack_events.at(0)->get_propensity()) return 0;

    #pragma omp parallel for shared(flag)
    for (long long int i=1;i<get_stack_size();i++)
    {
        if (flag) continue;
        if (searchvalue>stack_events.at(i-1)->get_propensity() && searchvalue<=stack_events.at(i)->get_propensity())
        {
            flag=true;
            value=i;
        }
    }

    if(flag) {return value;}
    else {cout<< "Watch out! Not found value in lineal search" << endl; throw exception();}
}


long long int Event_stack::set_initial_events_complex(Box * box, Control * control) //tiene que setearse después de que hayamos quitado los vacios
{

    vector<Event_definition*> events_definitions=control->get_list_event_definition();
    long long int events_definitions_size=events_definitions.size();
    double propensity=0.0;
    Atom * involvedatom;
    double rate_=0.0;
    long long int eventid=0;
    long long int limit = box->get_surface_size();

    vector<vector<Atom*>> linkeds;
    vector<long long int> type_linkeds;

     vector <string> types_event_definition;
     vector <double> distances_event_definition;
     long long int sizetypes_event_definition;

     vector<long long int> Nn;
     bool all0=true;
     double ffd=0.0;

     long long int sizeaffected=0;

    //puede pasar que un mismo atomo tenga varios eventos posibles, dependiendo como esten definidos los eventos Esto me puede traer problemas

    for (long long int i=0;i<limit;i++)
    {
        involvedatom=box->get_surface_atoms().at(i);
        linkeds=involvedatom->get_linked();
        type_linkeds=involvedatom->get_linked_type();
        vector<Atom*>  neighbours_involvedatom= involvedatom->get_neighbours();
        vector<double>  distances_involvedatom= involvedatom->get_distances_to_neighbours();

        for (long long int j=0; j<events_definitions_size;j++)
        {
            if (events_definitions.at(j)->get_involved_atom_type().compare(involvedatom->get_atom_type())==0)
            {

                types_event_definition=events_definitions.at(j)->get_type_neighbours();
                sizetypes_event_definition=types_event_definition.size() ;
                distances_event_definition=events_definitions.at(j)->get_distance_neighbours();

                for (long long int k =0; k<sizetypes_event_definition;k++)
                {
                    //tengo que pillar el tipo de los linkeds asociados, si lo hay.
                    long long int typelink_from_event_def=events_definitions.at(j)->get_list_linked_neighbour().at(k).get_type();
                    long long int truenumber=0;

                    truenumber=involvedatom->get_number_complex(types_event_definition.at(k),distances_event_definition.at(k),typelink_from_event_def);

                    //long long int truenumber=involvedatom->get_number(types_event_definition.at(k),distances_event_definition.at(k)); //esta pillando el doble
                    //no lo haria con el get number, lo haria a mano


                    //if (truenumber !=0) {all0=false;}
                    //a este numero habra que restarle los que no esten linkeados si es el caso
                    //si el tipo de vector<linked_neighbour>.at(j) es 0 no pasa nada, si es distinto de 0, habra que mirar los linked si estan


                    //cout<<" antes  " <<endl;
                    //cout<<" types_event_definition.at(k)  " <<types_event_definition.at(k)<<endl;
                    //cout<<" distances_event_definition.at(k)  " <<distances_event_definition.at(k)<<endl;

                    //cout<<"truenumber: "<<truenumber<<endl;

                    //cout<<" despues  " <<endl;
                    //cout<<" types_event_definition.at(k)  " <<types_event_definition.at(k)<<endl;
                    //cout<<" distances_event_definition.at(k)  " <<distances_event_definition.at(k)<<endl;

                    //cout<<"truenumber: "<<truenumber<<endl;


                    if(truenumber>0){all0=false;}
                    if(truenumber<0){cout<<"Error defining initial events, please contact the authors"<<endl; throw invalid_argument( " Error, " );}
                    Nn.push_back(truenumber);
                }

                //mejor igual no crear el evento si falta alguno de los linked no? la definicion del evento puede tener partes q no tengan linked
                // asi que lo suyo es que Nn vaya cambiando




                ffd=events_definitions.at(j)->get_ffd();
                //NUEVO
                rate_=Models::model42_dissolution(ffd,Nn, events_definitions.at(j)->get_Ed_neighbours(),events_definitions.at(j)->get_Ed_neighbours_direct());
                //rate_=Models::model42_dissolution(ffd,Nn, events_definitions.at(j)->get_Ed_neighbours());
                propensity+=rate_;

                long long int typeevento;
                if (all0) {typeevento=WITH_OUT_NEIGH_TYPE;} //atom without any neighbour
                if (!all0) {typeevento=DISSOLUTION_TYPE;}

                all0=true;

                Event * evento= new Event(eventid,typeevento,rate_,involvedatom); //!!!WATCH OUT

                eventid=increment_by1_totaleventsid();
                sizeaffected=involvedatom->get_affected().size();

                for (long long int l=0;l<sizeaffected;l++)
                {
                    evento->add_involved_atom(involvedatom->get_affected().at(l));
                }
                evento->set_involved_atom_type(events_definitions.at(j)->get_involved_atom_type());
                evento->set_distance_neighbours(events_definitions.at(j)->get_distance_neighbours());
                evento->set_type_neighbours(events_definitions.at(j)->get_type_neighbours());
                evento->set_list_linked_neighbour(events_definitions.at(j)->get_list_linked_neighbour());
                evento->set_Ed_neighbours(events_definitions.at(j)->get_Ed_neighbours());
                evento->set_Ep_neighbours(events_definitions.at(j)->get_Ep_neighbours());
                evento->set_Ed_neighbours_direct(events_definitions.at(j)->get_Ed_neighbours_direct());
                evento->set_Ep_neighbours_direct(events_definitions.at(j)->get_Ep_neighbours_direct());
                evento->set_deltaG(events_definitions.at(j)->get_deltaG());
                evento->set_ffd(events_definitions.at(j)->get_ffd());
                evento->set_ffp(events_definitions.at(j)->get_ffp());  //habia error
                evento->set_count_involved(Nn);

                evento->set_propensity(propensity);
                add_event(evento);
                Nn.clear();
            }
        }
    }
    stackpropensity=propensity;
    return eventid; //number of initial events;
}


long long int Event_stack::update_stackrate_propensity_complex()  //I remove atom after add_event_complex so the list of neighbour can change and
{                                                       //must be recheked
//habria que paralelizar

    vector<long long int> Nn;


    long long int sizetypesevent;

    Atom * involvedatom;
    Event * involvedevent;

    vector<vector<Atom*>> linkeds;
    vector<long long int> type_linkeds;

    bool all0=true;

    double rate_=0.0;
    double propensity=0.0;

    double ffd=0.0;
    double ffp=0.0;
    double dG=0.0;

    long long int limit=get_stack_size();

    for (long long int i=0;i<limit;i++)
    {
        involvedevent =  stack_events.at(i);
        involvedatom = stack_events.at(i)->get_involved_atoms().at(0);

        sizetypesevent= involvedevent->get_type_neighbours().size();

        linkeds=involvedatom->get_linked();
        type_linkeds=involvedatom->get_linked_type();

        vector<Atom*> neighbours_involvedatom= involvedatom->get_neighbours();
        vector<double> distances_involvedatom= involvedatom->get_distances_to_neighbours();

        //we have the atom, we look at the neighbor amount in each, and recalculate the rate.

        for (long long int k =0; k<sizetypesevent;k++)
        {

            long long int truenumber=involvedatom->get_number_complex(involvedevent->get_type_neighbours().at(k),involvedevent->get_distance_neighbours().at(k),involvedevent->get_list_linked_neighbour().at(k).get_type());


            if(truenumber>0){all0=false;}
            if(truenumber<0){cout<<"Error updating event rates, please contact the authors"<<endl; throw invalid_argument( " Error, " );}

            Nn.push_back(truenumber);

        }

        if (all0) involvedevent->set_type(WITH_OUT_NEIGH_TYPE);
        if (!all0) involvedevent->set_type(DISSOLUTION_TYPE);

        all0=true;

        ffd=involvedevent->get_ffd();
        ffp=involvedevent->get_ffp();
        dG=involvedevent->get_deltaG();


            //Nuevo

        rate_=Models::model50_dissolution(ffd, ffp, Nn, involvedevent->get_Ed_neighbours() , involvedevent->get_Ep_neighbours(),
                                              involvedevent->get_Ed_neighbours_direct(), involvedevent->get_Ep_neighbours_direct() ,dG);

        //Nuevo
        //rate_=Models::model42_dissolution(involvedevent->get_ffd(), Nn, involvedevent->get_Ed_neighbours(),involvedevent->get_Ed_neighbours_direct());
        //rate_=Models::model42_dissolution(involvedevent->get_ffd(), Nn, involvedevent->get_Ed_neighbours());
        involvedevent->set_rate(rate_);
        propensity+=rate_;
        involvedevent->set_propensity(propensity);
        involvedevent->set_count_involved(Nn);      //hemos actualizado tmb la cantidad de  vecinos
        Nn.clear();
     }

    stackpropensity=propensity;

    return 0;
}


long long int Event_stack::add_event_complex(Atom *atom, Control *control)  //cambiar esto
{

    vector<Event_definition*> events_definitions=control->get_list_event_definition();
    long long int events_definitions_size=events_definitions.size();
    double rate_=0.0;
    long long int eventid=increment_by1_totaleventsid();
    vector<Record> atomrecord;

     vector <string> types_event_definition;
     vector <double> distances_event_definition;
     long long int sizetypes_event_definition;

     vector<long long int> Nn;
     bool all0=true;
     double ffd=0.0;
     double ffp=0.0;
     double dG=0.0;

     long long int sizeaffected=0;

     double propensity=0.0;

    vector<vector<Atom*>> linkeds;
    vector<long long int> type_linkeds;



    //puede pasar que un mismo atomo tenga varios eventos posibles, dependiendo como esten definidos los eventos esto nos puede dar problemas

    for (long long int j=0; j<events_definitions_size;j++)
    {
        if (events_definitions.at(j)->get_involved_atom_type().compare(atom->get_atom_type())==0)
        {
            types_event_definition=events_definitions.at(j)->get_type_neighbours();
            sizetypes_event_definition=types_event_definition.size() ;
            distances_event_definition=events_definitions.at(j)->get_distance_neighbours();

            linkeds=atom->get_linked();
            type_linkeds=atom->get_linked_type();

            vector<Atom*>  neighbours_involvedatom= atom->get_neighbours();
            vector<double>  distances_involvedatom= atom->get_distances_to_neighbours();


            for (long long int k =0; k<sizetypes_event_definition;k++)
            {
                long long int typelink_from_event_def=events_definitions.at(j)->get_list_linked_neighbour().at(k).get_type();
                long long int truenumber=0;

                truenumber=atom->get_number_complex(types_event_definition.at(k),distances_event_definition.at(k),typelink_from_event_def);

                if(truenumber>0){all0=false;}
                if(truenumber<0){cout<<"Error adding an event, please contact the authors"<<endl; throw invalid_argument( " Error, " );}
                Nn.push_back(truenumber);
            }


            long long int typeevento;
            if (all0) {typeevento=WITH_OUT_NEIGH_TYPE;} //atom without any neighbour
            if (!all0) {typeevento=DISSOLUTION_TYPE;}

            all0=true;

            ffd=events_definitions.at(j)->get_ffd();
            ffp=events_definitions.at(j)->get_ffp();
            dG=events_definitions.at(j)->get_deltaG();


            //Nuevo

            rate_=Models::model50_dissolution(ffd, ffp, Nn, events_definitions.at(j)->get_Ed_neighbours() , events_definitions.at(j)->get_Ep_neighbours(),
                                              events_definitions.at(j)->get_Ed_neighbours_direct(), events_definitions.at(j)->get_Ep_neighbours_direct() ,dG);
            //rate_=Models::model42_dissolution(ffd,Nn, events_definitions.at(j)->get_Ed_neighbours(),events_definitions.at(j)->get_Ed_neighbours_direct());
            //rate_=Models::model42_dissolution(ffd,Nn, events_definitions.at(j)->get_Ed_neighbours());

            propensity+=rate_;
            Event * evento= new Event(eventid,typeevento,rate_,atom);  //!!!WATCH OUT
            eventid=increment_by1_totaleventsid();
            sizeaffected=atom->get_affected().size();

            for (long long int l=0;l<sizeaffected;l++)
            {
                evento->add_involved_atom(atom->get_affected().at(l));
            }
            evento->set_involved_atom_type(events_definitions.at(j)->get_involved_atom_type());
            evento->set_distance_neighbours(events_definitions.at(j)->get_distance_neighbours());
            evento->set_type_neighbours(events_definitions.at(j)->get_type_neighbours());
            evento->set_list_linked_neighbour(events_definitions.at(j)->get_list_linked_neighbour());
            evento->set_Ed_neighbours(events_definitions.at(j)->get_Ed_neighbours());
            evento->set_Ep_neighbours(events_definitions.at(j)->get_Ep_neighbours());
            evento->set_Ed_neighbours_direct(events_definitions.at(j)->get_Ed_neighbours_direct());
            evento->set_Ep_neighbours_direct(events_definitions.at(j)->get_Ep_neighbours_direct());
            evento->set_deltaG(events_definitions.at(j)->get_deltaG());
            evento->set_ffd(events_definitions.at(j)->get_ffd());
            evento->set_ffp(events_definitions.at(j)->get_ffp()); //habia error
            evento->set_count_involved(Nn);
            Nn.clear();

            evento->set_propensity(propensity);
            add_event(evento);

        }
    }
    //stackpropensity=propensity;   // no etiendo esto.. no deberia de estar no? lo unico que hemos hecho ha sido meter un evento la propensity del evento no es la total
    return eventid; //number of initial events;

}



