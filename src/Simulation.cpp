#ifdef __linux__
#define REMOVE "rm "
#endif // __linux__

#ifdef __WIN32__
#define REMOVE "del "
#endif // _WINDOWS

#ifdef __WIN64__
#define REMOVE "del "
#endif // _WINDOWS

#ifdef __APPLE__
#define REMOVE "rm "
#endif // _APPLE


#include "Simulation.h"

Simulation::Simulation()
{

}

Simulation::~Simulation()
{
    //dtor
}

void Simulation::update_surface_and_events(Atom * dissolved_atom,Box *box, Event_stack * eventstack, Control *control)
{
    long long int size=dissolved_atom->get_size_neighbour();
    for (long long int i=0;i<size;i++)
    {
        if (!box->check_in_surface(dissolved_atom->get_id_neighbour(i)))
        {
            box->add_surface_atom(dissolved_atom->get_neighbour(i));
            eventstack->add_event(dissolved_atom->get_neighbour(i), control);
        }
    }
}

void Simulation::update_surface_and_events_complex(Atom * dissolved_atom,Box *box, Event_stack * eventstack, Control *control)
{
    //puede pasar que al disolverse un atomo, un atomo, otro atomo en la surface ya no lo sea? si, puede pasar.. pero mala suerte
    // a la hora de definir los eventos se va a mirar.


    //puede pasar que  tenga vecinos pero que sus vecinos no le tengan a el como vecino..¿en serio? creo q si
    //pues vamos a suponer que los atomos que van a entrar en surface, son aquellos que son vecinos del atomo disuelto
    //habra que mirar lo de si sigue siendo  bulk o no a pesar del atomo que se ha ido!!!!!
    //tendremos que importer de control la definicion de evento para ver si es bulk o no, y por lo tanto si le vamos a meter en superficie o no.
    //tener en cuenta tambien si es insouble
    vector<Atom*> neighbors = dissolved_atom->get_neighbours();// tenemos los vecinos                Añadir que si tenemos un vecino que esta en dislocacion, ponerle como  Withoutneighbours
    long long int sizeneighbors=neighbors.size();

    //remove the atom and modify the neoghbors
    long long int dissolved_type= dissolved_atom->get_type();

    box->rm_atom_complex(dissolved_atom->get_id());

    vector<Event_definition*> list_definition=control->get_list_event_definition();
    long long int size_list_definition=list_definition.size();

    if (dissolved_type==NORMAL)
    {
        #pragma omp parallel for
        for (long long int i=0;i<sizeneighbors;i++)
        {
            Atom * targetatom=neighbors.at(i);

            vector<Atom*>  neighbours_involvedatom= targetatom->get_neighbours();
            vector<double>  distances_involvedatom= targetatom->get_distances_to_neighbours();

            if (!targetatom->get_insurface() && targetatom->get_type()==NORMAL) //if not in surface and normal
            {
                //we check if is still a bulk atom (the neighbor we want to enter in surface)

                for (long long int j=0;j<size_list_definition;j++)
                {
                    if (list_definition.at(j)->get_involved_atom_type().compare(targetatom->get_atom_type())==0)
                    {

                        vector<Linked_neighbour> linkvector=list_definition.at(j)->get_list_linked_neighbour();

                        bool bulk=true;
                        bulk=targetatom->is_bulk_with_links(list_definition.at(j)->get_max_to_bulk(),
                                                 list_definition.at(j)->get_type_neighbours(),list_definition.at(j)->get_distance_neighbours(),
                                                 linkvector);



                        if (!bulk)
                        {
                            #pragma omp critical
                            {
                                box->add_surface_atom(targetatom);
                                //targetatom->set_insurface(true); //no need, add_surface_atom do it for you
                                eventstack->add_event_complex(targetatom, control);
                                //tiene que ser los vecinos del que acaba de entrar en surface:
                            }
                            /**
                            for (long long int s=0;s<targetatom->get_size_neighbour();s++)
                            {
                                Atom * targetneighatom= targetatom->get_neighbour(s);

                                if ( targetneighatom->get_type()==INDISLOCATION) //SI EL VECINO ESTA EN UNA DISLOCACION, HACEMOS COMO SI SE ESTUVIERA DISOLVIENDO EL MISMO Y PONEMOS SUS VECINOS EN JUEGO
                                {
                                   //vector<Atom*> neighbors2 = targetneighatom->get_neighbours();
                                   //box->rm_atom_complex( targetneighatom->get_id());
                                   update_surface_and_events_complex(targetneighatom,box,eventstack, control);  //de momento esto no esta muy claro,
                                   //digamos que si lo penemos se abriría toda la cavidad, quizas en un futuro para un caso futuro..

                                   //deberiamos añadir aqui los vecinos de este que sean normales con el add_event_complex

                                   for (long long int l=0;l<neighbors2.size();l++)
                                   {
                                       if ( neighbors2.at(l)->get_type()==NORMAL) //SI EL VECINO ESTA EN UNA DISLOCACION, HACEMOS COMO SI SE ESTUVIERA DISOLVIENDO EL MISMO Y PONEMOS SUS VECINOS EN JUEGO
                                       {








                                       }
                                   }
                                   /

                                }
                            }
                            */
                            break;
                        }
                    }
                }
            }
        }
    }
}






/**

                            linkeds=targetatom->get_linked();
                            type_linkeds=targetatom->get_linked_type();

                            if (linkeds.size()==0) {in=true;}

                            for (long long int s=0; s<linkeds.size() ;s++)
                            {
                                long long int type = type_linkeds.at(s);

                                double distance_to_neight=distances_involvedatom.at(s);
                                string type_neight=neighbours_involvedatom.at(s)->get_atom_type();

                                in=false;

                                if (type==NO_LINKED_TYPE)
                                {
                                    in=true;
                                    break;
                                }
                                if (type == LINKED_TYPE_NORMAL && list_definition.at(j)->get_type_neighbours().compare(type_neight)==0  && Util::isEqualDistances(list_definition.at(j)->get_distance_neighbours(),distance_to_neight))
                                {
                                    //comprobamos que los linked que tiene el mismo id correspondiente, estan (todos ellos)
                                    vector<Atom*> inlinkeds= linkeds.at(s);
                                    for (long long int l=0;l<inlinkeds.size() ;l++)
                                    {
                                        if  (inlinkeds.at(l)->get_type()==NORMAL)
                                        {
                                            in=true;
                                        }
                                    }
                                }
                                if (type == LINKED_TYPE_DISSOLVED && list_definition.at(j)->get_type_neighbours().compare(type_neight)==0  && Util::isEqualDistances(list_definition.at(j)->get_distance_neighbours(),distance_to_neight))
                                {
                                    //comprobamos que los linked que tiene el mismo id correspondiente, estan (todos ellos)
                                    vector<Atom*> inlinkeds= linkeds.at(s);
                                    for (long long int l=0;l<inlinkeds.size() ;l++)
                                    {
                                        if  (inlinkeds.at(l)->get_type()==DISSOLVED)
                                        {
                                            in=true;
                                        }
                                    }
                                }
                            }
                            */
//tendre que mirar los linked para ver si lo meto a la superficie o no..


long long int Simulation::dissolve_atom_by_position_event( long long int position, Box * box, Event_stack * eventstack, Control * control)
{
    //Primero metemos los atomos en surface y sus correspondientes eventos events_stack
    Atom * dissolved_atom=eventstack->get_event_from_stack(position)->get_involved_atoms().at(0);
    update_surface_and_events(dissolved_atom, box, eventstack,  control);


    //primero borramos el evento que corresponde al atomo disuelto
    eventstack->rm_event_by_position(position);

    //Y Luego borramos el atomo (los vecinos de este ven su numero de vecinos actualizado, pero no su rate)
    if(box->rm_atom(dissolved_atom->get_id())!=0) {return 1;}

    // Y por ultimo actualizamos el rate de los átomos vecinos y ya que estamos la propensity
    eventstack->update_stackrate_propensity(control);

    return 0;
}


long long int Simulation::dissolve_atom_by_position_event_complex( long long int position, Box * box, Event_stack * eventstack, Control * control, Tracker * tracker)
{
    //Primero metemos los atomos en surface y sus correspondientes eventos events_stack, y el affected tmb
    Event * chosen_event=eventstack->get_event_from_stack(position);

    long long int size_involved_atoms=chosen_event->get_involved_atoms().size();

    Atom * dissolved_atom;

    //No paralelizar
    for (long long int i=0;i<size_involved_atoms;i++)
    {
            dissolved_atom=chosen_event->get_involved_atoms().at(i);

            //Le metemos al tracker si su tipo es normal..(si aun no estaba disuelto)

            if (dissolved_atom->get_type()==NORMAL) tracker->add_atom_dissolved(dissolved_atom);

            //Aqui aun tenemos la informacion del evento

            update_surface_and_events_complex(dissolved_atom, box, eventstack,  control);  //estoy intentando meter atomos que aun consideran que el atomo no esta disuelto

            //Y Luego borramos el atomo (los vecinos de este ven su numero de vecinos actualizado, pero no su rate)  HAY QUE BAJAR TMB LA CUENTA DEL RECORD
            //box->rm_atom_complex(dissolved_atom->get_id()); sobraria pq lo hemos metido dentro de update_surface_and_events_complex
    }

    //luego borramos el evento que corresponde al atomo disuelto
    eventstack->rm_event_by_position(position);

    // Y por ultimo actualizamos el rate de los átomos vecinos

    // y ya que estamos la propensity

    eventstack->update_stackrate_propensity_complex();

    return 0;
}


long long int Simulation::KMC_step_complex(Event_stack * mainstack, Control * control, Box * box, Tracker * tracker)  //ESTE Y EL SIGUIENTE SON LOS QUE SE USAN
{
    //Si el tamaño de los eventos del stack es 0.. ya no tenemos nada que hacer..
    if (mainstack->get_stack_size()==0) return 1; //fail

    double propensity;
    //importamos el inicializador del numero aleatorio (q capullos)
    mt19937 &mt = RandomGenerator::Instance().get();
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    mt.seed(seed);
    uniform_real_distribution<double> dist01(0.0, 1.0);
    long long int mainpos;
    double kmctime;
    //Set initial events tmb!! (tendra q ser antes, al comienzo del todo)
    //Seteamos el rate y la propensity en en main stack

    //si no se disuelve en la anterior el rate y la propensity es la misma..

    //mainstack->set_stackpropensity(); //no haria falta pq ya lo hemos definido
    //cout<< " -- mainstack->get_stackpropensity();  " << mainstack->get_stackpropensity() << endl;
    //Si llegamos a que alert es tan alto como unstuckeverystep, creamos el stack unstuck para agilizar
    if (control->get_unstuckalert()==control->get_unstuckeverystep())
    {
        Event_stack_unstuck unstuckstack(IDUNSTUCKSTACK);
        //Un vez creado, creamos los eventos oportunos
        unstuckstack.set_unstuck_events_complex(mainstack, control);
        //si el tamanno del unstuckstack es 0.. la liamos
        if (unstuckstack.get_stack_size()==0) {
            cout<< "unstuckstack without elements!, try to reduce TARGET_TIME or/and increase ACCURACY"<<endl;
            control->restart_unstuckalert();
            //reentramos en el metodo
            KMC_step_complex( mainstack, control, box, tracker);

        }
        //seteamos la propensity de cada evento y la total
        unstuckstack.set_stackpropensity();
        //generamos el numero aleatorio segun la propensity del stack...
        double unstuckstackpropensity= unstuckstack.get_stackpropensity();
        propensity=unstuckstackpropensity;

        //.....
        uniform_real_distribution<double> dist(0.0, unstuckstackpropensity);

        //buscamos segun el valor en el stack unstuck
        long long int unstuckpos;

        if (control->get_linealsearch()) {unstuckpos=unstuckstack.event_lineal_search(dist(mt));}
        else {unstuckpos=unstuckstack.event_binary_search(dist(mt));}

        //miramos el correspondiente posición en el stack principal
        mainpos=unstuckstack.get_correspondingposition(unstuckpos);
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del unstuckstack
        // y un nuevo numero aleatorio entre 0 y 1(que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/unstuckstackpropensity*log(1.0/uprima);

        //reseteamos la alert
        control->restart_unstuckalert();
    }
    else
    {
        double mainstackpropensity= mainstack->get_stackpropensity();
        propensity=mainstackpropensity;
        //cout<<"--propensity   "<<propensity<<endl;

        uniform_real_distribution<double> dist(0.0, mainstackpropensity);

        //miramos el correspondiente posición en el stack principal

        if (control->get_linealsearch()) {mainpos=mainstack->event_lineal_search(dist(mt));}
        else {mainpos=mainstack->event_binary_search(dist(mt));}

        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del mainstack
        // y un nuevo numero aleatorio entre 0 y 1 (que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/mainstackpropensity*log(1.0/uprima);
        //cout<<"--kmctime   "<<kmctime<<endl;

        //aumentamos la alert si el tiempo es muy pqueno

        if (kmctime < control->get_unstuck_compare_time()) control->increaseby1_unstuckalert();
    }
    Event * chosen_event= mainstack->get_event_from_stack(mainpos);

    //sin incremento de tiempo,
    if (chosen_event->get_type()==WITH_OUT_NEIGH_TYPE)
    {
            control->increaseby1_eventsaccepted();
            dissolve_atom_by_position_event_complex(mainpos, box, mainstack, control,tracker);
            control->increaseby1_step();
            return 0;
    }

    //Tenemos que aumentar el tiempo y  eventosaceptados
    control->increment_time_by(kmctime);
    control->increaseby1_eventsaccepted();

    //Tenemos el tiempo a aumentar.. y el evento elegido..
    //Miramos si el evento se lleva realmente a cabo, va a depender de deltaG

    //double ffp=chosen_event->get_ffp();
    //Nuevo
    //double deltaG=chosen_event->get_deltaG();

   // vector <double> Eps=chosen_event->get_Ep_neighbours();

    //Nuevo
    //vector<vector<double>> Ep_direct=chosen_event->get_Ep_neighbours_direct();

    //vector <long long int> Nn=chosen_event->get_count_involved();

    //double ratioevento= mainstack->get_eventrate(mainpos);

    //Nuevo
    //double ratiocaractreforma=Models::model42_reforming(ffp,Nn,Eps,Ep_direct,deltaG);
    //double ratiocaractreforma=Models::model42_reforming(ffp,Nn,Eps,control->get_deltag());

    //double tiempocaractdiss=-log(0.001)/(ratioevento);

    //double probabilidad=1.0-exp(-ratiocaractreforma*tiempocaractdiss);


    //sabemos la probabilidad.. generamos un numero aleatorio entre 0 y 1 y vemos si realmente se va o no
    //double uprimaprima=dist01(mt);


    //if (uprimaprima<probabilidad) //el atomo no se va
    //{
        //double uprima2=dist01(mt); //esto vuelve a  ser para el tiempo pero de la reformacion
        //No puede ser 0
        //if (uprima2<EPSILON*EPSILON*EPSILON*EPSILON) uprima2=EPSILON*EPSILON*EPSILON*EPSILON;

        //el tiempo aumentado va a ser el del resto del sistema total (propsensity) - el del evento de disolución
        //+ el de tiempo de reformación
       // double tiemporeformacion=1.0/(propensity-ratioevento+ratiocaractreforma)*log(1.0/uprima2);
      //  control->increment_time_by(tiemporeformacion);
    //}
    //else //elatomo se va
    //{
        //el tiempo no se aumenta (ya se ha aumentado antes)
        //actualizamos los eventos y la box y surface
        dissolve_atom_by_position_event_complex(mainpos, box, mainstack, control,tracker);
        //debemos actualizar tmb el rate, pq los vecinos han cambiado
    //}
    control->increaseby1_step();

    return 0;
}


long long int Simulation::KMC_step_complex(Event_stack * mainstack, Control * control, Box * box, Tracker * tracker,long long int seed)   //ESTE Y EL ANTERIOR SON LOS QUE SE USAN
{
    //Si el tamaño de los eventos del stack es 0.. ya no tenemos nada que hacer..
    if (mainstack->get_stack_size()==0) return 1; //fail

    double propensity;
    //importamos el inicializador del numero aleatorio (q capullos)
    mt19937 &mt = RandomGenerator::Instance().get();
    mt.seed(seed);
    uniform_real_distribution<double> dist01(0.0, 1.0);
    long long int mainpos;
    double kmctime;

    //Set initial events tmb!! (tendra q ser antes, al comienzo del todo)
    //Seteamos el rate y la propensity en en main stack

    //si no se disuelve en la anterior el rate y la propensity es la misma..

    //mainstack->set_stackpropensity(); //no haria falta pq ya lo hemos definido
    //cout<< " -- mainstack->get_stackpropensity();  " << mainstack->get_stackpropensity() << endl;
    //Si llegamos a que alert es tan alto como unstuckeverystep, creamos el stack unstuck para agilizar
    if (control->get_unstuckalert()==control->get_unstuckeverystep())
    {
        Event_stack_unstuck unstuckstack(IDUNSTUCKSTACK);
        //Un vez creado, creamos los eventos oportunos
        unstuckstack.set_unstuck_events_complex(mainstack, control);
        //si el tamanno del unstuckstack es 0.. la liamos
        if (unstuckstack.get_stack_size()==0) {
            cout<< "unstuckstack without elements!, try to reduce TARGET_TIME or/and increase ACCURACY"<<endl;
            control->restart_unstuckalert();
            //reentramos en el metodo
            KMC_step_complex( mainstack, control, box, tracker,seed);
        }

        //seteamos la propensity de cada evento y la total
        unstuckstack.set_stackpropensity();
        //generamos el numero aleatorio segun la propensity del stack...
        double unstuckstackpropensity= unstuckstack.get_stackpropensity();
        propensity=unstuckstackpropensity;

        //.....
        uniform_real_distribution<double> dist(0.0, unstuckstackpropensity);
        //buscamos segun el valor en el stack unstuck
        long long int unstuckpos;
        if (control->get_linealsearch()) {unstuckpos=unstuckstack.event_lineal_search(dist(mt));}
        else {unstuckpos=unstuckstack.event_binary_search(dist(mt));}
        //miramos el correspondiente posición en el stack principal
        mainpos=unstuckstack.get_correspondingposition(unstuckpos);
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del unstuckstack
        // y un nuevo numero aleatorio entre 0 y 1(que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/unstuckstackpropensity*log(1.0/uprima);

        //reseteamos la alert
        control->restart_unstuckalert();
    }
    else
    {
        double mainstackpropensity= mainstack->get_stackpropensity();
        propensity=mainstackpropensity;
        //cout<<"--propensity   "<<propensity<<endl;

        uniform_real_distribution<double> dist(0.0, mainstackpropensity);
        //miramos el correspondiente posición en el stack principal
        if (control->get_linealsearch()) { mainpos=mainstack->event_lineal_search(dist(mt));}
        else { mainpos=mainstack->event_binary_search(dist(mt));}

        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del mainstack
        // y un nuevo numero aleatorio entre 0 y 1 (que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/mainstackpropensity*log(1.0/uprima);

        //aumentamos la alert si el tiempo es muy pqueno
        if (kmctime < control->get_unstuck_compare_time()) control->increaseby1_unstuckalert();
    }
    Event * chosen_event= mainstack->get_event_from_stack(mainpos);

    //sin incremento de tiempo,
    if (chosen_event->get_type()==WITH_OUT_NEIGH_TYPE)
    {
            control->increaseby1_eventsaccepted();
            dissolve_atom_by_position_event_complex(mainpos, box, mainstack, control,tracker);
            control->increaseby1_step();
            return 0;
    }
    //Tenemos que aumentar el tiempo y  eventosaceptados
    control->increment_time_by(kmctime);
    control->increaseby1_eventsaccepted();

    //Tenemos el tiempo a aumentar.. y el evento elegido..
    //Miramos si el evento se lleva realmente a cabo, va a depender de deltaG

    //double ffp=chosen_event->get_ffp();
    //Nuevo
    //double deltaG=chosen_event->get_deltaG();
    //vector<vector<double>> Ep_direct= chosen_event->get_Ep_neighbours_direct();

    //vector <double> Eps=chosen_event->get_Ep_neighbours();
    //vector <long long int> Nn=chosen_event->get_count_involved();

    //double ratioevento= mainstack->get_eventrate(mainpos);

    //Nuevo

    //double ratiocaractreforma=Models::model42_reforming(ffp,Nn,Eps,Ep_direct,deltaG);
    //double ratiocaractreforma=Models::model42_reforming(ffp,Nn,Eps,control->get_deltag());

    //double tiempocaractdiss=-log(0.001)/(ratioevento);

    //double probabilidad=1.0-exp(-ratiocaractreforma*tiempocaractdiss);

    //sabemos la probabilidad.. generamos un numero aleatorio entre 0 y 1 y vemos si realmente se va o no
    //double uprimaprima=dist01(mt);

    //if (uprimaprima<probabilidad) //el atomo no se va
    //{
        //cout<<"--aqui no deberia llegar "<<endl;
        //double uprima2=dist01(mt); //esto vuelve a  ser para el tiempo pero de la reformacion
        //No puede ser 0
        //if (uprima2<EPSILON*EPSILON*EPSILON*EPSILON) uprima2=EPSILON*EPSILON*EPSILON*EPSILON;

        //el tiempo aumentado va a ser el del resto del sistema total (propsensity) - el del evento de disolución
        //+ el de tiempo de reformación
        //double tiemporeformacion=1.0/(propensity-ratioevento+ratiocaractreforma)*log(1.0/uprima2);
        //control->increment_time_by(tiemporeformacion);
    //}
    //else //elatomo se va
    //{
        //el tiempo no se aumenta (ya se ha aumentado antes)
        //actualizamos los eventos y la box y surface
        dissolve_atom_by_position_event_complex(mainpos, box, mainstack, control,tracker);
        //debemos actualizar tmb el rate, pq los vecinos han cambiado
    //}

    control->increaseby1_step();

    return 0;
}




long long int Simulation::KMC_step(Event_stack * mainstack, Control * control, Box * box)
{
    //Si el tamaño de los eventos del stack es 0.. ya no tenemos nada que hacer..
    if (mainstack->get_stack_size()==0) return 1; //fail

    double propensity;
    //importamos el inicializador del numero aleatorio (q capullos)
    mt19937 &mt = RandomGenerator::Instance().get();
    uniform_real_distribution<double> dist01(0.0, 1.0);
    long long int mainpos;
    double kmctime;
    //Set initial events tmb!! (tendra q ser antes, al comienzo del todo)
    //Seteamos el rate y la propensity en en main stack

    //si no se disuelve en la anterior el rate y la propensity es la misma..

    //mainstack->set_stackpropensity(); //no haria falta pq ya lo hemos definido
    //cout<< " -- mainstack->get_stackpropensity();  " << mainstack->get_stackpropensity() << endl;
    //Si llegamos a que alert es tan alto como unstuckeverystep, creamos el stack unstuck para agilizar
    if (control->get_unstuckalert()==control->get_unstuckeverystep())
    {
        Event_stack_unstuck unstuckstack(IDUNSTUCKSTACK);
        //Un vez creado, creamos los eventos oportunos
        unstuckstack.set_unstuck_events(mainstack, control);
        //seteamos la propensity de cada evento y la total
        unstuckstack.set_stackpropensity();
        //generamos el numero aleatorio segun la propensity del stack...
        double unstuckstackpropensity= unstuckstack.get_stackpropensity();
        propensity=unstuckstackpropensity;

        //.....
        uniform_real_distribution<double> dist(0.0, unstuckstackpropensity);
        //buscamos segun el valor en el stack unstuck
        long long int unstuckpos;
        if (control->get_linealsearch()) {unstuckpos=unstuckstack.event_lineal_search(dist(mt));}
        else {unstuckpos=unstuckstack.event_binary_search(dist(mt));}
        //miramos el correspondiente posición en el stack principal
        mainpos=unstuckstack.get_correspondingposition(unstuckpos);
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del unstuckstack
        // y un nuevo numero aleatorio entre 0 y 1(que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/unstuckstackpropensity*log(1.0/uprima);

        //reseteamos la alert
        control->restart_unstuckalert();
    }
    else
    {
        double mainstackpropensity= mainstack->get_stackpropensity();
        propensity=mainstackpropensity;
        //cout<<"--propensity   "<<propensity<<endl;

        uniform_real_distribution<double> dist(0.0, mainstackpropensity);
        //miramos el correspondiente posición en el stack principal
        if (control->get_linealsearch()) {mainpos=mainstack->event_lineal_search(dist(mt));}
        else {mainpos=mainstack->event_binary_search(dist(mt));}
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del mainstack
        // y un nuevo numero aleatorio entre 0 y 1 (que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/mainstackpropensity*log(1.0/uprima);
        //cout<<"--kmctime   "<<kmctime<<endl;

        //aumentamos la alert si el tiempo es muy pqueno
        if (kmctime < control->get_unstuck_compare_time()) control->increaseby1_unstuckalert();
    }
    //Tenemos que aumentar el tiempo y  eventosaceptados
    control->increment_time_by(kmctime);
    control->increaseby1_eventsaccepted();

    //Tenemos el tiempo a aumentar.. y el evento elegido..
    //Miramos si el evento se lleva realmente a cabo, va a depender de deltaG
    long long int vecinos = mainstack->get_involved_first_atom(mainpos)->get_size_neighbour();
    double ratioevento= mainstack->get_eventrate(mainpos);

    double ratiocaractreforma=Models::model42_reforming(control->get_fr(),vecinos,control->get_Er(),control->get_deltag());

    double tiempocaractdiss=-log(0.001)/(ratioevento);

    double probabilidad=1.0-exp(-ratiocaractreforma*tiempocaractdiss);

    //sabemos la probabilidad.. generamos un numero aleatorio entre 0 y 1 y vemos si realmente se va o no
    double uprimaprima=dist01(mt);


    if (uprimaprima<probabilidad) //el atomo no se va
    {
        cout<<"--aqui no deberia llegar "<<endl;
        double uprima2=dist01(mt); //esto vuelve a  ser para el tiempo pero de la reformacion
        //No puede ser 0
        if (uprima2<EPSILON*EPSILON*EPSILON*EPSILON) uprima2=EPSILON*EPSILON*EPSILON*EPSILON;

        //el tiempo aumentado va a ser el del resto del sistema total (propsensity) - el del evento de disolución
        //+ el de tiempo de reformación
        double tiemporeformacion=1.0/(propensity-ratioevento+ratiocaractreforma)*log(1.0/uprima2);
        control->increment_time_by(tiemporeformacion);
    }
    else //elatomo se va
    {
        //el tiempo no se aumenta (ya se ha aumentado antes)
        //actualizamos los eventos y la box y surface
        dissolve_atom_by_position_event(mainpos, box, mainstack, control);
        //debemos actualizar tmb el rate, pq los vecinos han cambiado
    }

    control->increaseby1_step();

    return 0;
}



long long int Simulation::KMC_step_resolution(Event_stack * mainstack, Control * control, Box * box)
{
    //Si el tamaño de los eventos del stack es 0.. ya no tenemos nada que hacer..
    if (mainstack->get_stack_size()==0) return 1; //fail

    double propensity;
    //importamos el inicializador del numero aleatorio (q capullos)
    mt19937 &mt = RandomGenerator::Instance().get();
    uniform_real_distribution<double> dist01(0.0, 1.0);
    long long int mainpos;
    double kmctime;
    //Set initial events tmb!! (tendra q ser antes, al comienzo del todo)
    //Seteamos la propensity en en main stack
    //mainstack->set_stackpropensity();

    //Si llegamos a que alert es tan alto como unstuckeverystep, creamos el stack unstuck para agilizar
    if (control->get_unstuckalert()==control->get_unstuckeverystep())
    {
        Event_stack_unstuck unstuckstack(IDUNSTUCKSTACK);
        //Un vez creado, creamos los eventos oportunos
        unstuckstack.set_unstuck_events(mainstack, control);
        //seteamos la propensity de cada evento y la total
        unstuckstack.set_stackpropensity();
        //generamos el numero aleatorio segun la propensity del stack...
        double unstuckstackpropensity= unstuckstack.get_stackpropensity();
        propensity=unstuckstackpropensity;

        //.....
        uniform_real_distribution<double> dist(0.0, unstuckstackpropensity);
        //buscamos segun el valor en el stack unstuck
        long long int unstuckpos;
        if (control->get_linealsearch()) {unstuckpos=unstuckstack.event_lineal_search(dist(mt));}
        else {unstuckpos=unstuckstack.event_binary_search(dist(mt));}
        //miramos el correspondiente posición en el stack principal
        mainpos=unstuckstack.get_correspondingposition(unstuckpos);
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del unstuckstack
        // y un nuevo numero aleatorio entre 0 y 1(que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/unstuckstackpropensity*log(1.0/uprima);

        //reseteamos la alert
        control->restart_unstuckalert();
    }
    else
    {
        double mainstackpropensity= mainstack->get_stackpropensity();
        propensity=mainstackpropensity;

        uniform_real_distribution<double> dist(0.0, mainstackpropensity);
        //miramos el correspondiente posición en el stack principal
        if (control->get_linealsearch()) {mainpos=mainstack->event_lineal_search(dist(mt));}
        else {mainpos=mainstack->event_binary_search(dist(mt));}
        //miramos el tiempo que se debería de aumentar el sistema de acuerdo con la propensity del mainstack
        // y un nuevo numero aleatorio entre 0 y 1 (que no puede ser 0.0)
        double uprima=dist01(mt);
        if (uprima<EPSILON*EPSILON*EPSILON*EPSILON) uprima=EPSILON*EPSILON*EPSILON*EPSILON;
        kmctime=1.0/mainstackpropensity*log(1.0/uprima);

        //aumentamos la alert si el tiempo es muy pqueno
        if (kmctime < control->get_unstuck_compare_time()) control->increaseby1_unstuckalert();
    }
    //Tenemos que aumentar el tiempo y  eventosaceptados
    control->increment_time_by(kmctime*(control->get_time_modifier()));
    control->increaseby1_eventsaccepted();

    //Tenemos el tiempo a aumentar.. y el evento elegido..
    //Miramos si el evento se lleva realmente a cabo, va a depender de deltaG
    long long int vecinos = mainstack->get_involved_first_atom(mainpos)->get_size_neighbour();
    double ratioevento= mainstack->get_eventrate(mainpos);

    double ratiocaractreforma=Models::model42_reforming(control->get_fr(),vecinos,control->get_Er(),control->get_deltag());

    double tiempocaractdiss=-log(0.001)/(ratioevento);

    double probabilidad=1.0-exp(-ratiocaractreforma*tiempocaractdiss);

    //sabemos la probabilidad.. generamos un numero aleatorio entre 0 y 1 y vemos si realmente se va o no
    double uprimaprima=dist01(mt);


    if (uprimaprima<probabilidad) //el atomo no se va
    {
        double uprima2=dist01(mt); //esto vuelve a  ser para el tiempo pero de la reformacion
        //No puede ser 0
        if (uprima2<EPSILON*EPSILON*EPSILON*EPSILON) uprima2=EPSILON*EPSILON*EPSILON*EPSILON;

        //el tiempo aumentado va a ser el del resto del sistema total (propsensity) - el del evento de disolución
        //+ el de tiempo de reformación
        double tiemporeformacion=1.0/(propensity-ratioevento+ratiocaractreforma)*log(1.0/uprima2);
        control->increment_time_by(tiemporeformacion*(control->get_time_modifier()));
    }
    else //elatomo se va
    {
        //el tiempo no se aumenta (ya se ha aumentado antes)
        //actualizamos los eventos y la box y surface
        dissolve_atom_by_position_event(mainpos, box, mainstack, control);
    }

    control->increaseby1_step();

    return 0;
}

long long int Simulation::complete_simulation(string name_input_file)
{
        cout << "Welcome to Kimera Program" << endl;

        string * pointer_name_input_file=& name_input_file;

        //Neceistamos: Celda inicial, data_printer, el data_reader, el control, la box, el stack de eventos?, y el tracker?
        Event_stack event_stack1=Event_stack(1);
        Event_stack * pointer_event_stack1= &event_stack1;

        Tracker tracker1= Tracker();
        Tracker * pointer_tracker1=& tracker1;

        Cell firstcell=Cell(-1);
        Cell * pointer_firstcell=& firstcell;

        Control control1=Control();
        Control * pointer_control1= & control1;

        Data_reader reader1=Data_reader();

        reader1.read_input_file(pointer_name_input_file,pointer_control1);

        //No paralelizamos la simulacion si es

        if(!control1.get_parallelize())    {omp_set_num_threads(1);}//  PARA QUITAR LA PARALELIZACION
        else
        {
            long long int concurentThreadsSupported = thread::hardware_concurrency();
            long long int cores=control1.get_parallelizenumber();

            if (cores>concurentThreadsSupported)
            {
                cout<< endl <<"WARNING: You have asked for  " <<cores<< " cores, nevertheless the maximum number of cores available in your system are  "<< concurentThreadsSupported <<endl;
                cores=concurentThreadsSupported;
            }

            cout<< endl <<"The simulation is allocated in "<< cores <<" cores"<<endl<< endl;
            omp_set_num_threads(cores);
        }

        //necesitabamos el nombre del work
        Data_printer printer1(control1.get_workname());

        //Ya puedo crear la Box

        Box box1(control1.get_dimension_x(), control1.get_dimension_y(), control1.get_dimension_z(),
        control1.get_cell_a(), control1.get_cell_b(), control1.get_cell_c(),
        control1.get_angle_alpha(), control1.get_angle_beta(), control1.get_angle_gamma());

        Box * pointer_box1 = &box1;

        //Si hemos dicho que lo leemos de un archivo kimera..

        if (control1.get_kimerafile())
        {
            string file_kimera_box_name=control1.get_pathtokimerafile();
            string * pointer_file_kimera_box_name=& file_kimera_box_name;

            reader1.read_kimera_box(pointer_file_kimera_box_name,pointer_box1);
            box1.clean_gaps_before_start_complex_from_kimerafile();

            //HAY QUE HACER MAS COSAS
            //Hay que saber los eventos

            reader1.read_input_file_dis_events(pointer_name_input_file,pointer_control1);

            //surface inicial / no hace falta pq ya debería estar definida
            //box1.set_initial_surface_atoms(control1.get_list_event_definition());

            //se modifica la surface initial
            //reader1.read_remove_add_from_surface(pointer_name_input_file,pointer_box1);

            //eventos iniciales

            event_stack1.set_initial_events_complex(pointer_box1,pointer_control1);
        }
        else
        {
            //puede pasar que haya leido un archivo xyz, y le quiera annadir unas posiciones
            if (control1.get_xyzfile())
            {
                string file_xyz_name=control1.get_pathtoxyzfile();
                string * pointer_file_xyz_name=& file_xyz_name;

                reader1.read_create_xyz_file_without_corners(pointer_file_xyz_name,pointer_firstcell,pointer_control1);
            }

            //leo la celda para definir la box
            reader1.read_input_file_cell(pointer_name_input_file,pointer_firstcell);
            //Si el tamano de las posiciones de la celda es 0 despues de leer el fichero.. se lia
            if (firstcell.get_cell_positions().size()==0){cout<<endl<< "Initial cell empty! " <<endl; exit(1);}


            //creo la box a partir de la celda, mirar el seedboxbool
            if(control1.get_seedboxbool()) {box1.create_box_by_cell(pointer_firstcell,control1.get_seedbox());}
            else {box1.create_box_by_cell(pointer_firstcell);}

            //AQUI VENDRIA MODIFICAR ELEMENTOS
            reader1.read_modify_element(pointer_name_input_file,pointer_box1);


            //Depende de las condicioner periodicas, los vecinos van a ser unos u otros
            box1.set_cell_neighbours_26(control1.get_plane_x(),control1.get_plane_y(),control1.get_plane_z());


            //leemos los eventos..son necesarios para definir los vecinos
            reader1.read_input_file_dis_events(pointer_name_input_file,pointer_control1);


            //definimos los vecinos segun la cantidad de event Definition
            vector <Event_definition*> events= control1.get_list_event_definition();
            long long int sizeevents=events.size();

            for (long long int i=0; i<sizeevents;i++)
            {

                if (control1.get_linked())
                {

                    box1.set_atom_neighbourhood_linked(events.at(i)->get_involved_atom_type(), events.at(i)->get_type_neighbours(), events.at(i)->get_distance_neighbours(), events.at(i)->get_list_linked_neighbour()) ;

                }
                else
                {
                    box1.set_atom_neighbourhood(events.at(i)->get_involved_atom_type(), events.at(i)->get_type_neighbours(), events.at(i)->get_distance_neighbours());
                }

                if (control1.get_affected())
                {
                    box1.set_atom_affected(events.at(i)->get_involved_atom_type(),
                       events.at(i)->get_type_affected(),
                       events.at(i)->get_distance_affected());
                }

            }

            //una vez tenemos todos los vecinos, definimos esu record
            long long int sizebox=box1.get_box_atoms().size();
            #pragma omp parallel for
            for (long long int  i=0;i<sizebox;i++)
            {
                box1.get_box_atoms().at(i)->set_neigh_record();
            }

            //ahora que tenemos la box y sus vecinos, vamos a modificarla
            reader1.read_input_file_topography(pointer_name_input_file,pointer_box1);

/**
            cout<<"VAMOS A VER AQUI QUE PASA CON LOS VECINOS Y LOS LINKED"<<endl;

            for (long long int  i=0;i<sizebox;i++)
            {
                cout<<"atomo numero: "<< i <<endl;
                cout<<"cantidad de vecinos: "<<box1.get_box_atoms().at(i)->get_size_neighbour()<<endl;
                cout<<"cantidad de linked: "<<box1.get_box_atoms().at(i)->get_linked().size()<<endl;
                for (long long int j=0;j<box1.get_box_atoms().at(i)->get_linked().size();j++)
                {
                    cout<< box1.get_box_atoms().at(i)->get_linked().at(j).size() <<endl;
                }
            }
*/
                            //printer1.print_all_data(pointer_box1,pointer_control1,pointer_tracker1); //QUITAR
                            //printer1.print_initial_box(pointer_box1,"initial");//QUITAR



            //liberamos los espacios de la box  que quedan vacios. Si es distinto de 0, en alguno ha habido error
            box1.clean_gaps_before_start_complex();

            //printer1.print_complete_system(pointer_box1,pointer_control1);



            //Ahora definimos la surface inicial

            box1.set_initial_surface_atoms(control1.get_list_event_definition());

            reader1.read_remove_add_from_surface(pointer_name_input_file,pointer_box1);


            //AQUI UN BUEN REMOVE FROM SURFACE!!!!

            //printer1.print_all_data(pointer_box1,pointer_control1,pointer_tracker1); //BORRAR

            //ahora los eventos iniciales

            event_stack1.set_initial_events_complex(pointer_box1,pointer_control1);
        }


        //Se define la masa de los atomos

        reader1.read_input_file_mass(pointer_name_input_file,pointer_box1);


        //si hemos dicho que se estime el tiempo..

        if (control1.get_estimatetime())
        {
            //¿porque se queda corto el tiempo? pq estamos calculando al final el tiempo de un atomo para coord 1? seguro?
            //podemos jugar con el rate del stack?
            //deberia poner una funcion en control que sea llamada cada vez que defino un evento

            //vamos a intentar lo siguiente: buscamos el segundo evento que menos rate tenga, y ese le usamos multiplicado por 5, por ejemplo


            vector<double> list_rates;

            bool invector=false;


            for (long long int h=0; h<(long long int)event_stack1.get_stack_size();h++)
            {
                for (long long int f=0; f<(long long int)list_rates.size(); f++)
                {
                    if (Util::isEqualhigh(event_stack1.get_eventrate(h),list_rates.at(f)))
                    {
                        invector=true;
                        break;
                    }
                }

                if (!invector)
                {
                    list_rates.push_back(event_stack1.get_eventrate(h));
                }
                invector=false;
            }

            sort(list_rates.begin(), list_rates.end());

            int pos=(int)(list_rates.size()/2.0);  //si coge el mas pequenno mejor, es el que mas tarda
            control1.set_targettime((double)box1.get_box_atoms().size()/(list_rates.at(pos)*(double)event_stack1.get_stack_size()));
            cout << endl<<"Estimated time for dissolution: "<<(double)box1.get_box_atoms().size()/(list_rates.at(pos)*(double)event_stack1.get_stack_size())<<" seconds"<<endl;


            /**
            control1.set_targettime(box1.get_box_atoms().size()*event_stack1.get_stack_size()/event_stack1.get_stackpropensity());
            cout << endl<<"Estimated time for dissolution: "<<box1.get_box_atoms().size()*event_stack1.get_stack_size()/event_stack1.get_stackpropensity()<<" seconds"<<endl;
            */
        }

        cout<<endl<<endl<<"You managed to reach this point!!!, Starting the simulation..."<<endl<<endl;

        //Creo que tenemos todo, imprimimos las cosas:



        if (control1.get_printinitialkimerafile()) printer1.print_initial_box(pointer_box1,"initial",control1.get_affected(), control1.get_linked());

        if (!control1.get_print_by_steps())
        {
            printer1.print_all_data(pointer_box1,pointer_control1,pointer_tracker1);

            while(control1.get_time()<control1.get_targettime())
            {
                //try{
                    //mirar el seedsimbool
                    if  (control1.get_seedsimbool())
                    {
                         if(KMC_step_complex(pointer_event_stack1, pointer_control1,pointer_box1,pointer_tracker1,control1.get_seedsim())  !=0) break;
                    }
                    else
                    {
                        if(KMC_step_complex(pointer_event_stack1, pointer_control1,pointer_box1,pointer_tracker1)  !=0) break;
                    }

                    printer1.print_all_data(pointer_box1,pointer_control1,pointer_tracker1);
                //}

                //catch (exception& e)
                //{
                //   cout<<endl<<"There has been an error during the simulation :'( "<<endl;
                //   cout<<endl<<"You can send a copy of your input file to the author of this \"not as wonderful as expected\" code"<<endl;
                //   if (control1.get_printfinalkimerafile()) printer1.print_initial_box(pointer_box1,"final",control1.get_affected(), control1.get_linked());
                //   exit(1);
                //}
            }
        }

        else
        {
            printer1.print_all_data_step(pointer_box1,pointer_control1,pointer_tracker1);

            while(control1.get_step()<control1.get_targetsteps())
            {
                //try{
                    //mirar el seedsimbool
                    if  (control1.get_seedsimbool())
                    {
                         if(KMC_step_complex(pointer_event_stack1, pointer_control1,pointer_box1,pointer_tracker1,control1.get_seedsim())  !=0) break;
                    }
                    else
                    {
                        if(KMC_step_complex(pointer_event_stack1, pointer_control1,pointer_box1,pointer_tracker1)  !=0) break;
                    }

                    printer1.print_all_data_step(pointer_box1,pointer_control1,pointer_tracker1);
                //}

                //catch (exception& e)
                //{
                //   cout<<endl<<"There has been an error during the simulation :'( "<<endl;
                //   cout<<endl<<"You can send a copy of your input file to the author of this \"not as wonderful as expected\" code"<<endl;
                //   if (control1.get_printfinalkimerafile()) printer1.print_initial_box(pointer_box1,"final",control1.get_affected(), control1.get_linked());
                //   exit(1);
                //}
            }
        }

        //imprimir kimerabox si el sistema tiene algo y lo hemos pedido

        cout<<endl<<"Simulation finished!!"<<endl<<endl;

       if (control1.get_printfinalkimerafile()) printer1.print_initial_box(pointer_box1,"final",control1.get_affected(), control1.get_linked());

       cout<<endl<<"Compressing output files..."<<endl<<endl;
       //zip the output files

       string phrase="zip kimeraoutput.zip "+control1.get_workname()+".*";

        system(phrase.c_str());

        string phrase2=REMOVE+control1.get_workname()+".*";

        //si exite el fichero kimeraoutput.zip, borramos el resto de archivos

        ifstream filezip("kimeraoutput.zip");

        if (filezip.is_open())
        {
            system(phrase2.c_str());
        }
        else
        {
            cout<<endl<<"Unable to zip the output files. You need zip command in your system!"<<endl<<endl;
        }


        return 0;
}
