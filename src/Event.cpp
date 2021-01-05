#include "Event.h"

Event::Event(long long int id_, long long int type_,double rate_, Atom* involvedatom)
{
    id=id_;
    type=type_;
    rate=rate_;
    add_involved_atom(involvedatom);
}

Event::~Event()
{

}

long long int Event::get_id()
{
    return id;
}

long long int Event::get_type()
{
    return type;
}

void Event::set_type(long long int type_)
{
    type=type_;
}

double Event::get_propensity()
{
    return propensity;
}

void Event::set_propensity(double propensity_)
{
    propensity=propensity_;
}

double Event::get_rate()
{
    return rate;
}

void Event::set_rate(double rate_)
{
    rate=rate_;
}

void Event::add_involved_atom(Atom *atom)
{
    involved_atoms.push_back(atom);
}

vector<Atom*> Event::get_involved_atoms()
{
    return involved_atoms;
}

long long int Event::get_size_involved_atom()
{
    return involved_atoms.size();
}



vector <double> Event::get_distance_neighbours()
{
    return distance_neighbours;
}

long long int Event::set_distance_neighbours(vector <double> distance_neighbours_)
{
    distance_neighbours=distance_neighbours_;
    return distance_neighbours.size();
}

vector <string> Event::get_type_neighbours()
{
    return type_neighbours;
}

long long int Event::set_type_neighbours(vector <string> type_neighbours_)
{
    type_neighbours=type_neighbours_;
    return type_neighbours.size();

}

vector <Linked_neighbour>  Event::get_list_linked_neighbour()
{
    return list_linked_neighbour;
}

long long int Event::set_list_linked_neighbour(vector <Linked_neighbour> list_linked_neighbour_)
{
    list_linked_neighbour=list_linked_neighbour_;
    return list_linked_neighbour.size();
}

vector <double> Event::get_Ed_neighbours()
{
    return Ed_neighbours;
}

long long int Event::set_Ed_neighbours(vector <double> Ed_neighbours_)
{
    Ed_neighbours=Ed_neighbours_;
    return Ed_neighbours.size();
}

vector <double> Event::get_Ep_neighbours()
{
    return Ep_neighbours;
}

long long int Event::set_Ep_neighbours(vector <double> Ep_neighbours_)
{
    Ep_neighbours=Ep_neighbours_;
    return Ep_neighbours.size();
}


long long int Event::set_Ed_neighbours_direct(vector<vector<double>> Ed_neighbours_direct_)
{
    Ed_neighbours_direct=Ed_neighbours_direct_;
    return Ed_neighbours_direct.size();
}

vector<vector<double>> Event::get_Ed_neighbours_direct()
{
    return Ed_neighbours_direct;
}

vector<double> Event::get_Ed_neighbours_direct_by_pos(long long int pos_)
{
    return Ed_neighbours_direct.at(pos_);
}

long long int Event::get_size_Ed_neighbours_direct_by_pos(long long int pos_)     //vamos a usar este como criterio, si el tamaño de este es 0, quiere decir que es Lineal en vez de Directo
{
    return Ed_neighbours_direct.at(pos_).size();
}

double Event::get_Ed_neighbours_direct_by_posij(long long int posi_,long long int posj_)
{
    return Ed_neighbours_direct.at(posi_).at(posj_);
}



long long int Event::set_Ep_neighbours_direct(vector<vector<double>> Ep_neighbours_direct_)
{
    Ep_neighbours_direct=Ep_neighbours_direct_;
    return Ep_neighbours_direct.size();
}

vector<vector<double>> Event::get_Ep_neighbours_direct()
{
    return Ep_neighbours_direct;
}

vector<double> Event::get_Ep_neighbours_direct_by_pos(long long int pos_)
{
    return Ep_neighbours_direct.at(pos_);
}

long long int Event::get_size_Ep_neighbours_direct_by_pos(long long int pos_)    //vamos a usar este como criterio, si el tamaño de este es 0, quiere decir que es Lineal en vez de Directo
{
    return Ep_neighbours_direct.at(pos_).size();
}

double Event::get_Ep_neighbours_direct_by_posij(long long int posi_,long long int posj_)
{
    return Ep_neighbours_direct.at(posi_).at(posj_);
}



string Event::get_involved_atom_type()
{
    return involved_atom_type;
}

void Event::set_involved_atom_type(string involved_atom_type_)
{
    involved_atom_type=involved_atom_type_;
}

vector <long long int> Event::get_count_involved()
{
    return count_involved;
}

long long int Event::set_count_involved(vector <long long int> count_involved_)
{
    count_involved=count_involved_;
    return count_involved.size();
}

void Event::set_ffd(double ffd_)
{
    ffd=ffd_;
}

double Event::get_ffd()
{
    return ffd;
}

void Event::set_ffp(double ffp_)
{
    ffp=ffp_;
}

double Event::get_ffp()
{
    return ffp;
}

void Event::set_deltaG(double deltaG_)
{
    deltaG=deltaG_;
}

double Event::get_deltaG()
{
    return deltaG;
}
