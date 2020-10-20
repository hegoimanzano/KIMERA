#include "Event_definition.h"

Event_definition::Event_definition()
{

}

Event_definition::~Event_definition()
{

}

string Event_definition::get_involved_atom_type()
{
    return involved_atom_type;
}

void Event_definition::set_involved_atom_type(string involved_atom_type_)
{
    involved_atom_type=involved_atom_type_;
}

long long int Event_definition::add_distance_neighbours(double distance_)           //return vector size
{
    distance_neighbours.push_back(distance_);
    return distance_neighbours.size();
}

double Event_definition::get_distance_neighbours_by_pos(long long int pos_)         //return distance to neighbor in pos
{
    return distance_neighbours.at(pos_);
}

long long int Event_definition::add_type_neighbours(string type_neighbours_)
{
    type_neighbours.push_back(type_neighbours_);
    return type_neighbours.size();
}
string Event_definition::get_type_neighbours_by_pos(long long int pos_)
{
    return type_neighbours.at(pos_);
}

long long int Event_definition::add_Ed_neighbours(double Ed_neighbours_)
{
    Ed_neighbours.push_back(Ed_neighbours_);
    return Ed_neighbours.size();
}

double Event_definition::get_Ed_neighbours_by_pos(long long int pos)
{
    return Ed_neighbours.at(pos);
}

long long int Event_definition::add_Ep_neighbours(double Ep_neighbours_)
{
    Ep_neighbours.push_back(Ep_neighbours_);
    return Ep_neighbours.size();
}

double Event_definition::get_Ep_neighbours_by_pos(long long int pos)
{
    return Ep_neighbours.at(pos);
}

vector<vector<double>> Event_definition::get_Ed_neighbours_direct()
{
    return Ed_neighbours_direct;
}

long long int Event_definition::add_Ed_neighbours_direct(vector<double> Ed_neighbour_direct_)
{
    Ed_neighbours_direct.push_back(Ed_neighbour_direct_);
    return Ed_neighbours_direct.size();
}

vector<double> Event_definition::get_Ed_neighbours_direct_by_pos(long long int pos_)
{
    return Ed_neighbours_direct.at(pos_);
}

long long int Event_definition::get_size_Ed_neighbours_direct_by_pos(long long int pos_)
{
    return Ed_neighbours_direct.at(pos_).size();
}

double Event_definition::get_Ed_neighbours_direct_by_posij(long long int posi_,long long int posj_)
{
    return Ed_neighbours_direct.at(posi_).at(posj_);
}

vector<vector<double>> Event_definition::get_Ep_neighbours_direct()
{
    return Ep_neighbours_direct;
}

long long int Event_definition::add_Ep_neighbours_direct(vector<double> Ep_neighbour_direct_)
{
    Ep_neighbours_direct.push_back(Ep_neighbour_direct_);
    return Ep_neighbours_direct.size();
}

vector<double> Event_definition::get_Ep_neighbours_direct_by_pos(long long int pos_)
{
    return Ep_neighbours_direct.at(pos_);
}

long long int Event_definition::get_size_Ep_neighbours_direct_by_pos(long long int pos_)
{
    return Ep_neighbours_direct.at(pos_).size();
}

double Event_definition::get_Ep_neighbours_direct_by_posij(long long int posi_,long long int posj_)
{
    return Ep_neighbours_direct.at(posi_).at(posj_);
}


long long int Event_definition::add_distance_affected(double distance_affected_)  //return vector size
{
    distance_affected.push_back(distance_affected_);
    return distance_affected.size();
}

double Event_definition::get_distance_affected_by_pos(long long int pos_)      //returna distance to  affected in n position
{
    return distance_affected.at(pos_);
}
vector<double> Event_definition::get_distance_affected()
{
    return distance_affected;
}

long long int Event_definition::add_type_affected(string type_affected_)   //return vector size
{
    type_affected.push_back(type_affected_);
    return type_affected.size();
}
string Event_definition::get_type_affected_by_pos(long long int pos_)   //return distance to affected
{
    return type_affected.at(pos_);
}
vector<string> Event_definition::get_type_affected()
{
    return type_affected;
}

long long int Event_definition::add_to_list_linked_neighbour(Linked_neighbour Linked_neighbour_)
{
    list_linked_neighbour.push_back(Linked_neighbour_);
    return list_linked_neighbour.size();
}

vector <Linked_neighbour> Event_definition::get_list_linked_neighbour()
{
    return list_linked_neighbour;
}

void Event_definition::set_deltaG(double deltaG_)
{
    deltaG=deltaG_;
}

double Event_definition::get_deltaG()
{
    return deltaG;
}


void Event_definition::set_ffd(double ffd_)
{
    ffd=ffd_;
}
double Event_definition::get_ffd()
{
    return ffd;
}

void Event_definition::set_ffp(double ffp_)
{
    ffp=ffp_;
}
double Event_definition::get_ffp()
{
    return ffp;
}

long long int Event_definition::add_max_to_bulk(long long int max_to_bulk_)      //return vector size
{
    max_to_bulk.push_back(max_to_bulk_);
    return max_to_bulk.size();
}

vector <long long int> Event_definition::get_max_to_bulk()
{
    return max_to_bulk;
}

vector <double> Event_definition::get_distance_neighbours()
{
    return distance_neighbours;
}

vector <string> Event_definition::get_type_neighbours()
{
    return type_neighbours;
}

vector <double> Event_definition::get_Ed_neighbours()
{
    return Ed_neighbours;
}

vector <double> Event_definition::get_Ep_neighbours()
{
    return Ep_neighbours;
}




