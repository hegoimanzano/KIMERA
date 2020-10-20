#include "Linked_neighbour.h"


Linked_neighbour::Linked_neighbour()
{
    //ctor
}

Linked_neighbour::~Linked_neighbour()
{
    //dtor
}

long long int Linked_neighbour::add_distance_origin_linked(double distance_origin_linked_)
{
    distance_origin_linked.push_back(distance_origin_linked_);
    return distance_origin_linked.size();
}


long long int Linked_neighbour::add_type_linked(string type_linked_)
{
    type_linked.push_back(type_linked_);
    return type_linked.size();
}

long long int Linked_neighbour::add_distance_target_linked(double distance_target_linked_)
{
    distance_target_linked.push_back(distance_target_linked_);
    return distance_target_linked.size();
}

long long int Linked_neighbour::add_linked_neighbour(string type_linked_,double distance_origin_linked_, double distance_target_linked_)
{
    type_linked.push_back(type_linked_);
    distance_target_linked.push_back(distance_target_linked_);
    distance_origin_linked.push_back(distance_origin_linked_);
    return distance_target_linked.size();
}


