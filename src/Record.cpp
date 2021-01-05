#include "Record.h"

Record::Record(string type_, double distance_)
{
    number=1;
    type=type_;
    distance=distance_;
}

Record::~Record()
{
    //dtor
}


long long int Record::get_number()
{
    return number;
}
void Record::set_number(long long int number_)
{
    number=number_;
}
long long int Record::increment_number_by_1()
{
    number++;
    return number;
}

long long int Record::reduce_number_by_1()
{
    number--;
    return number;
}

string Record::get_type()
{
    return type;
}

void Record::set_type(string type_)
{
    type=type_;
}

double Record::get_distance()
{
    return distance;
}
void Record::set_distance(double distance_)
{
    distance=distance_;
}
