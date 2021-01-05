#include "Position.h"

Position::Position(double x_,double y_,double z_, long long int position_id_,string type_, double prob_)
{
    x=x_;
    y=y_;
    z=z_;
    position_id=position_id_;
    types.push_back(type_);
    probs.push_back(prob_);
}

Position::~Position()
{
    //dtor
}

double Position::get_x()
{
    return x;
}

double Position::get_y()
{
    return y;
}

double Position::get_z()
{
    return z;
}

long long int Position::get_position_id()
{
    return position_id;
}

long long int Position::get_number_possible_types()
{
    return types.size();
}

vector<string> Position::get_types()
{
    return types;
}

vector<double> Position::get_probs()
{
    return probs;
}


string Position::get_type(long long int pos)
{
    return types.at(pos);
}

double Position::get_prob(long long int pos)
{
    return probs.at(pos);
}

long long int Position::add_type_prob(string type, double prob)
{
    types.push_back(type);
    probs.push_back(prob);
    return 0;
}

long long int Position::rm_type_prob(string type)
{
    vector<string>::iterator it = find(types.begin(), types.end(), type);


    if (it == types.end())
    {
        return 1; //not found
    }

    long long int i = distance(types.begin(), it);

    types.erase(types.begin()+i);
    probs.erase(probs.begin()+i);

    return 0;
}


bool Position::compare_position(Position * position_)
{
    double x_=position_->get_x();
    double y_=position_->get_y();
    double z_=position_->get_z();

    if (Util::isEqual(x_,x) && Util::isEqual(y_,y) && Util::isEqual(z_,z)) return true;

    return false;
}

double Position::get_sum_probs()
{
    double sum=0.0;
    long long int probsize=probs.size();
    for (long long int i=0;i<probsize;i++)
    {
        sum+=probs.at(i);
    }
    return sum;
}
