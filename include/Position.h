#ifndef POSITION_H
#define POSITION_H

#include <algorithm>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include "Util.h"

using namespace std;

class Position
{
    public:
        Position(double x_,double y_,double z_, long long int position_id_, string type_, double prob_);
        virtual ~Position();

        double get_x();
        double get_y();
        double get_z();
        long long int get_position_id();
        string get_type(long long int pos);
        double get_prob(long long int pos);

        long long int get_number_possible_types();

        vector<string> get_types();
        vector<double> get_probs();

        long long int add_type_prob(string type, double prob);
        long long int rm_type_prob(string type);

        bool compare_position(Position * position_);

        double get_sum_probs();

    protected:

    private:
        long long int position_id;
        double x;
        double y;
        double z;

        vector<string> types;
        vector<double> probs;

};

#endif // POSITION_H
