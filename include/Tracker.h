#ifndef TRACKER_H
#define TRACKER_H

#include <string>
#include <vector>
#include "Atom.h"
#include "Record.h"
#include "Util.h"

using namespace std;

class Tracker
{
    public:
        Tracker();
        virtual ~Tracker();

        long long int add_atom_dissolved(Atom *atom); //annade a type si no lo hay,si lo hay, sube la cuenta del type, y annade el record al vector de records

        vector<string> get_types_dissolved_atoms();
        vector<long long int> get_number_dissolved_atoms();
        vector<Record> get_records_dissolved(long long int position);


        double get_mean(string targettype_, string neightype_, double neightdistance_);

        long long int get_number(string targettype_);



    protected:

    private:

        vector <string> types_dissolved_atoms;
        vector <long long int> number_dissolved_atoms;

        vector <vector<Record>> records_dissolved;
        vector <string> types_of_record;
};

#endif // TRACKER_H
