#include "Tracker.h"

Tracker::Tracker()
{
    //ctor
}

Tracker::~Tracker()
{
    //dtor
}

long long int Tracker::add_atom_dissolved(Atom *atom)
{
    records_dissolved.push_back(atom->get_neigh_record());
    types_of_record.push_back(atom->get_atom_type());

    string targettype=atom->get_atom_type();

    long long int sizetypesdissolved=types_dissolved_atoms.size();
    bool invector=false;
    long long int posinvector=0;

    for (long long int i=0; i<sizetypesdissolved; i++)
    {
        if (targettype.compare(types_dissolved_atoms.at(i))==0) {invector=true; posinvector=i; break;}
    }

    if (invector)
    {
        number_dissolved_atoms.at(posinvector)++;
    }
    else
    {
        types_dissolved_atoms.push_back(targettype);
        number_dissolved_atoms.push_back(1);
    }


    return types_dissolved_atoms.size();
}


double Tracker::get_mean(string targettype_, string neightype_, double neightdistance_)
{
    long long int sizerecordsdis=records_dissolved.size();
    long long int counter=0;
    long long int amount=0;


    #pragma omp parallel for shared (counter, amount)
    for (long long int i=0;i<sizerecordsdis; i++)
    {

        if (types_of_record.at(i).compare(targettype_)==0)
        {
            long long int targetrecordsize=records_dissolved.at(i).size();
            for (long long int j=0;j<targetrecordsize; j++)
            {
               Record targetrecord=records_dissolved.at(i).at(j);
               if(targetrecord.get_type().compare(neightype_)==0 && Util::isEqualDistances(targetrecord.get_distance(),neightdistance_))
               {
                   #pragma omp critical
                   {
                        counter++;
                        amount+=targetrecord.get_number() ;
                   }
               }
            }
        }
    }

    if (counter==0) return 0;
    return (double)amount/(double)counter;
}



vector<string> Tracker::get_types_dissolved_atoms()
{
    return types_dissolved_atoms;
}
vector<long long int> Tracker::get_number_dissolved_atoms()
{
    return number_dissolved_atoms;
}
vector<Record> Tracker::get_records_dissolved(long long int position)
{
    return records_dissolved.at(position);
}

long long int Tracker::get_number(string targettype_)
{
     long long int value;
    volatile bool flag=false;

     long long int sizetypesdis=types_dissolved_atoms.size();

    //#pragma omp parallel for shared(flag)
    for ( long long int i=0; i<sizetypesdis;i++)
    {
        if (flag) continue;
        if (targettype_.compare(types_dissolved_atoms.at(i))==0)
        {
            value=i;
            flag=true;
        }
    }
    if (flag) return number_dissolved_atoms.at(value);
    return 0;
}
