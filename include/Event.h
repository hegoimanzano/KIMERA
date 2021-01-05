#ifndef EVENT_H
#define EVENT_H

#include "Linked_neighbour.h"
#include "Util.h"
#include <vector>
#include <string>


#define DISSOLUTION_TYPE 1
#define WITH_OUT_NEIGH_TYPE 2

using namespace std;

class Atom;

class Event
{
    public:
        Event(long long int id_, long long int type_, double rate, Atom* involvedatom);
        virtual ~Event();

        long long int get_id();
        long long int get_type();
        void set_type(long long int type_);

        void set_propensity(double propensity_);
        double get_propensity();

        double get_rate();
        void set_rate(double rate_);

        vector<Atom*> get_involved_atoms();
        void add_involved_atom(Atom *atom);
        long long int get_size_involved_atom();


        vector <double> get_distance_neighbours();
        long long int set_distance_neighbours(vector <double> distance_neighbours_);

        vector <string> get_type_neighbours();
        long long int set_type_neighbours(vector <string> type_neighbours_);

        vector <Linked_neighbour> get_list_linked_neighbour();
        long long int set_list_linked_neighbour(vector<Linked_neighbour> list_linked_neighbour_);

        vector <double> get_Ed_neighbours();
        long long int set_Ed_neighbours(vector <double> Ed_neighbours_);

        vector <double> get_Ep_neighbours();
        long long int set_Ep_neighbours(vector <double> Ep_neighbours_);


        long long int set_Ed_neighbours_direct(vector<vector<double>> Ed_neighbours_direct_);
        vector<vector<double>> get_Ed_neighbours_direct();
        vector<double> get_Ed_neighbours_direct_by_pos(long long int pos_);
        long long int get_size_Ed_neighbours_direct_by_pos(long long int pos_);     //vamos a usar este como criterio, si el tamaño de este es 0, quiere decir que es Lineal en vez de Directo
        double get_Ed_neighbours_direct_by_posij(long long int posi_,long long int posj_);


        long long int set_Ep_neighbours_direct(vector<vector<double>> Ep_neighbours_direct_);
        vector<vector<double>> get_Ep_neighbours_direct();
        vector<double> get_Ep_neighbours_direct_by_pos(long long int pos_);
        long long int get_size_Ep_neighbours_direct_by_pos(long long int pos_);
        double get_Ep_neighbours_direct_by_posij(long long int posi_,long long int posj_);


        string get_involved_atom_type();
        void set_involved_atom_type(string involved_atom_type_);

        vector <long long int> get_count_involved();
        long long int set_count_involved(vector <long long int> count_involved_);


        void set_ffd(double ffd_);
        double get_ffd();

        void set_ffp(double ffp_);
        double get_ffp();

        void set_deltaG(double deltaG_);
        double get_deltaG();

    protected:

    private:

        long long int id;
        long long int type;
        vector<Atom*> involved_atoms;
        double propensity;
        double rate;

        string involved_atom_type;              //Tipo de átomo que se ve involucrado en el evento
        vector <double> distance_neighbours;    //Distancia a Vecinos  que condicionan en evento
        vector <string> type_neighbours;        //tipo de vecinos que condicionan el evento
        vector <Linked_neighbour> list_linked_neighbour;

        vector <double> Ed_neighbours;           //Energia de disolucion respecto a los vecinos
        vector <double> Ep_neighbours;           //Energia de precipitacion respecto a los vecinos

        vector <vector<double>> Ed_neighbours_direct;                 //en caso de que sea modelo directo. cada vecino tiene que ir aumentando conforme lo hace Ed_neighbours
        vector <vector<double>> Ep_neighbours_direct;                 //en caso de que sea modelo directo. cada vecino tiene que ir aumentando conforme lo hace Ep_neighbours


        vector<long long int> count_involved;

        double ffd;                              //frecuencia fundamental de disolcion del evento
        double ffp;
        double deltaG;

};

#endif // EVENT_H
