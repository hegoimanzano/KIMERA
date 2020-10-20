#ifndef ATOM_H
#define ATOM_H


#include "Util.h"
#include "Record.h"
#include "Linked_neighbour.h"
#include <vector>
#include <string>

#define INDISLOCATION  2
#define NORMAL 1
#define FREEPOSITION -1  //-1 type of atom is reserved for atoms that are going to be removed before start the simulation
#define INSOLUBLE -2  //-2 type reserved for insoluble atoms
#define DISSOLVED -3  //-3 type reserved for dissolved atoms
#define REMOVED -4  //-4 type reserved for removed before starting the simulation





using namespace std;

class Atom
{
    public:
        Atom( long long int id_,double x_,double y_,double z_,long long int type_,bool insurface_);
        Atom( long long int id_,double x_,double y_,double z_,string atom_type_,bool insurface_);
        virtual ~Atom();

        void set_id(long long int id_);
        long long int get_id();

        void set_x(double x_);
        double get_x();

        void set_y(double y_);
        double get_y();

        void set_z(double z_);
        double get_z();

        void set_type(long long int type_);
        long long int get_type();

        void set_atom_type(string atom_type_);
        string get_atom_type();

        void set_mass(double mass_);
        double get_mass();

        void set_insurface(bool insurface_);
        bool get_insurface();

        double get_distance(Atom *otheratom);

        long long int add_neighbour(Atom *neighbour);
        long long int get_size_neighbour();
        Atom* get_neighbour(long long int pos_);
        long long int rm_neighbour(long long int id);
        long long int rm_allneighbour();

        long long int get_id_neighbour(long long int pos_);
        bool get_insurface_neighbour(long long int pos_);
        double get_x_neighbour(long long int pos_);
        double get_y_neighbour(long long int pos_);
        double get_z_neighbour(long long int pos_);
        long long int get_type_neighbour(long long int pos_);

        vector<Atom*> get_neighbours();

        vector<double> get_distances_to_neighbours();
        long long int add_distance_to_neighbours(double distance_);
        double get_distance_to_neighbours_by_pos(long long int pos_);

        long long int set_neigh_record();          //retorna size atom_record
        vector<Record> get_neigh_record();

        long long int rm_neighbour_and_record(long long int id_neigh);

        long long int rm_neighbour_and_record_both_directions(long long int id_neigh);

        long long int reduce_record_count_by_1(string type_, double distance_);

        long long int get_number(string type_, double distance_);


        vector<Atom*> get_affected();
        long long int add_affected(Atom * afected_);
        long long int rm_affected(long long int id);

        vector<vector<Atom*>> get_linked();
        long long int add_linked_to(Atom * linked_,long long int pos_);
        long long int rm_linked(long long int pos_);

        vector<long long int> get_linked_type(){return linked_type;}
        long long int add_linked_type(long long int linked_type_){linked_type.push_back(linked_type_);return linked_type.size();}

        void add_void_vector_in_linked();
        void remove_neighbours_with_out_linkeds();

        bool is_bulk(vector<long long int> max_to_bulk, vector<string>types_, vector<double>distances_);
        bool is_bulk_with_links(vector<long long int> max_to_bulk, vector<string>types_, vector<double>distances_, vector<Linked_neighbour> link_vector);

        long long int get_number_complex(  string type_event_definition, double distance_event_definition, long long int typelink_from_event_def);


    protected:

    private:
        long long int * id;
        double * x;
        double * y;
        double * z;
        long long int * type;
        string * atom_type;
        bool * insurface;
        vector<Atom*> neighbours;
        //en pareja sincronizado con los neinghbours
        vector<vector<Atom*>> linked;
        vector<long long int> linked_type;

        vector<Atom*> affected;
        vector<double> distances_to_neighbours;

        vector<Record> neigh_record;

        double * mass;
};

#endif // ATOM_H
