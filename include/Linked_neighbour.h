#ifndef LINKED_NEIGHBOUR_H
#define LINKED_NEIGHBOUR_H

#include <string>
#include <vector>

#define NO_LINKED_TYPE 0
#define LINKED_TYPE_NORMAL 1
#define LINKED_TYPE_DISSOLVED 2


using namespace std;

class Linked_neighbour
{
    public:
        Linked_neighbour();
        virtual ~Linked_neighbour();

        long long int get_id(){return id;}
        void set_id(long long int id_){id=id_;}

        long long int get_type() { return type; }
        void set_type(long long int val) { type = val; }
        vector <string> get_type_linked() { return type_linked; }
        void set_type_linked(vector <string> val) { type_linked = val; }
        vector <double> get_distance_origin_linked() { return distance_origin_linked; }
        void set_distance_origin_linked(vector <double> val) { distance_origin_linked = val; }
        vector <double> get_distance_target_linked() { return distance_target_linked; }
        void set_distance_target_linked(vector <double> val) { distance_target_linked = val; }

        long long int add_type_linked(string type_linked_);
        long long int add_distance_origin_linked(double distance_origin_linked_);
        long long int add_distance_target_linked(double distance_target_linked_);

        //retorna la cantidad actual de los atomos linkeados
        long long int add_linked_neighbour(string type_linked_, double distance_origin_linked_, double distance_target_linked_); //retorna la cantidad actual de los atomos linkeados

    protected:

    private:
        long long int id;
        long long int type;
        vector <string> type_linked;
        vector <double> distance_origin_linked;
        vector <double> distance_target_linked;
};

#endif // LINKED_NEIGHBOUR_H
