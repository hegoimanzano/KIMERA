#ifndef EVENT_DEFINITION_H
#define EVENT_DEFINITION_H

#include <vector>
#include <string>
#include "Linked_neighbour.h"


using namespace std;

class Event_definition
{
    public:
        Event_definition();
        virtual ~Event_definition();

        string get_involved_atom_type();
        void set_involved_atom_type(string involved_atom_type_);

        long long int add_distance_neighbours(double distance_);           //retorna tamanno del vector
        double get_distance_neighbours_by_pos(long long int pos_);         //retorna distancia a vecino en posicion n

        long long int add_type_neighbours(string type_neighbours_);
        string get_type_neighbours_by_pos(long long int pos_);

        long long int add_Ed_neighbours(double Ed_neighbours_);
        double get_Ed_neighbours_by_pos(long long int pos);

        long long int add_Ep_neighbours(double Ep_neighbours_);
        double get_Ep_neighbours_by_pos(long long int pos);

        vector<vector<double>> get_Ed_neighbours_direct();
        long long int add_Ed_neighbours_direct(vector<double> Ed_neighbour_direct_);
        vector<double> get_Ed_neighbours_direct_by_pos(long long int pos_);
        long long int get_size_Ed_neighbours_direct_by_pos(long long int pos_);     //vamos a usar este como criterio, si el tamaño de este es 0, quiere decir que es Lineal en vez de Directo
        double get_Ed_neighbours_direct_by_posij(long long int posi_,long long int posj_);

        vector<vector<double>> get_Ep_neighbours_direct();
        long long int add_Ep_neighbours_direct(vector<double> Ep_neighbour_direct_);
        vector<double> get_Ep_neighbours_direct_by_pos(long long int pos_);
        long long int get_size_Ep_neighbours_direct_by_pos(long long int pos_);
        double get_Ep_neighbours_direct_by_posij(long long int posi_,long long int posj_);



        long long int add_distance_affected(double distance_affected_);   //retorna tamanno del vector
        double get_distance_affected_by_pos(long long int pos_);         //retorna distancia a afectado en posicion n
        vector<double> get_distance_affected();         //retorna distancia a afectado en posicion n

        long long int add_type_affected(string type_affected_);   //retorna tamanno del vector
        string get_type_affected_by_pos(long long int pos_);
        vector<string> get_type_affected();


        ////////////// LINK PART

        vector <Linked_neighbour> get_list_linked_neighbour();
        long long int add_to_list_linked_neighbour(Linked_neighbour linked_neighbour_);

        /////////////

        void set_deltaG(double deltaG_);
        double get_deltaG();

        void set_ffd(double ffd_);
        double get_ffd();

        void set_ffp(double ffp_);
        double get_ffp();

        long long int add_max_to_bulk(long long int max_to_bulk_);            //retorna tamanno, si esta en 0, no influye
        vector <long long int> get_max_to_bulk();

        vector <double> get_distance_neighbours();
        vector <string> get_type_neighbours();

        vector <double> get_Ed_neighbours();
        vector <double> get_Ep_neighbours();






    protected:

    private:
        //¿La definicion de un evento puede estar compuesta por combinacion de vecinos de forma directa + vecinos de forma lineal?

        //long long int event_definition_type;    // Tipo del evento del que se trata. Puede ser Lineal (Con el numero de vecinos) O modelo directo | esto no vale, vamos a seguir el criterio de tamaño de Ed_neighbours_direct = 0

        string involved_atom_type;              //Tipo de átomo que se ve involucrado en el evento


        vector <double> distance_neighbours;    //Distancia a Vecinos  que condicionan en evento

        vector <string> type_neighbours;        //tipo de vecinos que condicionan el evento

        vector <double> Ed_neighbours;           //Energia de disolucion respecto a los vecinos

        vector <double> Ep_neighbours;           //Energia de precipitacion respecto a los vecinos


        vector <vector<double>> Ed_neighbours_direct;                 //en caso de que sea modelo directo. cada vecino tiene que ir aumentando conforme lo hace Ed_neighbours

        vector <vector<double>> Ep_neighbours_direct;                 //en caso de que sea modelo directo. cada vecino tiene que ir aumentando conforme lo hace Ep_neighbours


        vector <double> distance_affected;       //Distancia a cosas que se ven afectados por el evento

        vector <string> type_affected;           //tipo de cosas que se ven afectado por el evento

        vector <Linked_neighbour> list_linked_neighbour;




        double ffd;                              //frecuencia fundamental de disolcion del evento
        double ffp;                              //frecuencia fundamental de reformacion del evento
        double deltaG;                           //deltaG para el evento

        vector <long long int> max_to_bulk;                //maximos vecinos con los que se considera que esta en bulk,

};

#endif // EVENT_DEFINITION_H
