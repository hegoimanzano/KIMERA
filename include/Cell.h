#ifndef CELL_H
#define CELL_H

#include "Atom.h"
#include "Position.h"
#include "Util.h"
#include <vector>
#include <string>
#include <math.h>
#include <iostream>

using namespace std;

class Cell
{
    public:
        Cell(long long int cell_id_);
        virtual ~Cell();

        long long int get_cell_id();

        long long int add_atom(Atom * incellatom);
        long long int rm_atom(long long int id);


        long long int add_position(double x_, double y_, double z_,string type_, double prob_);
        long long int rm_position(long long int pos_id_);

        //metodo que compara posiciones con las ya existentes y si ya existe, le añade probabilidad y tipo, y si no, la crea directamente
        long long int compare_and_add(double x_, double y_, double z_,string type_, double prob_);        //retorna 0 si no existia
                                                                                                                    //retorna 1 si ya existia
        Position * get_pos_in(double x_,double y_, double z_);

        vector<Atom*> get_cell_atoms();

        vector<Position*> get_cell_positions();

        vector<Position*> get_cell_removed_positions();

        void set_cornerx(long long int cornerx);
        long long int get_cornerx(); //1, 0 o -1

        void set_cornery(long long int cornery);
        long long int get_cornery();

        void set_cornerz(long long int cornerz);
        long long int get_cornerz();


        long long int add_cell_neighbour(Cell* cell);

        vector<Cell*> get_cell_neighbour();


    protected:

    private:

        long long int cell_id;

        long long int cornerx=0;
        long long int cornery=0;
        long long int cornerz=0;

        //String spatial_group;
        vector<Atom*> cell_atoms;
        vector<Position*> cell_positions;
        vector<Position*> cell_removed_positions;

        vector<Cell*> cell_neighbours;

};

#endif // CELL_H
