#ifndef DATA_READER_H
#define DATA_READER_H


#include <stdexcept> //Errors
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <vector>
#include "Control.h"
#include "Sym_equation.h"
#include "Cell.h"
#include "Position.h"
#include "Event_definition.h"
#include "Box.h"
#include "Atom.h"


using namespace std;

class Data_reader
{
    public:
        Data_reader();
        virtual ~Data_reader();

        long long int read_create_xyz_file(string * filename_, Cell * cell);
        long long int read_create_xyz_file_without_corners(string * filename_, Cell * cell, Control * control);
        long long int read_create_cif_file(string * filename_,Control * control, Cell * cell);

        long long int read_input_file(string * filename_, Control * control); //modifica control (lo primero), al final de leer llamariamos un metodo que lo setee

        long long int read_input_file_cell(string * filename_, Cell * cell); //define la primera celda con la que luego haremos la box


        long long int read_input_file_topography(string * filename_, Box * box); // una vez definida la box, la cambiamos

                                                                        // luego llamariamos al metodo que define la superficie

                                                                        // y luego se definirian los eventos inciales

        long long int read_input_file_mass(string * filename_, Box * box); // una vez definida la box, se otorga masa a los atomos


        long long int read_input_file_dis_events(string * filename_, Control * control);

        long long int read_kimera_box(string * filename_, Box * box);

        long long int read_remove_add_from_surface(string * filename_, Box * box);  // una vez definido la surface, removemos distinas cosas de la superficie

        long long int read_modify_element(string * filename_, Box * box);  //una vez definida la box, se pueden cambiar las particulas de la box definida


    protected:

    private:


};

#endif // DATA_READER_H
