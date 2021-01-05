#ifndef DATA_PRINTER_H
#define DATA_PRINTER_H

#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include "Control.h"
#include "Box.h"
#include "Cell.h"
#include "Tracker.h"

using namespace std;

class Data_printer
{
    public:
        Data_printer(string filename_);
        virtual ~Data_printer();

        string get_filename();
        //set_filename(char * filename_);

        long long int print_all_data(Box * box, Control *control, Tracker * tracker);
        long long int print_all_data_step(Box * box, Control *control, Tracker * tracker);

        long long int print_state(Box * box, Control *control); //Para guardar el estado... TODO
        long long int print_initial_box(Box * box, string name, bool affected_, bool links_);  //kimera format
        long long int print_complete_system(Box * box, Control *control);
        long long int print_data(Box * box, Control *control,Tracker * tracker);
        long long int print_surface(Box * box, Control *control);
        long long int print_coord(Box * box, Control *control); //de momento nos olvidamos de este
        long long int print_mean_dissolved_coordination( Control *control, Tracker * tracker);

        long long int print_xlayer(Box * box, Control *control);
        long long int print_ylayer(Box * box, Control *control);
        long long int print_zlayer(Box * box, Control *control);

        double get_surface_dispersion_x(Box * box);
        double get_surface_dispersion_y(Box * box);
        double get_surface_dispersion_z(Box * box);

        double get_grain_gyradious(Box * box);
        double get_grain_alpha(Box * box, Control *control);
        double get_grain_alpha(Box * box, Control *control, double gyradious);

        long long int set_initial_type_number(Box * box);      //return 0 si el tamanno total concuerda con el tamanno de la box

        long long int get_number_inbox_initial(string targettype_);

    private:
        string filename;

        long long int previousremoved;
        long long int previouspreviousremoved;

        //double previousfraction; creo que no hace falta
        //double previouspreviousfraction;

        double previousdispersionx;
        double previouspreviousdispersionx;

        double previousdispersiony;
        double previouspreviousdispersiony;

        double previousdispersionz;
        double previouspreviousdispersionz;

        double previousalpha;
        double previouspreviousalpha;

        double previousgyradious;
        double previouspreviousgyradious;

        double previoustime;
        double previousprevioustime;

        //Para imprimir el tiempo cada un cierto numero de pasos
        double controldata;
        double previouscontroldata;

        double controlbox;
        double previouscontrolbox;

        double controllayer;
        double previouscontrollayer;

        double controlsurface;
        double previouscontrolsurface;

        double controlmean;
        double previouscontrolmean;

        //Para imprimir con cada un cierto numero de pasos
        long long int controldatastep;
        long long int previouscontroldatastep;

        long long int controlboxstep;
        long long int previouscontrolboxstep;

        long long int controllayerstep;
        long long int previouscontrollayerstep;

        long long int controlsurfacestep;
        long long int previouscontrolsurfacestep;

        long long int controlmeanstep;
        long long int previouscontrolmeanstep;


        vector<string> typesinbox_initial;
        vector<long long int> numberinbox_initial;

};

#endif // DATA_PRINTER_H
