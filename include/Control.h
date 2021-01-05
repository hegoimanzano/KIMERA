#ifndef CONTROL_H
#define CONTROL_H


#define BOLTZCONSTANT 1.38064852e-23
#define PLANCK 6.62607004e-34
#define COORDINATIONDISLOCA 4.0

#include <math.h>
#include <string>
#include "Util.h"
#include "Event_definition.h"


using namespace std;

class Control
{
    public:
        Control();
        virtual ~Control();


     long long int get_dimension_x();
    void set_dimension_x( long long int dimension_x_);

     long long int get_dimension_y();
    void set_dimension_y( long long int dimension_y_);

     long long int get_dimension_z();
    void set_dimension_z( long long int dimension_z_);

    double get_cell_a();
    void set_cell_a(double cell_a_);

    double get_cell_b();
    void set_cell_b(double cell_b_);

    double get_cell_c();
    void set_cell_c(double cell_c_);

    double get_angle_gamma();
    void set_angle_gamma(double angle_gamma_);

    double get_angle_alpha();
    void set_angle_alpha(double angle_alpha_);

    double get_angle_beta();
    void set_angle_beta(double angle_beta_);

    //XYZ UVW OPERATIONS

    double from_uvw_to_x(double u, double v, double w);
    double from_uvw_to_y(double v, double w);
    double from_uvw_to_z(double w);

    double from_xyz_to_u(double x, double y, double z);
    double from_xyz_to_v(double y, double z);
    double from_xyz_to_w(double z);


    double get_box_boundx();
    void set_box_boundx();

    double get_box_boundx0();
    void set_box_boundx0();

    double get_box_boundy();
    void set_box_boundy();

    double get_box_boundy0();
    void set_box_boundy0();

    double get_box_boundz();
    void set_box_boundz();

    double get_box_boundz0();
    void set_box_boundz0();

    double get_box_boundxl();
    void set_box_boundxl();

    double get_box_boundxl0();
    void set_box_boundxl0();

    double get_box_boundyl();
    void set_box_boundyl();

    double get_box_boundyl0();
    void set_box_boundyl0();

    double get_box_boundzl();
    void set_box_boundzl();

    double get_box_boundzl0();
    void set_box_boundzl0();

    void set_tilt_factors();
    double get_tilt_factors_xy();
    double get_tilt_factors_xz();
    double get_tilt_factors_yz();

    void set_tilt_bounds_factors();

    double get_temperature();
    void set_temperature(double temperature_);

    double get_deltag();
    void set_deltag(double deltag_);

    double get_fd();
    void set_fd(double fd_);
    double set_fd_from_T();

    double get_fr();
    void set_fr(double fr_);

    double get_Ed();
    void set_Ed(double Ed_);

    double get_Er();
    void set_Er(double Er_);


    double get_time_modifier();
    void set_time_modifier(double time_modifier_);

    double get_scale_parameter_m1();
    void set_scale_parameter_m1(double scale_parameter_m1_);

    double get_scale_parameter_m2();
    void set_scale_parameter_m2(double scale_parameter_m2_);


    double get_time();
    void set_time(double time_);
    void increment_time_by(double increment);

    double get_targettime();
    void set_targettime(double targettime_);

    double get_resolution();
    void   set_resolution(double resolution_);

    long long int add_list_event_definition(Event_definition * Event_definition_); //return size
    Event_definition* get_event_definition_by_pos(long long int pos_);
    vector<Event_definition*> get_list_event_definition();


    long long int get_eventsaccepted();
    long long int increaseby1_eventsaccepted();

    double change_distance_accuracy(double distance_accuracy_);
    double get_distance_accuracy();

    double get_unstuckaccuracy();
    void set_unstuckaccuracy(double unstuckaccuracy_);

    long long int get_unstuckeverystep();
    void set_unstuckeverystep(long long int unstuckeverystep_);


    long long int get_unstuckalert();
    void set_unstuckalert(long long int unstuckalert_);
    long long int increaseby1_unstuckalert();
    long long int restart_unstuckalert();


    long long int get_filesdatanumber();
    void set_filesdatanumber(long long int filesdatanumber_);
    //TRUE FALSE
    //---------------

    bool get_estimatetime();
    void set_estimatetime(bool estimatetime_);

    bool get_randomize();//esto demomento siempre es true
    void set_randomize(bool randomize_); //esto demomento siempre es true

    bool get_dataanalisys();
    void set_dataanalisys(bool dataanalisys_);

    bool get_layerxanalysis();
    void set_layerxanalysis(bool layerxanalysis_);

    bool get_layeryanalysis();
    void set_layeryanalysis(bool layeryanalysis_);

    bool get_layerzanalysis();
    void set_layerzanalysis(bool layerzanalysis_);

    bool get_coordinationanalysis();
    void set_coordinationanalysis(bool coordinationanalysis_);

    bool get_boxframe();
    void set_boxframe(bool boxframe_);

    bool get_surfaceframe();
    void set_surfaceframe(bool surfaceframe_);

    bool get_initialdata();
    void set_initialdata(bool initialdata_);

    bool get_xyzfile();
    void set_xyzfile(bool xyzfile_);

    string get_pathtoxyzfile();
    void set_pathtoxyzfile(string pathtoxyzfile_);

    bool get_kimerafile();
    void set_kimerafile(bool kimerafile_);

    string get_pathtokimerafile();
    void set_pathtokimerafile(string pathtokimerafile_);

    bool get_printinitialkimerafile();
    void set_printinitialkimerafile(bool printinitialkimerafile_);

    bool get_printfinalkimerafile();
    void set_printfinalkimerafile(bool printfinalkimerafile_);

    bool get_print_by_steps();
    void set_print_by_steps(bool print_by_steps_);

    long long int get_targetsteps();
    void set_targetsteps(long long int targetsteps_);

    bool get_parallelize();
    void set_parallelize(bool parallelize_);

    long long int get_parallelizenumber();
    void set_parallelizenumber(long long int  parallelizenumber_);

    bool get_linealsearch();
    void set_linealsearch(bool linealsearch_);

    bool get_meandiscoord();
    void set_meandiscoord(bool meandiscoord_);

    bool get_affected();
    void set_affected(bool affected_);

    bool get_linked();
    void set_linked(bool linked_);

    string get_workname();
    void set_workname(string workname_);

    //---------------
    //---------------

    double get_unstuck_compare_time();

    double estimate_time();

    //----Auxiliar printing parameters----//

     long long int get_step();
     long long int increaseby1_step();


    bool get_grain();
    void set_grain(bool grain_);

    bool get_plane_x();
    void set_plane_x(bool plane_x);
    bool get_plane_y();
    void set_plane_y(bool plane_y);
    bool get_plane_z();
    void set_plane_z(bool plane_z);



    double get_gyradiousini();
    void set_gyradiousini(double gyradiousini_);

     long long int get_maxcoordini();
    void set_maxcoordini( long long int maxcoordini_);

     long long int get_totalnumberpositions();
    long long int increaseby1_totalnumberpositions();
    long long int decreaseby1_totalnumberpositions();

    long long int get_datasteps();
    void set_datastepts(long long int datasteps_);

    long long int get_layersteps();
    void set_layersteps(long long int layersteps_);

    long long int get_boxsteps();
    void set_boxsteps(long long int boxsteps_);

    long long int get_surfacesteps();
    void set_surfacesteps(long long int surfacesteps_);

    long long int get_meansteps();
    void set_meansteps(long long int meansteps_);

    bool get_seedboxbool();
    void set_seedboxbool(bool seedboxbool_);

    long long int get_seedbox();
    void set_seedbox(long long int seedbox_);

    bool get_seedsimbool();
    void set_seedsimbool(bool seedsimbool_);

    long long int get_seedsim();
    void set_seedsim(long long int seedsim_);

     long long int get_id_linked();
     void set_id_linked(long long int id_);


    protected:

    private:

    //----Box parameters----//
     long long int dimension_x;
     long long int dimension_y;
     long long int dimension_z;
    double cell_a;
    double cell_b;
    double cell_c;
    double angle_gamma;   //angle xy
    double angle_alpha;      //angle yz
    double angle_beta;      //angle xz


    double box_boundx0;
    double box_boundy0;
    double box_boundz0;
    double box_boundx;
    double box_boundy;
    double box_boundz;

    double box_boundxl0;
    double box_boundyl0;
    double box_boundzl0;
    double box_boundxl;
    double box_boundyl;
    double box_boundzl;


    double tilt_factor_xy;
    double tilt_factor_xz;
    double tilt_factor_yz;

    //----Physics parameters----//
    double temperature;
    double deltag;
    double fd;
    double fr;
    double Ed;
    double Er;

    double time_modifier;
    double scale_parameter_m1;
    double scale_parameter_m2;

    double time;

    double targettime;

    double resolution;

    vector<Event_definition*> list_event_definition;

    //----Auxiliar parameters----//
    string workname;

    long long int eventsaccepted;

    long long int unstuckeverystep;
    double unstuckaccuracy;
    long long int unstuckalert;

    bool estimatetime;
    bool randomize; //de momento simpre true

    long long int filesdatanumber;

    bool dataanalisys;
    long long int datasteps;

    bool layerxanalysis;
    bool layeryanalysis;
    bool layerzanalysis;
    long long int layersteps;

    bool coordinationanalysis; //creo que quitar

    bool boxframe;
    long long int boxsteps;

    bool surfaceframe;
    long long int surfacesteps;

    bool meandiscoord;
    long long int meansteps;

    bool initialdata;           //tenemos esto implementado?

    bool parallelize;
    long long int parallelizenumber;

    bool xyzfile;
    string pathtoxyzfile;

    bool kimerafile;
    string pathtokimerafile;

    bool linealsearch;

    bool printinitialkimerafile;
    bool printfinalkimerafile;

    bool print_by_steps;
    long long int targetsteps;

    //----Auxiliar printing parameters----//

    long long int step;
    bool grain;
    bool planex;
    bool planey;
    bool planez;

    double gyradiousini; //creo que quitar

     long long int maxcoordini; //creo que quitar

     long long int initialnumberatoms; //creo que quitar

     long long int totalnumberpositions; //Para el id de las posiciones creo que quitar

    vector <string> types_initial;
    vector <long long int> number_initial;

    bool seedboxbool;
    long long int seedbox;

    bool seedsimbool;
    long long int seedsim;

    bool linked;
    bool affected;

    long long int id_linked;
};

#endif // CONTROL_H
