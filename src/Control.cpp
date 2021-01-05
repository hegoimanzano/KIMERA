#include "Control.h"

Control::Control()
{


        //----Box parameters----//
    dimension_x=10;
    dimension_y=10;
    dimension_z=10;
    cell_a=2.5;
    cell_b=2.5;
    cell_c=2.5;
    angle_gamma=90.0;   //angle xy
    angle_alpha=90.0;      //angle yz
    angle_beta=90.0;      //angle xz

    box_boundx=dimension_x*cell_a;
    box_boundy=dimension_y*cell_b;
    box_boundz=dimension_z*cell_c;

    box_boundxl=dimension_x*cell_a;
    box_boundyl=dimension_y*cell_b;
    box_boundzl=dimension_z*cell_c;

    box_boundx0=0.0;
    box_boundy0=0.0;
    box_boundz0=0.0;

    box_boundxl0=0.0;
    box_boundyl0=0.0;
    box_boundzl0=0.0;

    tilt_factor_xy=0.0;
    tilt_factor_xz=0.0;
    tilt_factor_yz=0.0;


    //----Physics parameters----//
    temperature=300.0;
    deltag=-200.0;
    fd=temperature*BOLTZCONSTANT/PLANCK;
    fr=fd;
    Ed=7.47;
    Er=4.95;

    scale_parameter_m1=53.0;
    scale_parameter_m2=1.68;

    time_modifier=1.0;

    time=0.0;

    targettime=1.0;
    targetsteps=1000;

    resolution=1.0;

    //----Auxiliar parameters----//
    //uniform_real_distribution<double> dis
    eventsaccepted=0;

    unstuckeverystep=100;
    unstuckaccuracy=1.0e38;
    unstuckalert=0;



    print_by_steps=false;
    estimatetime=false;
    randomize=true;


    filesdatanumber=20;

    dataanalisys=false;
    layerxanalysis=false;
    layeryanalysis=false;
    layerzanalysis=false;
    coordinationanalysis=false;
    boxframe=false;
    surfaceframe=false;
    initialdata=false;
    meandiscoord=false;
    linealsearch=false;

    linked=false;
    affected=false;

    datasteps=20;
    layersteps=20;
    boxsteps=20;
    surfacesteps=20;
    meansteps=20;

    //Auxiliar Printing parameters

    step=0;
    grain=true; //The simulation is a grain case
    planex=false; //The simulation is a infinity plane case
    planey=false;
    planez=false;

    gyradiousini=0.0;

    maxcoordini=0;

    totalnumberpositions=0;

    initialnumberatoms=0;

    seedboxbool=false;
    seedbox=0;

    seedsimbool=false;
    seedsim=0;


    parallelize=false;
    parallelizenumber=0;

    xyzfile=false;
    pathtoxyzfile="";

    kimerafile=false;
    pathtokimerafile="";

    printinitialkimerafile=false;
    printfinalkimerafile=false;
    id_linked=0;

    workname="Kimera_Simulation";


}

Control::~Control()
{
/**
    //destructor de los Event_definition
    long long int limit=list_event_definition.size();

    for (long long int i=0;i<limit;i++)
    {
        Event_definition * todelete = list_event_definition.at(i);
        delete todelete;
    }//dtor
    */
}
//-----------------------------------------
 long long int Control::get_dimension_x()
{
    return dimension_x;
}
void Control::set_dimension_x( long long int dimension_x_)
{
    dimension_x=dimension_x_;
}
//-----------------------------------------
 long long int Control::get_dimension_y()
{
    return dimension_y;
}
void Control::set_dimension_y( long long int dimension_y_)
{
    dimension_y=dimension_y_;
}
//-----------------------------------------
 long long int Control::get_dimension_z()
{
    return dimension_z;
}
void Control::set_dimension_z( long long int dimension_z_)
{
    dimension_z=dimension_z_;
}
//-----------------------------------------
double Control::get_cell_a()
{
    return cell_a;
}
void Control::set_cell_a(double cell_a_)
{
    cell_a=cell_a_;
}
//-----------------------------------------
double Control::get_cell_b()
{
    return cell_b;
}
void Control::set_cell_b(double cell_b_)
{
    cell_b=cell_b_;
}
//-----------------------------------------
double Control::get_cell_c()
{
    return cell_c;
}
void Control::set_cell_c(double cell_c_)
{
    cell_c=cell_c_;
}
//-----------------------------------------
double Control::get_angle_gamma()
{
    return angle_gamma;
}
void Control::set_angle_gamma(double angle_gamma_)
{
    angle_gamma=angle_gamma_;
}
//-----------------------------------------
double Control::get_angle_alpha()
{
    return angle_alpha;
}
void Control::set_angle_alpha(double angle_alpha_)
{
    angle_alpha=angle_alpha_;
}
//-----------------------------------------
double Control::get_angle_beta()
{
    return angle_beta;
}
void Control::set_angle_beta(double angle_beta_)
{
    angle_beta=angle_beta_;
}
//-----------------------------------------
double Control::from_uvw_to_x(double u, double v, double w)
{
    return Util::from_uvw_to_x(u,v,w,cell_a,cell_b,cell_c,angle_beta,angle_gamma);
}
double Control::from_uvw_to_y(double v, double w)
{
    return Util::from_uvw_to_y(v,w, cell_b, cell_c, angle_alpha, angle_beta, angle_gamma);
}
double Control::from_uvw_to_z(double w)
{
    return Util::from_uvw_to_z(w,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}

double Control::from_xyz_to_u(double x, double y, double z)
{
    return Util::from_xyz_to_u(x,y,z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}
double Control::from_xyz_to_v(double y, double z)
{
    return Util::from_xyz_to_v(y,z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}
double Control::from_xyz_to_w(double z)
{
    return Util::from_xyz_to_w(z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}

//-----------------------------------------
double Control::get_box_boundx()
{
    return box_boundx;
}
void Control::set_box_boundx()
{
    box_boundx=Util::from_uvw_to_x(get_dimension_x(),0.0,0.0,get_cell_a(),get_cell_b(),get_cell_c(),get_angle_beta(),get_angle_gamma());
    //box_boundx=Util::from_uvw_to_x(get_dimension_x()*get_cell_a(),0.0,0.0,get_cell_a(),get_cell_b(),get_cell_c(),get_angle_beta(),get_angle_gamma()

}
//-----------------------------------------
double Control::get_box_boundx0()
{
    return box_boundx0;
}
void Control::set_box_boundx0()
{
    box_boundx0=Util::from_uvw_to_x(0.0,0.0,0.0,get_cell_a(),get_cell_b(),get_cell_c(),get_angle_beta(),get_angle_gamma());
}
//-----------------------------------------
double Control::get_box_boundxl0()
{
    return box_boundxl0;
}
void Control::set_box_boundxl0()
{
    box_boundxl0=get_box_boundx0()+min(min(0.0, get_tilt_factors_xy()), min(get_tilt_factors_xz(), get_tilt_factors_xy()+get_tilt_factors_xz()));
}
//-----------------------------------------
double Control::get_box_boundxl()
{
    return box_boundxl;
}
void Control::set_box_boundxl()
{
    box_boundxl=get_box_boundx()+max(max(0.0, get_tilt_factors_xy()), max(get_tilt_factors_xz(), get_tilt_factors_xy()+get_tilt_factors_xz()));
}
//-----------------------------------------

double Control::get_box_boundy()
{
    return box_boundy;
}
void Control::set_box_boundy()
{
    box_boundy=Util::from_uvw_to_y(get_dimension_y(),0.0,get_cell_b(),get_cell_c(),get_angle_alpha(),get_angle_beta(),get_angle_gamma());
}
//-----------------------------------------
double Control::get_box_boundy0()
{
    return box_boundy0;
}
void Control::set_box_boundy0()
{
    box_boundy0=Util::from_uvw_to_y(0.0,0.0,get_cell_b(),get_cell_c(),get_angle_alpha(),get_angle_beta(),get_angle_gamma());
}
//-----------------------------------------
double Control::get_box_boundyl0()
{
    return box_boundyl0;
}
void Control::set_box_boundyl0()
{
    box_boundyl0=get_box_boundy0()+min(0.0, get_tilt_factors_yz());
}
//-----------------------------------------
double Control::get_box_boundyl()
{
    return box_boundyl;
}
void Control::set_box_boundyl()
{
    box_boundyl=get_box_boundy()+max(0.0, get_tilt_factors_yz());
}
//-----------------------------------------

double Control::get_box_boundz()
{
    return box_boundz;
}
void Control::set_box_boundz()
{
    box_boundz=Util::from_uvw_to_z(get_dimension_z(), get_cell_a(), get_cell_b(), get_cell_c(), get_angle_alpha(),get_angle_beta(),get_angle_gamma());
}
//-----------------------------------------
double Control::get_box_boundz0()
{
    return box_boundz0;
}
void Control::set_box_boundz0()
{
    box_boundz0=Util::from_uvw_to_z(0.0, get_cell_a(), get_cell_b(), get_cell_c(), get_angle_alpha(),get_angle_beta(),get_angle_gamma());
}
//-----------------------------------------
double Control::get_box_boundzl0()
{
    return box_boundzl0;
}
void Control::set_box_boundzl0()
{
    box_boundzl0=get_box_boundz0();
}
//-----------------------------------------
double Control::get_box_boundzl()
{
    return box_boundzl;
}
void Control::set_box_boundzl()
{
    box_boundzl=get_box_boundz();
}
//-----------------------------------------

void Control::set_tilt_factors()
{
    tilt_factor_xy=get_cell_b()*get_dimension_y()*cos(get_angle_gamma()*PI/180.0);
    tilt_factor_xz=get_cell_c()*get_dimension_z()*cos(get_angle_beta()*PI/180.0);
    tilt_factor_yz=(get_cell_b()*get_dimension_y()*get_cell_c()*get_dimension_z()*cos(get_angle_alpha()*PI/180.0) - tilt_factor_xy*tilt_factor_xz)/(sqrt(pow(get_cell_b()*get_dimension_y(),2.0)-pow(tilt_factor_xy,2.0)));

}
//-----------------------------------------
void Control::set_tilt_bounds_factors() // los l son lo y hi respectivamente (en la pagina web)
{
    set_tilt_factors();

    set_box_boundx();
    set_box_boundx0();
    set_box_boundy();
    set_box_boundy0();
    set_box_boundz();
    set_box_boundz0();
    set_box_boundxl();
    set_box_boundxl0();
    set_box_boundyl();
    set_box_boundyl0();
    set_box_boundzl();
    set_box_boundzl0();
}

double Control::get_tilt_factors_xy()
{
  return tilt_factor_xy;
}
//-----------------------------------------
double Control::get_tilt_factors_xz()
{
  return tilt_factor_xz;
}
//-----------------------------------------
double Control::get_tilt_factors_yz()
{
  return tilt_factor_yz;
}
//-----------------------------------------
double Control::get_temperature()
{
    return temperature;
}
void Control::set_temperature(double temperature_)
{
    temperature=temperature_;
}
//-----------------------------------------
double Control::get_deltag()
{
    return deltag;
}
void Control::set_deltag(double deltag_)
{
    deltag=deltag_;
}
//-----------------------------------------
double Control::get_fd()
{
    return fd;
}
void Control::set_fd(double fd_)
{
    fd=fd_;
}
double Control::set_fd_from_T()
{
    fd=temperature*BOLTZCONSTANT/PLANCK;
    return fd;
}
//-----------------------------------------
double Control::get_fr()
{
    return fr;
}
void Control::set_fr(double fr_)
{
    fr=fr_;
}
//-----------------------------------------
double Control::get_Ed()
{
    return Ed;
}
void Control::set_Ed(double Ed_)
{
    Ed=Ed_;
}
//-----------------------------------------
double Control::get_Er()
{
    return Er;
}
void Control::set_Er(double Er_)
{
    Er=Er_;
}
//-----------------------------------------
double Control::get_time_modifier()
{
    return time_modifier;
}
void Control::set_time_modifier(double time_modifier_)
{
    time_modifier=time_modifier_;
}
//-----------------------------------------
double Control::get_scale_parameter_m1()
{
    return scale_parameter_m1;
}
void Control::set_scale_parameter_m1(double scale_parameter_m1_)
{
    scale_parameter_m1=scale_parameter_m1_;
}
//-----------------------------------------
double Control::get_scale_parameter_m2()
{
    return scale_parameter_m2;
}
void Control::set_scale_parameter_m2(double scale_parameter_m2_)
{
    scale_parameter_m2=scale_parameter_m2_;
}
//-----------------------------------------
double Control::get_time()
{
    return time;
}
void Control::set_time(double time_)
{
    time=time_;
}
void Control::increment_time_by(double increment)
{
    time+=increment;
}
//-----------------------------------------
double Control::get_targettime()
{
    return targettime;
}
void Control::set_targettime(double targettime_)
{
    targettime=targettime_;
}
//-----------------------------------------
double Control::get_resolution()
{
    return resolution;
}
void Control::set_resolution(double resolution_)
{
    resolution=resolution_;
    set_fd(get_fd()/resolution_);
    set_fr(get_fr()/resolution_);
}

//-----------------------------------------
long long int Control::add_list_event_definition(Event_definition * Event_definition_)
{
    list_event_definition.push_back(Event_definition_);
    return list_event_definition.size();
}

Event_definition * Control::get_event_definition_by_pos(long long int pos_)
{
    return list_event_definition.at(pos_);
}

vector <Event_definition*> Control::get_list_event_definition()
{
    return list_event_definition;
}
//-----------------------------------------

long long int Control::get_eventsaccepted()
{
    return eventsaccepted;
}
long long int Control::increaseby1_eventsaccepted()
{
    eventsaccepted++;
    return eventsaccepted;
}
//-----------------------------------------
double Control::change_distance_accuracy(double distance_accuracy_)
{
    Util::distance_accuracy=distance_accuracy_;
    return Util::distance_accuracy;
}

double Control::get_distance_accuracy()
{
    return Util::distance_accuracy;
}
//-----------------------------------------
double Control::get_unstuckaccuracy()
{
    return unstuckaccuracy;
}
void Control::set_unstuckaccuracy(double unstuckaccuracy_)
{
    unstuckaccuracy=unstuckaccuracy_;
}
//-----------------------------------------
 long long int Control::get_unstuckeverystep()
 {
    return unstuckeverystep;
 }
void Control::set_unstuckeverystep(long long int unstuckeverystep_)
{
    unstuckeverystep=unstuckeverystep_;
}
//-----------------------------------------
 long long int Control::get_unstuckalert()
 {
    return unstuckalert;
 }
void Control::set_unstuckalert(long long int unstuckalert_)
{
    unstuckalert=unstuckalert_;
}
long long int Control::increaseby1_unstuckalert()
{
    unstuckalert++;
    return unstuckalert;
}
long long int Control::restart_unstuckalert()
{
    unstuckalert=0;
    return unstuckalert;
}
//-----------------------------------------
long long int Control::get_filesdatanumber()
{
    return filesdatanumber;
}
void Control::set_filesdatanumber(long long int filesdatanumber_)
{
    filesdatanumber=filesdatanumber_;
}


//TRUE FALSE
//---------------
bool Control::get_estimatetime()
{
    return estimatetime;
}
void Control::set_estimatetime(bool estimatetime_)
{
    estimatetime=estimatetime_;
}
//-----------------------------------------
bool Control::get_randomize()
{
    return randomize;
}
void Control::set_randomize(bool randomize_)
{
    randomize=randomize_;
}
//-----------------------------------------
bool Control::get_dataanalisys()
{
    return dataanalisys;
}
void Control::set_dataanalisys(bool dataanalisys_)
{
    dataanalisys=dataanalisys_;
}
//-----------------------------------------
bool Control::get_layerxanalysis()
{
    return layerxanalysis;
}
void Control::set_layerxanalysis(bool layerxanalysis_)
{
    layerxanalysis=layerxanalysis_;
}
//-----------------------------------------
bool Control::get_layeryanalysis()
{
    return layeryanalysis;
}
void Control::set_layeryanalysis(bool layeryanalysis_)
{
    layeryanalysis=layeryanalysis_;
}
//-----------------------------------------
bool Control::get_layerzanalysis()
{
    return layerzanalysis;
}
void Control::set_layerzanalysis(bool layerzanalysis_)
{
    layerzanalysis=layerzanalysis_;
}
//-----------------------------------------
bool Control::get_coordinationanalysis()
{
    return coordinationanalysis;
}
void Control::set_coordinationanalysis(bool coordinationanalysis_)
{
    coordinationanalysis=coordinationanalysis_;
}
//-----------------------------------------
bool Control::get_boxframe()
{
    return boxframe;
}
void Control::set_boxframe(bool boxframe_)
{
    boxframe=boxframe_;
}
//-----------------------------------------
bool Control::get_surfaceframe()
{
    return surfaceframe;
}
void Control::set_surfaceframe(bool surfaceframe_)
{
    surfaceframe=surfaceframe_;
}
//-----------------------------------------
bool Control::get_initialdata()  //esto no se usa
{
    return initialdata;
}
void Control::set_initialdata(bool initialdata_) //esto no se usa
{
    initialdata=initialdata_;
}
//-----------------------------------------
string Control::get_workname()
{
    return workname;
}
void Control::set_workname(string workname_)
{
    workname=workname_;
}
//-----------------------------------------

double Control::get_unstuck_compare_time()
{
    //the higher unstuckaccuracy value, the smaller the comparation time, and the simulation accuracy is higher
    return get_targettime()/((double)get_unstuckaccuracy()*(double)(get_filesdatanumber()));
}


double Control::estimate_time()
{
    double estimatedtime=3.0/( get_fd() *  exp(-get_Ed()*COORDINATIONDISLOCA));
    set_targettime(estimatedtime);
    return estimatedtime;
}

//----Auxiliar printing parameters----//
 long long int Control::get_step()
{
    return step;
}

 long long int Control::increaseby1_step()
{
    step++;
    return step;
}

//-----------------------------------------
bool Control::get_grain()
{
    return grain;
}
void Control::set_grain(bool grain_)
{
    grain=grain_;
}
//-----------------------------------------
bool Control::get_plane_x()
{
    return planex;
}
void Control::set_plane_x(bool plane_x)
{
    planex=plane_x;
}
//-----------------------------------------
bool Control::get_plane_y()
{
    return planey;
}
void Control::set_plane_y(bool plane_y)
{
    planey=plane_y;
}
//-----------------------------------------
bool Control::get_plane_z()
{
    return planez;
}
void Control::set_plane_z(bool plane_z)
{
    planez=plane_z;
}
//-----------------------------------------

double Control::get_gyradiousini()
{
    return gyradiousini;
}
void Control::set_gyradiousini(double gyradiousini_)
{
    gyradiousini=gyradiousini_;
}
//-----------------------------------------
 long long int Control::get_maxcoordini()
{
    return maxcoordini;
}
void Control::set_maxcoordini( long long int maxcoordini_)
{
    maxcoordini=maxcoordini_;
}
//-----------------------------------------
 long long int Control::get_totalnumberpositions()
{
    return totalnumberpositions;
}
long long int Control::increaseby1_totalnumberpositions()
{
    totalnumberpositions++;
    return totalnumberpositions;
}
long long int Control::decreaseby1_totalnumberpositions()
{
    totalnumberpositions--;
    return totalnumberpositions;
}
//-----------------------------------------
long long int Control::get_datasteps()
{
    return datasteps;
}
void Control::set_datastepts(long long int datasteps_)
{
    datasteps=datasteps_;
}
//-----------------------------------------
long long int Control::get_layersteps()
{
    return layersteps;
}
void Control::set_layersteps(long long int layersteps_)
{
    layersteps=layersteps_;
}
//-----------------------------------------
long long int Control::get_boxsteps()
{
    return boxsteps;
}
void Control::set_boxsteps(long long int boxsteps_)
{
    boxsteps=boxsteps_;
}
//-----------------------------------------
long long int Control::get_surfacesteps()
{
    return surfacesteps;
}
void Control::set_surfacesteps(long long int surfacesteps_)
{
    surfacesteps=surfacesteps_;
}
//-----------------------------------------
long long int Control::get_meansteps()
{
    return meansteps;
}
void Control::set_meansteps(long long int meansteps_)
{
    meansteps=meansteps_;
}
//-----------------------------------------
bool Control::get_xyzfile()
{
    return xyzfile;
}
void Control::set_xyzfile(bool xyzfile_)
{
    xyzfile=xyzfile_;
}
//-----------------------------------------
string Control::get_pathtoxyzfile()
{
    return pathtoxyzfile;
}
void Control::set_pathtoxyzfile(string pathtoxyzfile_)
{
    pathtoxyzfile=pathtoxyzfile_;
}
//-----------------------------------------
bool Control::get_kimerafile()
{
    return kimerafile;
}
void Control::set_kimerafile(bool kimerafile_)
{
    kimerafile=kimerafile_;
}
//-----------------------------------------
string Control::get_pathtokimerafile()
{
    return pathtokimerafile;
}
void Control::set_pathtokimerafile(string pathtokimerafile_)
{
    pathtokimerafile=pathtokimerafile_;
}
//-----------------------------------------
bool Control::get_printinitialkimerafile()
{
    return printinitialkimerafile;
}
void Control::set_printinitialkimerafile(bool printinitialkimerafile_)
{
    printinitialkimerafile=printinitialkimerafile_;
}
//-----------------------------------------
bool Control::get_printfinalkimerafile()
{
    return printfinalkimerafile;
}
void Control::set_printfinalkimerafile(bool printfinalkimerafile_)
{
    printfinalkimerafile=printfinalkimerafile_;
}
//-----------------------------------------
bool Control::get_print_by_steps()
{
    return print_by_steps;
}

void Control::set_print_by_steps(bool print_by_steps_)
{
    print_by_steps=print_by_steps_;
}
//-----------------------------------------
long long int Control::get_targetsteps()
{
    return targetsteps;
}

void Control::set_targetsteps(long long int targetsteps_)
{
    targetsteps=targetsteps_;
}

bool Control::get_parallelize()
{
    return parallelize;
}
void Control::set_parallelize(bool parallelize_)
{
    parallelize=parallelize_;
}
//-----------------------------------------
long long int Control::get_parallelizenumber()
{
    return parallelizenumber;
}
void Control::set_parallelizenumber(long long int parallelizenumber_)
{
    parallelizenumber=parallelizenumber_;
}
//-----------------------------------------
bool Control::get_linealsearch()
{
    return linealsearch;
}
void Control::set_linealsearch(bool linealsearch_)
{
    linealsearch=linealsearch_;
}
//-----------------------------------------
bool Control::get_meandiscoord()
{
    return meandiscoord;
}
void Control::set_meandiscoord(bool meandiscoord_)
{
    meandiscoord=meandiscoord_;
}
//-----------------------------------------
bool Control::get_affected()
{
    return affected;
}
void Control::set_affected(bool affected_)
{
    affected=affected_;
}
//-----------------------------------------
bool Control::get_linked()
{
    return linked;
}
void Control::set_linked(bool linked_)
{
    linked=linked_;
}
//-----------------------------------------
bool Control::get_seedboxbool()
{
    return seedboxbool;
}
void Control::set_seedboxbool(bool seedboxbool_)
{
    seedboxbool=seedboxbool_;
}
//-----------------------------------------
long long int Control::get_seedbox()
{
    return seedbox;
}
void Control::set_seedbox(long long int seedbox_)
{
    seedbox=seedbox_;
}
//-----------------------------------------
bool Control::get_seedsimbool()
{
    return seedsimbool;
}
void Control::set_seedsimbool(bool seedsimbool_)
{
    seedsimbool=seedsimbool_;
}
//-----------------------------------------
long long int Control::get_seedsim()
{
    return seedsim;
}
void Control::set_seedsim(long long int seedsim_)
{
    seedsim=seedsim_;
}
//-----------------------------------------
 long long int Control::get_id_linked()
 {
     return id_linked;
 }
 void Control::set_id_linked(long long int id_)
 {
     id_linked=id_;
 }
//-----------------------------------------
