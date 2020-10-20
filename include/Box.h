#ifndef BOX_H
#define BOX_H


#include <omp.h>
#include "Event_definition.h"
#include "Record.h"
#include "Cell.h"
#include "Position.h"
#include "RandomGenerator.h"
#include "Atom.h"
#include "Util.h"
#include <vector>
#include <math.h>
#include <iostream>

#define KSNEIGH 6





using namespace std;

class Box
{
    public:
        Box( long long int dimension_x_,  long long int dimension_y_,  long long int dimension_z_,
            double cell_a_, double cell_b_, double cell_c,
            double angle_gamma_, double angle_beta_, double angle_alpha_);
        Box();

        virtual ~Box();

        double get_dimension_x();
        double get_dimension_y();
        double get_dimension_z();



        vector<Atom*> get_box_atoms();    // returns a pointer to a vector which constains a list of atoms.
        void set_box_atoms();

        void add_atom(Atom *atom);
        Atom * add_atom( long long int id, string atom_type, long long int type ,double x, double y, double z);
        long long int rm_atom(long long int id);
        long long int rm_atom_complex(long long int id);

        void set_neighbour(bool periodicity_x,bool periodicity_y,bool periodicity_z);
        void set_neighbour(double distance_);
        //void set_neighbour(vector * list_neighbours); // set the neighbour from a list |id atom

        void mv_atom(Atom *atom, double x, double y, double z);

        //void define_surface();

        //SURFACE METHODS
        long long int add_surface_atom(Atom *atom);
        long long int rm_surface_atom(long long int id);

        vector<Atom*> get_surface_atoms();
        void set_initial_surface_atoms();

        long long int get_surface_size();

        void add_new_atom(long long int dimension_a,long long int dimension_b, long long int dimension_c, long long int type);


        //IRREGULARITIES
        //------------------------------------------------------------------------
        long long int remove_sym_position(long long int dimension_a, long long int dimension_b, long long int dimension_c);
        long long int add_sym_position(long long int dimension_a, long long int dimension_b, long long int dimension_c,long long int type);
        long long int add_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez,long long int type);
        long long int rm_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);
        long long int add_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez,long long int type);
        long long int rm_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int add_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey,long long int type);
        long long int rm_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int add_xy_heli_dislocation( long long int dimension_a, long long int dimension_b, long long int radiousdisloca);
        long long int add_xz_heli_dislocation( long long int dimension_a, long long int dimension_c, long long int radiousdisloca);
        long long int add_yz_heli_dislocation( long long int dimension_b, long long int dimension_c, long long int radiousdisloca);
        long long int define_insoluble_atom(long long int dimension_a, long long int dimension_b, long long int dimension_c);
        long long int define_insoluble_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int define_insoluble_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int define_insoluble_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);

        //------------------------------------------------------------------------


        //Good IRREGULARITIES
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        long long int change_type_sym_cell(long long int dim_a, long long int dim_b, long long int dim_c,long long int type);
        long long int change_type_cell_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, long long int type);
        long long int change_type_cell_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, long long int type);
        long long int change_type_cell_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, long long int type);

        long long int add_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);//x,z postion at z=0
        long long int add_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca);//x,z postion at z=0
        long long int add_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0
        long long int add_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0

        //Definimos planos
        //Ax+By+Cz+D=0
        long long int change_type_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_, long long int type);

        //definimos gap
        long long int change_type_gap(double x_gap,double y_gap, double z_gap,double radious_gap,long long int type);

        long long int change_type_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere,long long int type);
        long long int change_type_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, long long int type);
        long long int change_type_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,long long int type);
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //IRREGULARITIES SELECTIVAS
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        long long int change_type_sym_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c,long long int type);
        long long int change_type_cell_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, long long int type);
        long long int change_type_cell_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, long long int type);
        long long int change_type_cell_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, long long int type);

        long long int add_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);//x,z postion at z=0
        long long int add_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca);//x,z postion at z=0
        long long int add_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0
        long long int add_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0

        //Definimos planos
        //Ax+By+Cz+D=0
        long long int change_type_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_, long long int type);

        //definimos gap
        long long int change_type_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap,long long int type);

        long long int change_type_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere,long long int type);
        long long int change_type_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, long long int type);
        long long int change_type_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,long long int type);
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------

        //--------------------------REMOVE_FROM_SURFACE--------------------------



        long long int remove_from_surface_cell(long long int dim_a, long long int dim_b, long long int dim_c);
        long long int remove_from_surface_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);
        long long int remove_from_surface_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int remove_from_surface_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int remove_from_surface_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int remove_from_surface_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca);
        long long int remove_from_surface_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);
        long long int remove_from_surface_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca);
        long long int remove_from_surface_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0
        long long int remove_from_surface_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca);
        long long int remove_from_surface_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_);

        long long int remove_from_surface_gap(double x_gap,double y_gap, double z_gap,double radious_gap);
        long long int remove_from_surface_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere);
        long long int remove_from_surface_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z);
        long long int remove_from_surface_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K);
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        long long int add_to_surface_cell(long long int dim_a, long long int dim_b, long long int dim_c);
        long long int add_to_surface_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);
        long long int add_to_surface_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int add_to_surface_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int add_to_surface_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_to_surface_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca);
        long long int add_to_surface_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);
        long long int add_to_surface_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca);
        long long int add_to_surface_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);//x,z postion at z=0
        long long int add_to_surface_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca);
        long long int add_to_surface_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_);

        long long int add_to_surface_gap(double x_gap,double y_gap, double z_gap,double radious_gap);
        long long int add_to_surface_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere);
        long long int add_to_surface_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z);
        long long int add_to_surface_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K);



        //REVOME_FROM_SURFACE SELECTIVAS
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        long long int remove_from_surface_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c);
        long long int remove_from_surface_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);
        long long int remove_from_surface_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int remove_from_surface_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int remove_from_surface_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int remove_from_surface_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca);
        long long int remove_from_surface_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);
        long long int remove_from_surface_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca);
        long long int remove_from_surface_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);
        long long int remove_from_surface_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca);

        long long int remove_from_surface_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_);
        long long int remove_from_surface_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap);
        long long int remove_from_surface_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere);
        long long int remove_from_surface_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z);
        long long int remove_from_surface_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K);


        //ADD_TO_SURFACE SELECTIVAS
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        long long int add_to_surface_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c);
        long long int add_to_surface_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez);
        long long int add_to_surface_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez);
        long long int add_to_surface_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey);
        long long int add_to_surface_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca);//x,y postion at z=0
        long long int add_to_surface_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca);
        long long int add_to_surface_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca);
        long long int add_to_surface_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca);
        long long int add_to_surface_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca);
        long long int add_to_surface_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca);

        long long int add_to_surface_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_);
        long long int add_to_surface_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap);
        long long int add_to_surface_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere);
        long long int add_to_surface_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z);
        long long int add_to_surface_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K);


          //--------------------------CHANGE_ELEMENT--------------------------

        long long int change_element_cell(long long int dim_a, long long int dim_b, long long int dim_c,string new_atom_type);
        long long int change_element_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, string new_atom_type);
        long long int change_element_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, string new_atom_type);
        long long int change_element_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, string new_atom_type);

        long long int change_element_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca, string new_atom_type);//x,y postion at z=0
        long long int change_element_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca, string new_atom_type);//x,y postion at z=0
        long long int change_element_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca, string new_atom_type);//x,z postion at z=0
        long long int change_element_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca, string new_atom_type);//x,z postion at z=0
        long long int change_element_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca, string new_atom_type);//x,z postion at z=0
        long long int change_element_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca, string new_atom_type);//x,z postion at z=0

        //Definimos planos
        //Ax+By+Cz+D=0
        long long int change_element_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_, string new_atom_type);

        long long int change_element_gap(double x_gap,double y_gap, double z_gap,double radious_gap, string new_atom_type);

        long long int change_element_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere, string new_atom_type);
        long long int change_element_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, string new_atom_type);
        long long int change_element_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,string new_atom_type);


         //--------------------------CHANGE_ELEMENT_TO_TYPE--------------------------

        long long int change_element_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c, string new_atom_type);
        long long int change_element_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, string new_atom_type);
        long long int change_element_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, string new_atom_type);
        long long int change_element_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, string new_atom_type);

        long long int change_element_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca,string new_atom_type);//x,y postion at z=0
        long long int change_element_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca,string new_atom_type);//x,y postion at z=0
        long long int change_element_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca,string new_atom_type);//x,z postion at z=0
        long long int change_element_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca,string new_atom_type);//x,z postion at z=0
        long long int change_element_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca,string new_atom_type);//x,z postion at z=0
        long long int change_element_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca,string new_atom_type);//x,z postion at z=0

        //Definimos planos
        //Ax+By+Cz+D=0
        long long int change_element_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_, string new_atom_type);

        long long int change_element_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap,string new_atom_type);

        long long int change_element_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere,string new_atom_type);
        long long int change_element_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, string new_atom_type);
        long long int change_element_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K, string new_atom_type);


        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------



        long long int set_mass(string type, double mass);

        long long int clean_gaps_before_start();
        long long int clean_gaps_before_start_complex();
        long long int clean_gaps_before_start_complex_from_kimerafile();

        void update_surface(Atom * dissolved_atom); //introuduce los  vecinos del atomo disuelto si no estan ya introducidos //QUITAR

        long long int dissolve_atom(Atom * dissolvedatom);

        bool check_in_surface(long long int id);

        Atom* select_atom_by_id(long long int id);
        Atom* select_atom_by_position(long long int dimension_a, long long int dimension_b, long long int dimension_c);
        double from_dimension_to_x(long long int dimension_a,long long int dimension_b,long long int dimension_c);
        double from_dimension_to_y(long long int dimension_b,long long int dimension_c);
        double from_dimension_to_z(long long int dimension_c);

        long long int from_possition_to_dimension_a(double x_, double y_, double z_);
        long long int from_possition_to_dimension_b(double y_, double z_);
        long long int from_possition_to_dimension_c(double z_);

        vector<Atom*> get_box_removed_atoms();

        long long int create_box_by_cell(Cell * first_cell,long long int seed);
        long long int create_box_by_cell(Cell * first_cell); //le metemos la cell con las posiciones y eso y las duplicamos segun las dimensiones
                                                    //tenemos en cuenta las occupancies para meter el átomo que sea

        vector<Cell*> get_box_cells();
        Cell * get_cell_by_id(long long int cell_id);
        //TODO

        long long int add_cell_to_box(Cell * cell_);
        Cell * add_cell_to_box(long long int id_cell);


        double from_uvw_to_x(double u, double v, double w);
        double from_uvw_to_y(double v, double w);
        double from_uvw_to_z(double w);
        double from_xyz_to_u(double x, double y, double z);
        double from_xyz_to_v(double y, double z);
        double from_xyz_to_w(double z);

        long long int set_cell_neighbours_26(bool periodicity_x,bool periodicity_y,bool periodicity_z);
        long long int set_atom_neighbourhood(string targettype_, vector<string> neigh_type,vector<double> distances_); //se construye al final el record de los vecinos
        long long int set_atom_neighbourhood_linked(string targettype_, vector<string> neigh_type,vector<double> distances_, vector<Linked_neighbour> linked_neighbours);
        long long int set_atom_affected(string targettype_, vector<string> affected_type,vector<double> affected_distances_); //quizas en un futuro unir metodos

        long long int set_initial_surface_atoms(vector<Event_definition*> defined_events);

        Cell* select_cell_by_position(long long int dimension_a, long long int dimension_b, long long int dimension_c);

    protected:

    private:

         long long int dimension_x;
         long long int dimension_y;
         long long int dimension_z;

        double cell_a;
        double cell_b;
        double cell_c;
        double angle_alpha;
        double angle_beta;
        double angle_gamma;


        vector<Atom*> box_atoms;
        vector<Atom*> box_removed_atoms;
        vector<Atom*> box_removed_atoms_before_start;
        vector<Atom*> surface_atoms;

        vector<Cell*> box_cells;   //Borrar las cells al final




};

#endif // BOX_H
