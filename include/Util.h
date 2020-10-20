#ifndef UTIL_H
#define UTIL_H

#define EPSILON 1.0e-5    //isequal
#define EPSILON2 0.0001  //distance  //quartz distancia oxigeno oxigeno 2.585 +- 0.01     0.0001 antes
#define EPSILON3 1.0e-30
#define PI 3.14159265359

#include <iostream>
#include <math.h>



using namespace std;


class Util
{



    public:
        Util();
        virtual ~Util();

        static double distance_accuracy;

        static bool isEqualhigh(double a,double b);

        static bool isEqual(double a,double b);

        static bool isEqualDistances(double a,double b);

        static double get_omega(double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);

        static double from_uvw_to_x(double u, double v, double w, double cell_a, double cell_b, double cell_c, double cell_beta, double cell_gamma);

        static double from_uvw_to_y(double v, double w, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);

        static double from_uvw_to_z(double w, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);

        static double from_xyz_to_u(double x, double y, double z, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);

        static double from_xyz_to_v( double y, double z,  double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);

        static double from_xyz_to_w(double z, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma);


        static double distance_point_to_plane(double x_point, double y_point, double z_point, double A_plane, double B_plane, double C_plane, double D_plane);



};

#endif // UTIL_H
