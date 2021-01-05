#include "Util.h"


Util::Util()
{

}

Util::~Util()
{
    //dtor
}

double Util::distance_accuracy=EPSILON2;


bool Util::isEqual(double a,double b)
{
    if (a-b< EPSILON && a-b> -EPSILON)  return true;
    else return false;
}

bool Util::isEqualhigh(double a,double b)
{
    if (a-b< EPSILON3 && a-b> -EPSILON3)  return true;
    else return false;
}


bool Util::isEqualDistances(double a,double b)
{
    if (a-b< distance_accuracy && a-b> -distance_accuracy)  return true;
    else return false;
}

double Util::get_omega(double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    return cell_a*cell_b*cell_c*sqrt(1.0-pow(cos(cell_alpha*PI/180.0),2.0)-pow(cos(cell_beta*PI/180.0),2.0)-pow(cos(cell_gamma*PI/180.0),2.0)
            + 2.0*cos(cell_alpha*PI/180.0)*cos(cell_beta*PI/180.0)*cos(cell_gamma*PI/180.0));
}

double Util::from_uvw_to_x(double u, double v, double w, double cell_a, double cell_b, double cell_c, double cell_beta, double cell_gamma)
{
    double x = cell_a*u +
               cell_b*cos(cell_gamma*PI/180.0) * v +
               cell_c*cos(cell_beta*PI/180.0) * w;

    return x;
}

double Util::from_uvw_to_y(double v, double w, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    double y =  cell_b*sin(cell_gamma*PI/180.0) * v +
                cell_c*(cos(cell_alpha*PI/180.0) - cos(cell_beta*PI/180.0)*cos(cell_gamma*PI/180.0))/sin(cell_gamma*PI/180.0) * w;
    return y;
}

double Util::from_uvw_to_z(double w, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    double omega= get_omega(cell_a,cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);
    double z =  omega/(cell_a*cell_b*sin(cell_gamma*PI/180.0)) *w;
    return z;
}



double Util::from_xyz_to_u(double x, double y, double z, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    double omega= get_omega(cell_a,cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);
    double u=   1.0/cell_a * x
                - cos(cell_gamma*PI/180.0)/(cell_a*sin(cell_gamma*PI/180.0)) * y
                + cell_b*cell_c*(cos(cell_alpha*PI/180.0)*cos(cell_gamma*PI/180.0)-cos(cell_beta*PI/180.0))/(omega*sin(cell_gamma*PI/180.0)) * z;
    return u;
}

double Util::from_xyz_to_v( double y, double z,  double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    double omega= get_omega(cell_a,cell_b, cell_c, cell_alpha, cell_beta, cell_gamma);
    double v=   1.0/(cell_b*sin(cell_gamma*PI/180.0)) * y
                + cell_a*cell_c*(cos(cell_beta*PI/180.0)*cos(cell_gamma*PI/180.0)-cos(cell_alpha*PI/180.0))/(omega*sin(cell_gamma*PI/180.0)) * z;
    return v;
}

double Util::from_xyz_to_w(double z, double cell_a, double cell_b, double cell_c, double cell_alpha, double cell_beta, double cell_gamma)
{
    double omega=get_omega(cell_a,cell_b, cell_c, cell_alpha,cell_beta, cell_gamma);
    double w=   cell_a * cell_b * sin(cell_gamma*PI/180.0) / omega * z;
    return w;
}

double Util::distance_point_to_plane(double x_point, double y_point, double z_point, double A_plane, double B_plane, double C_plane, double D_plane)
{
    double result=(fabs(A_plane*x_point+B_plane*y_point+C_plane*z_point+D_plane)/(sqrt(A_plane*A_plane+B_plane*B_plane+C_plane*C_plane)));
    return result;
}
