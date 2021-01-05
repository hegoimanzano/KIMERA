#ifndef MODELS_H
#define MODELS_H

#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

class Models
{
    public:
        Models();
        virtual ~Models();

        static double model42_dissolution(double fd, long long int neighbour, double Ed );
        static double model42_reforming(double fr, long long int neighbour, double Er, double deltaG);

        static double model42_dissolution(double fd, vector <long long int> Nneighbour, vector <double> Eds);
        static double model42_reforming(double fr, vector <long long int> Nneighbour,  vector <double> Eps, double deltaG);

        static double model42_dissolution(double fd, vector <long long int> Nneighbour, vector <double> Eds , vector<vector<double>> Eds_direct);
        static double model42_reforming(double fr, vector <long long int> Nneighbour,  vector <double> Eps, vector<vector<double>> Eps_direct, double deltaG);

        static double model50_dissolution(double fd, double fr, vector <long long int> Nneighbour, vector <double> Eds , vector <double> Eps,  vector<vector<double>> Eds_direct, vector<vector<double>> Eps_direct ,double  deltaG);  // poner llegado el momento, no haria falta mas que los eventos de disolucion

    protected:

    private:
};

#endif // MODELS_H
