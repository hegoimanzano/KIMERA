#include "Models.h"

Models::Models()
{
    //ctor
}

Models::~Models()
{
    //dtor
}

double Models::model42_dissolution(double fd, long long int neighbour, double Ed )
{
    double rate= fd*exp(-1.0*(double)neighbour*Ed);
    return rate;
}

double Models::model42_reforming(double fr, long long int neighbour, double Er, double deltaG)
{
    double rate= fr *exp(-1.0*(Er*(double)neighbour-deltaG));
    return rate;
}

double Models::model42_dissolution(double fd, vector <long long int> Nneighbour, vector <double> Eds )
{
    long long int sizeN=  Nneighbour.size();
    double exponent=0.0;

    for (long long int i = 0; i< sizeN; i++)
    {
        exponent+=(double)Nneighbour.at(i) * Eds.at(i);
    }

    double rate= fd*exp(-1.0*exponent);


    return rate;
}

double Models::model42_reforming(double fr, vector <long long int> Nneighbour,  vector <double> Eps, double deltaG)
{
    long long int sizeN= Nneighbour.size();
    double exponent=0.0;

    for (long long int i = 0; i< sizeN; i++)
    {
        exponent+=(double)Nneighbour.at(i) * Eps.at(i);
    }

    double rate= fr *exp(-1.0*(exponent-deltaG));
    return rate;
}

double Models::model42_dissolution(double fd, vector <long long int> Nneighbour, vector <double> Eds , vector<vector<double>> Eds_direct)
{
    long long int sizeN=  Nneighbour.size();
    double exponent=0.0;
    long long int size_Eds_direct=0;
    long long int j=0;


    for (long long int i = 0; i< sizeN; i++)
    {
        size_Eds_direct=Eds_direct.at(i).size();
        if (size_Eds_direct==0)  //si el tamaño es 0, si se considera el modelo lineal
        {
            exponent+=(double)Nneighbour.at(i) * Eds.at(i);
        }
        else   //de lo contrario se coge la posicion j de i, donde j es la cantidad de vecinos-1  (Nneighbour.at(i)-1)    puede pasar que el numero de vecinos sea mayor que el tamaños que hemos definido
        {

            j = Nneighbour.at(i)-1;
            if (Nneighbour.at(i)==0) {exponent+= 0.0; } //luego no se incrementa el tiempo
            else
            {
                if(j>=size_Eds_direct){cout << "Error! the direct definition do not cover the amount of neighbours, check the neighbour list"<<endl; exit(1); }


                exponent+=Eds_direct.at(i).at(j); //El usuario debe tener en cuenta que se suman los anteriores

            }
        }

    }



    double rate= fd*exp(-1.0*exponent);


    return rate;
}

double Models::model42_reforming(double fr, vector <long long int> Nneighbour,  vector <double> Eps, vector<vector<double>> Eps_direct, double deltaG)
{
    long long int sizeN= Nneighbour.size();
    double exponent=0.0;
    long long int size_Eps_direct=0;
    long long int j=0;

    for (long long int i = 0; i< sizeN; i++)
    {
        size_Eps_direct=Eps_direct.at(i).size();
        if (size_Eps_direct==0)  //si el tamaño es 0, si se considera el modelo lineal
        {
            exponent+=(double)Nneighbour.at(i) * Eps.at(i);
        }
        else   //de lo contrario se coge la posicion j de i, donde j es la cantidad de vecinos-1  (Nneighbour.at(i)-1)    puede pasar que el numero de vecinos sea mayor que el tamaños que hemos definido
        {
            j = Nneighbour.at(i)-1;
            if (Nneighbour.at(i)==0) {exponent+= 0.0; } //luego no se incrementa el tiempo
            else
            {
                if(j>=size_Eps_direct){cout << "Error! the direct definition do not cover the amount of neighbours, check the neighbour list"<<endl; exit(1);  }

                exponent+=Eps_direct.at(i).at(j); //El usuario debe tener en cuenta que se suman los anteriores
            }
        }
    }


    double rate= fr *exp(-1.0*(exponent-deltaG));
    return rate;
}


double Models::model50_dissolution(double fd, double fr, vector <long long int> Nneighbour, vector <double> Eds , vector <double> Eps,  vector<vector<double>> Eds_direct, vector<vector<double>> Eps_direct ,double  deltaG)
{


    long long int sizeN=  Nneighbour.size();
    double exponentd=0.0;
    long long int size_Eds_direct=0;
    long long int j=0;

    for (long long int i = 0; i< sizeN; i++)
    {
        size_Eds_direct=Eds_direct.at(i).size();
        if (size_Eds_direct==0)  //si el tamaño es 0, si se considera el modelo lineal
        {
            exponentd+=(double)Nneighbour.at(i) * Eds.at(i);
        }
        else   //de lo contrario se coge la posicion j de i, donde j es la cantidad de vecinos-1  (Nneighbour.at(i)-1)    puede pasar que el numero de vecinos sea mayor que el tamaños que hemos definido
        {

            j = Nneighbour.at(i)-1;
            if (Nneighbour.at(i)==0) {exponentd+= 0.0; } //luego no se incrementa el tiempo
            else
            {
                if(j>=size_Eds_direct){cout << "Error! the direct definition do not cover the amount of neighbours, check the neighbour list. Amount of neighbours: "<< Nneighbour.at(i) <<endl; exit(1); }


                exponentd+=Eds_direct.at(i).at(j); //El usuario debe tener en cuenta que se suman los anteriores

            }
        }

    }


    sizeN= Nneighbour.size();
    double exponentp=0.0;
    long long int size_Eps_direct=0;
    j=0;

    for (long long int i = 0; i< sizeN; i++)
    {
        size_Eps_direct=Eps_direct.at(i).size();
        if (size_Eps_direct==0)  //si el tamaño es 0, si se considera el modelo lineal
        {
            exponentp+=(double)Nneighbour.at(i) * Eps.at(i);
        }
        else   //de lo contrario se coge la posicion j de i, donde j es la cantidad de vecinos-1  (Nneighbour.at(i)-1)    puede pasar que el numero de vecinos sea mayor que el tamaños que hemos definido
        {
            j = Nneighbour.at(i)-1;
            if (Nneighbour.at(i)==0) {exponentp+= 0.0; } //luego no se incrementa el tiempo
            else
            {
                if(j>=size_Eps_direct){cout << "Error! the direct definition do not cover the amount of neighbours, check the neighbour list. Amount of neighbours: "<< Nneighbour.at(i) <<endl; exit(1);  }

                exponentp+=Eps_direct.at(i).at(j); //El usuario debe tener en cuenta que se suman los anteriores
            }
        }
    }



    double rate= fd*exp(-exponentd) *exp(-(-log(0.001))*exp(-exponentp+deltaG)*fr/(exp(-exponentd)*fd) );
    return rate;




}
