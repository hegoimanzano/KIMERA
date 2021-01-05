#ifndef SYM_EQUATION_H
#define SYM_EQUATION_H

#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include "Position.h"

#define NEWPOSID -1


using namespace std;

class Sym_equation
{
    public:
        Sym_equation(string expresion_);
        virtual ~Sym_equation();

        long long int read_expresions_set_parameters(string expresion_);
        Position get_sym_pos(Position* position_);


    protected:

    private:

        string expresion; //A+ax+B+by+C+cz

        double a;
        double b;
        double c;
        double A;
        double B;
        double C;
};

#endif // SYM_EQUATION_H
