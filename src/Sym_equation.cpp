#include "Sym_equation.h"


Sym_equation::Sym_equation(string expresion_)
{
    expresion=expresion_;
    read_expresions_set_parameters(expresion_);
}

Sym_equation::~Sym_equation()
{
    //dtor
}

long long int Sym_equation::read_expresions_set_parameters(string expresion_)
{
    //A+ax+B+by+C+cz
    A=0.0; a=1.0; B=0.0; b=1.0; C=0.0; c=1.0;
    if (expresion_.compare("'x,y,z'")==0)
    {
        return 0;
    }

    long long int index = expresion_.find('x');
    long long int indey = expresion_.find('y');
    long long int indez = expresion_.find('z');

    //buscamos + o - antes que la x
    string prex = expresion_.substr (index-1,index);


    if (prex.compare("-")==0)   a=-1.0;

    if (prex.compare("+")==0 || prex.compare("-")==0)
    {
        string prexx = expresion_.substr (index-2,index-1);
        string prexxxx = expresion_.substr (index-4,index-3);
        double den=stod(prexx);
        double num=stod(prexxxx);
        A=num/den;
    }

    string prey = expresion_.substr (indey-1,indey);

    if (prey.compare("-")==0)   b=-1.0;

    if (prey.compare("+")==0 || prey.compare("-")==0)
    {
        string preyy = expresion_.substr (indey-2,indey-1);
        string preyyyy = expresion_.substr (indey-4,indey-3);
        double den=stod(preyy);
        double num=stod(preyyyy);
        B=num/den;
    }

    string prez = expresion_.substr (indez-1,indez);

    if (prez.compare("-")==0)   c=-1.0;

    if (prez.compare("+")==0 || prez.compare("-")==0)
    {
        string prezz = expresion_.substr (indez-2,indez-1);
        string prezzzz = expresion_.substr (indez-4,indez-3);
        double den=stod(prezz);
        double num=stod(prezzzz);
        C=num/den;
    }
        return 0;
}



    /**
    string prex = expresion_.substr (0,index);
    long long int indexplusminusx=prex.find("+");
    //si no se encuentra, a=1.0
    if (indexplusminusx==npos) a=1.0



         string str = "0,07";
     long long int index = str.find(',');
     str.replace(index, index+1, '.');

     double number = stod(str);

    A=0.0;
    B=0.0;
    C=0.0;

    a=0.0;
    b=0.0;
    c=0.0;
*/

Position Sym_equation::get_sym_pos(Position * position_)
{
    double x=position_->get_x();
    double y=position_->get_y();
    double z=position_->get_z();
    long long int newid=NEWPOSID;
    vector<string> types=position_->get_types();
    vector<double> probs=position_->get_probs();


    double newx=A+a*x;
    double newy=B+b*y;
    double newz=C+c*z;
    string firsttype=types.at(0);
    double firstprob=probs.at(0);

    Position auxposition=Position(newx,newy,newz,newid,firsttype,firstprob);

    //add the rest of the types and probs
    long long int typessize=types.size();
    if (typessize>1)
    {
        for ( long long int i=1; i<(long long int)types.size(); i++)
        {
            auxposition.add_type_prob(types.at(i),probs.at(i));
        }
    }


    return auxposition;
}
