#include "Box.h"



Box::Box( long long int dimension_x_,  long long int dimension_y_,   long long int dimension_z_, double cell_a_, double cell_b_, double cell_c_,
             double angle_alpha_, double angle_beta_,double angle_gamma_)
{
    dimension_x=dimension_x_;
    dimension_y=dimension_y_;
    dimension_z=dimension_z_;

    cell_a=cell_a_;
    cell_b=cell_b_;
    cell_c=cell_c_;

    angle_gamma=angle_gamma_;
    angle_alpha=angle_alpha_;
    angle_beta=angle_beta_;
}

Box::Box()
{
    dimension_x=0.0;
    dimension_y=0.0;
    dimension_z=0.0;

    cell_a=0.0;
    cell_b=0.0;
    cell_c=0.0;

    angle_gamma=0.0;
    angle_alpha=0.0;
    angle_beta=0.0;
}

Box::~Box()
{
 for (long long int i=0; i<(long long int)box_atoms.size() ;i++)
 {
    delete box_atoms.at(i);
 }
 box_atoms.clear();

 for (long long int i=0; i<(long long int)box_removed_atoms.size() ;i++)
 {
    delete box_removed_atoms.at(i);
 }
 box_removed_atoms.clear();

for (long long int i=0; i<(long long int)box_removed_atoms_before_start.size() ;i++)
 {
    delete box_removed_atoms_before_start.at(i);
 }
 box_removed_atoms_before_start.clear();

 for (long long int i=0; i<(long long int)box_cells.size() ;i++)
 {
    delete box_cells.at(i);

 }
    box_cells.clear();
}

double Box::get_dimension_x()
{
    return dimension_x;
}
double Box::get_dimension_y()
{
    return dimension_y;
}
double Box::get_dimension_z()
{
    return dimension_z;
}


vector<Atom*> Box::get_box_atoms()
{
    return box_atoms;
}

void Box::set_box_atoms()
{
     long long int id=1;

    for ( long long int i= 0;i < dimension_x;i++){
        for ( long long int j= 0;j < dimension_y;j++ ){
            for ( long long int k= 0;k < dimension_z;k++ ){
                double x=j*cell_b*cos(angle_gamma*PI/180.0)
                        +k*cell_c*cos(angle_beta*PI/180.0)
                        +i*cell_a;
                double y=j*cell_b*sin(angle_gamma*PI/180.0)
                        +k*cell_c*cos(angle_alpha*PI/180.0);
                double z=k*cell_c*sin(angle_alpha*PI/180.0)*sin(angle_beta*PI/180.0);

                Atom * atombox= new Atom(id,x,y,z,1,false);  //WACTH OUT
                add_atom(atombox);


                id++;
            }
        }
    }
}

void Box::add_atom(Atom *atom)
{
    box_atoms.push_back(atom);
}

Atom * Box::add_atom( long long int id, string atom_type, long long int type ,double x, double y, double z)
{
    Atom * pointatom= new Atom(id,x,y,z,atom_type,false);
    pointatom->set_type(type);
    add_atom(pointatom);
    return pointatom;
}

long long int Box::rm_atom(long long int id)
{

    rm_surface_atom(id);  // maybe it is not in the surface if it is affected

    //Se puede paralelizar.. ya se hara si eso
    for ( long long int i=0; i<(long long int)box_atoms.size(); i++)
    {
        if (id == box_atoms.at(i)->get_id())
        {
            if (box_atoms.at(i)->get_type()==NORMAL)//si es type normal
            {
                Atom * todelete = box_atoms.at(i);
                //atoms don't have it as  neighbor anymore
                //-----------------------------
                for (long long int j=0; j<box_atoms.at(i)->get_size_neighbour(); j++)
                {
                    for (long long int k=0; k<box_atoms.at(i)->get_neighbour(j)->get_size_neighbour(); k++)
                    {

                        if (box_atoms.at(i)->get_neighbour(j)->get_id_neighbour(k) == id)
                        {
                            long long int a=box_atoms.at(i)->get_neighbour(j)->rm_neighbour(id); //a should be 0
                            cout << "i   " << i << "  j   " << j << "  k   " << k <<"  a   " << a <<endl;
                        }
                    }
                }
                //-----------------------------
                todelete->set_type(DISSOLVED);
                box_removed_atoms.push_back(todelete);
                box_atoms.erase(box_atoms.begin()+i);
                //delete todelete;  //deberia liberar la memoria guardada para el atomo Ya no lo borramos
                return 0;
            }

        }
    }

    return 1;//fail  //puede pasar que haya un atomo en dos listas de affected, si ya esta quitado.. ¿mala suerte? El atomo no sabe que esta affected(solo conoce sus vecinos y sus affected)
}


long long int Box::rm_atom_complex(long long int id)
{
     long long int value;
    volatile bool flag=false;


    rm_surface_atom(id);  // maybe it is not in the surface if it is affected

    #pragma omp parallel for shared(flag)
    for ( long long int i=0; i<(long long int)box_atoms.size(); i++)
    {
        if (flag) continue;
        if (id == box_atoms.at(i)->get_id())
        {
            if (box_atoms.at(i)->get_type()==NORMAL)//if the type is normal
            {
                Atom * todelete = box_atoms.at(i);
                vector<Atom*> neighs=todelete->get_neighbours();

                //atoms don't have it as  neighbor anymore
                //-----------------------------
                for ( long long int j=0; j<(long long int)neighs.size(); j++)
                {
                    neighs.at(j)->rm_neighbour_and_record_both_directions(id);   //we are removing records from neighbors
                }
                //-----------------------------
                todelete->set_type(DISSOLVED);
                box_removed_atoms.push_back(todelete);

                value=i;
                flag=true;
            }
            /**
            if (box_atoms.at(i)->get_type()==INDISLOCATION)//if the type is normal
            {
                Atom * todelete = box_atoms.at(i);
                vector<Atom*> neighs=todelete->get_neighbours();
                vector<Atom*> affecteds=todelete->get_affected();

                //atoms don't have it as  neighbor anymore
                //-----------------------------
                for ( long long int j=0; j<neighs.size(); j++)
                {
                    if(neighs.at(j)->get_type()==NORMAL){

                        neighs.at(j)->rm_neighbour_and_record_both_directions(id);   //we are removing records from neighbors

                    }

                }

                for (long long int k=0;k<affecteds.size();k++)
                {
                        rm_atom_complex(affecteds.at(k)->get_id());
                }
                //-----------------------------
                todelete->set_type(DISSOLVED);
                box_removed_atoms.push_back(todelete);

                value=i;
                flag=true;
            }*/

        }
    }
    if (flag) {box_atoms.erase(box_atoms.begin()+value); return 0;}
    return 1;//fail
}

void Box::set_neighbour(bool periodicity_x,bool periodicity_y,bool periodicity_z)
{
    //se puede paralelizar.. y seguramente se haga
    for ( long long int i=0; i < dimension_x*dimension_y*dimension_z ;i++)
    {
        //atomo al que queremos fijar vecinos

    // vecinos arriba

        //cout << "box_atoms.at(0)->get_z()  " << box_atoms.at(0)->get_z() << "   (dimension_z-1)*cell_c   " << (dimension_z-1)*cell_c << endl;
        if (Util::isEqual(box_atoms.at(i)->get_z(),(dimension_z-1)*cell_c))  //si esta en la parte superior, consideramos vecinos ,pero luego no eventos
        {
            //box1.get_box_atoms().at(2)->add_neighbour(box1.get_box_atoms().at(1));
            //cout << "i   " << i << "dimension_z-1   " << dimension_z-1 << endl;
            //cout << "box_atoms.size()   " << box_atoms.size() << "i-(dimension_z-1)   " << i-(dimension_z-1) << endl;
            if (periodicity_z) {box_atoms.at(i)->add_neighbour(box_atoms.at(i-(dimension_z-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i+1));
        }

        //vecinos abajo
        if (Util::isEqual(box_atoms.at(i)->get_z(),0.0) )
        {
            if (periodicity_z) {box_atoms.at(i)->add_neighbour(box_atoms.at(i+(dimension_z-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i-1));
        }



    //vecinos derecha y


        if (Util::isEqual(box_atoms.at(i)->get_y(),(dimension_y-1)*cell_b))  //si esta en la parte derecha en y, tendrá vecinos en la parte izquierda del box
        {
            if (periodicity_y) {box_atoms.at(i)->add_neighbour(box_atoms.at(i-dimension_z*(dimension_y-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i+dimension_z));
        }


    //vecinos izquierda y
        if (Util::isEqual(box_atoms.at(i)->get_y(),0.0))  //si esta en la parte izquierda en y, tendrá vecinos en la parte derecha del box
        {
            if (periodicity_y) {box_atoms.at(i)->add_neighbour(box_atoms.at(i+dimension_z*(dimension_y-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i-dimension_z));
        }

    //vecinos derecha x

        if (Util::isEqual(box_atoms.at(i)->get_x(),(dimension_x-1)*cell_a))  //si esta en la parte derecha en x, tendrá vecinos en la parte izquierda del box
        {
            if (periodicity_x) {box_atoms.at(i)->add_neighbour(box_atoms.at(i-dimension_z*dimension_y*(dimension_x-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i+dimension_y*dimension_z));
        }

    //vecinos izquierda x
        if (Util::isEqual(box_atoms.at(i)->get_x(),0.0))  //si esta en la parte izquierda en y, tendrá vecinos en la parte derecha del box
        {
            if (periodicity_x) {box_atoms.at(i)->add_neighbour(box_atoms.at(i+dimension_z*dimension_y*(dimension_x-1)));}
        }
        else
        {
            box_atoms.at(i)->add_neighbour(box_atoms.at(i-dimension_y*dimension_z));
        }
    }
}

void Box::set_neighbour(double distance_) //imcompleto quitar en el futuro
{
    long long int limit=get_box_atoms().size();
    for (long long int i=0;i<limit;i++)
    {
        for (long long int j=0;j<limit;j++)
        {
            if ( Util::isEqual(distance_,box_atoms.at(i)->get_distance(box_atoms.at(j))))
            {
                box_atoms.at(i)->add_neighbour(box_atoms.at(j));
            }

        }
    }
}

void Box::mv_atom(Atom *atom, double x, double y, double z)
{

}


long long int Box::add_surface_atom(Atom *atom)
{
    if (!atom->get_insurface())
    {
        surface_atoms.push_back(atom);
        atom->set_insurface(true);          //es parte de la superficie
        return 0;
    }

    return 1;//fail
}


long long int Box::rm_surface_atom(long long int id)
{
    long long int value;
    volatile bool flag=false;

    #pragma omp parallel for shared(flag)
    for (long long int i=0; i<get_surface_size(); i++)
    {
        if (flag) continue;
        if (id == surface_atoms.at(i)->get_id())
        {
            Atom * toremovesurface = surface_atoms.at(i);
            value=i;
            toremovesurface->set_insurface(false);
            flag=true;
        }
    }
    if (flag) {surface_atoms.erase(surface_atoms.begin()+value); return 0;}
    return 1;//fail
}


vector<Atom*> Box::get_surface_atoms()
{
    return surface_atoms;
}

void Box::set_initial_surface_atoms()
{
    //ALL ATOMS WITH LESS NEIGHBOUR THAN 6 ARE IN SURFACE

    //se hara #pragma omp parallel for
    for ( long long int i=0;i<(long long int)box_atoms.size();i++)
    {
        if (box_atoms.at(i)->get_size_neighbour()!= KSNEIGH && box_atoms.at(i)->get_type()!= INSOLUBLE)
        {
            add_surface_atom(box_atoms.at(i));
        }
    }
}

//introduce the neighbors of the dissolved atoms if they are not already introduce
void Box::update_surface(Atom * dissolved_atom)
{
     long long int size=dissolved_atom->get_size_neighbour();
    for ( long long int i=0;i<size;i++)
    {
        if (!check_in_surface(dissolved_atom->get_id_neighbour(i)))
        {
            add_surface_atom(dissolved_atom->get_neighbour(i));
        }
    }

}


bool Box::check_in_surface(long long int id)
{
    volatile bool flag=false;
    #pragma omp parallel for shared(flag)
    for (long long int i=0;i<get_surface_size();i++)
    {
        if (flag) continue;
        if (surface_atoms.at(i)->get_id() == id)
        {
            flag=true;

        }
    }
    if (flag) return true;
    return false;
}

long long int Box::get_surface_size()
{
    return surface_atoms.size();
}



void Box::add_new_atom(long long int dimension_a,long long int dimension_b, long long int dimension_c, long long int type)
{
     long long int id=box_atoms.size()+1; //WACTH OUT the id -->> DON'T combine  the 2 create system methods!
    double x=from_dimension_to_x(dimension_a,dimension_b,dimension_c);
    double y=from_dimension_to_y(dimension_b,dimension_c);
    double z=from_dimension_to_z(dimension_c);


    Atom * atombox= new Atom(id,x,y,z,1,false);  //WACTH OUT
    add_atom(atombox);
    //Add neighbourhood

}





//IRREGULARITIES
//-1 type of atom is reserved for atoms that are going to be removed before start the simulation

long long int Box::remove_sym_position(long long int dimension_a, long long int dimension_b, long long int dimension_c)
{
    if (select_atom_by_position(dimension_a, dimension_b, dimension_c)!=nullptr)
    {
        select_atom_by_position(dimension_a, dimension_b, dimension_c)->set_type(FREEPOSITION);
        return 0;
    }
    return 1;
}

long long int Box::add_sym_position(long long int dimension_a, long long int dimension_b, long long int dimension_c,long long int type)
{
    if (select_atom_by_position(dimension_a, dimension_b, dimension_c)!=nullptr)
    {
        select_atom_by_position(dimension_a, dimension_b, dimension_c)->set_type(type);
        return 0;
    }
    return 1;
}

long long int Box::add_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez,long long int type)
{
    bool all=true;
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(add_sym_position(dimension_a, dimension_b+i, dimension_c+j, type)!=0){ all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::rm_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(remove_sym_position(dimension_a, dimension_b+i, dimension_c+j)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::add_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez,long long int type)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(add_sym_position(dimension_a+i, dimension_b, dimension_c+j, type)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::rm_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(remove_sym_position(dimension_a+i, dimension_b, dimension_c+j)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::add_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey,long long int type)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(add_sym_position(dimension_a+i, dimension_b+j, dimension_c, type)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::rm_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(remove_sym_position(dimension_a+i, dimension_b+j, dimension_c)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::add_xy_heli_dislocation( long long int dimension_a, long long int dimension_b, long long int radiousdisloca)
{
    if (radiousdisloca==0) return 1;
    bool all=true;

    for ( long long int i=0;i<dimension_z;i++)
    {
        for (long long int j=0;j<radiousdisloca;j++)
        {
            for (long long int k=0;k<radiousdisloca;k++)
            {
                if(select_atom_by_position(dimension_a+j, dimension_b+k, i)==nullptr) {all=false;}
                select_atom_by_position(dimension_a+j, dimension_b+k, i)->set_type(FREEPOSITION);
            }
        }
    }
    if (!all) {return 1;}
    return 0;
}

long long int Box::add_xz_heli_dislocation( long long int dimension_a, long long int dimension_c, long long int radiousdisloca)
{
    if (radiousdisloca==0) return 1;
    bool all=true;

    for ( long long int i=0;i<dimension_y;i++)
    {
        for (long long int j=0;j<radiousdisloca;j++)
        {
            for (long long int k=0;k<radiousdisloca;k++)
            {
                if(select_atom_by_position(dimension_a+j, i, dimension_c+k)==nullptr) {all=false;}
                select_atom_by_position(dimension_a+j, i,dimension_c+k)->set_type(FREEPOSITION);
            }
        }
    }
    if (!all) {return 1;}
    return 0;
}

long long int Box::add_yz_heli_dislocation( long long int dimension_b, long long int dimension_c, long long int radiousdisloca)
{
    if (radiousdisloca==0) return 1;
    bool all=true;

    for ( long long int i=0;i<dimension_x;i++)
    {
        for (long long int j=0;j<radiousdisloca;j++)
        {
            for (long long int k=0;k<radiousdisloca;k++)
            {
                if(select_atom_by_position(i, dimension_b+j, dimension_c+k)==nullptr) {all=false;}
                select_atom_by_position(i, dimension_b+j,dimension_c+k)->set_type(FREEPOSITION);
            }
        }
    }
    if (!all) {return 1;}
    return 0;
}

//type -2 reserved for insoluble atoms
long long int Box::define_insoluble_atom(long long int dimension_a, long long int dimension_b, long long int dimension_c)
{
    if (select_atom_by_position( dimension_a,  dimension_b,  dimension_c)!=nullptr)
    {
        select_atom_by_position(dimension_a, dimension_b, dimension_c)->set_type(INSOLUBLE);
        return 0;
    }
    return 1;
}

long long int Box::define_insoluble_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(define_insoluble_atom(dimension_a+i, dimension_b+j, dimension_c)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::define_insoluble_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(define_insoluble_atom(dimension_a+i, dimension_b, dimension_c+j)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}

long long int Box::define_insoluble_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(define_insoluble_atom(dimension_a, dimension_b+i, dimension_c+j)!=0) {all=false;}
        }
    }
    if(all) {return 0;}
    else return 1;
}


long long int Box::set_mass(string type, double mass)
{
     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        string atomtype=box_atoms.at(i)->get_atom_type();

        if (atomtype.compare(type)==0)
        {
            box_atoms.at(i)->set_mass(mass);
        }
    }
    return 0;
}


long long int Box::clean_gaps_before_start()
{
    bool all=true;
    //paralelizar
    for ( long long int i=0; i<(long long int)box_atoms.size();i++)
    {
        if (box_atoms.at(i)->get_type()==FREEPOSITION)
        {
            if(rm_atom(box_atoms.at(i)->get_id())!=0) {all=false;}
            else {i--;}
        }
    }
    if (!all) {return 1;}
    return 0;
}

//https://stackoverflow.com/questions/27252778/parallelise-if-else-statement-using-openmp
long long int Box::clean_gaps_before_start_complex()
{
    //bool all=true;
    //NO paralelizar
    for (  long long int i=0; i<(long long int)box_atoms.size();i++)
    {
        if (box_atoms.at(i)->get_type()==FREEPOSITION)
        {
            Atom * todelete = box_atoms.at(i);
            long long int id=todelete->get_id();
            vector<Atom*> neighs=todelete->get_neighbours();

            //atoms don't have it as  neighbor anymore
            //-----------------------------
            for ( long long int j=0; j<(long long int)neighs.size(); j++)
            {
                neighs.at(j)->rm_neighbour_and_record_both_directions(id);   //desde los vecinos, estamos borrando los records y los linked
            }
            // quitamos tambien los affected


            //-----------------------------
            todelete->set_type(REMOVED);
            box_removed_atoms_before_start.push_back(todelete);
            box_atoms.erase(box_atoms.begin()+i);
            //delete todelete;  //deberia liberar la memoria guardada para el atomo Ya no lo borramos
            i--; // hay que restar 1 pq hemos disminuido el tamanno del array
        }
    }
    return 0;
}


long long int Box::clean_gaps_before_start_complex_from_kimerafile()
{
    //bool all=true;
    //NO paralelizar
    for (  long long int i=0; i<(long long int)box_atoms.size();i++)
    {
        if (box_atoms.at(i)->get_type()==REMOVED)
        {
            Atom * todelete = box_atoms.at(i);

            box_removed_atoms_before_start.push_back(todelete);
            box_atoms.erase(box_atoms.begin()+i);
            //delete todelete;  //deberia liberar la memoria guardada para el atomo Ya no lo borramos
            i--; // hay que restar 1 pq hemos disminuido el tamanno del array
        }
    }
    return 0;
}
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//-1 type of atom is reserved for atoms that are going to be removed before start the simulation




//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//////////////////////////////////////////IRREGULARITIES////////////////////////////////////////
//-1 type of atom is reserved for atoms that are going to be removed before start the simulation




Atom* Box::select_atom_by_position(long long int dimension_a, long long int dimension_b, long long int dimension_c)
{
    double x=from_dimension_to_x(dimension_a,dimension_b,dimension_c);
    double y=from_dimension_to_y(dimension_b,dimension_c);
    double z=from_dimension_to_z(dimension_c);

    //se podria paralelizar en un futuro
    for ( long long int i=0;i<(long long int)box_atoms.size();i++)
    {
        if (Util::isEqual(x, box_atoms.at(i)->get_x() )
            && Util::isEqual(y, box_atoms.at(i)->get_y() )
            && Util::isEqual(z, box_atoms.at(i)->get_z() ))
            {
                return box_atoms.at(i);
            }
    }
    return nullptr; //Not found
}


Atom* Box::select_atom_by_id(long long int id)
{
    volatile bool flag=false;
     long long int value;

    #pragma omp parallel for shared(flag)
    for ( long long int i=0;i<(long long int)box_atoms.size();i++)
    {
        if (flag) continue;
        if (box_atoms.at(i)->get_id()==id)
        {
            flag=true;
            value=i;
        }
    }
    if (flag) return box_atoms.at(value);
    return nullptr; //Not found
}




long long int Box::dissolve_atom(Atom * dissolvedatom)
{
    //Primero metemos los atomos en surface
    update_surface(dissolvedatom);

    //Y Luego borramos el atomo
    if(rm_atom(dissolvedatom->get_id())!=0) {return 1;}
    return 0;
}



double Box::from_dimension_to_x(long long int dimension_a,long long int dimension_b,long long int dimension_c)
{
    double x=dimension_b*cell_b*cos(angle_gamma*PI/180.0)
                        +dimension_c*cell_c*cos(angle_beta*PI/180.0)
                        +dimension_a*cell_a;
    return x;
}

double Box::from_dimension_to_y(long long int dimension_b,long long int dimension_c)
{
    double y=dimension_b*cell_b*sin(angle_gamma*PI/180.0)
                +dimension_c*cell_c*cos(angle_alpha*PI/180.0);
    return y;
}

double Box::from_dimension_to_z(long long int dimension_c)
{
    double z=dimension_c*cell_c*sin(angle_alpha*PI/180.0)*sin(angle_beta*PI/180.0);
    return z;
}



long long int Box::from_possition_to_dimension_c(double z_)
{
    long long int dimensionc=lrint(z_/cell_c/sin(angle_alpha*PI/180.0)/sin(angle_beta*PI/180.0));
    return dimensionc;
}

long long int Box::from_possition_to_dimension_b(double y_, double z_)
{
    long long int dimensionb=lrint((y_-((double)from_possition_to_dimension_c(z_))*cell_c*cos(angle_alpha*PI/180.0) )/(cell_b*sin(angle_gamma*PI/180.0)));
    return dimensionb;
}

long long int Box::from_possition_to_dimension_a(double x_, double y_, double z_)
{
    long long int dimensiona=lrint((x_-((double)from_possition_to_dimension_b(y_, z_))*cell_b*cos(angle_gamma*PI/180.0)-((double)from_possition_to_dimension_c(z_))*cell_c*cos(angle_beta*PI/180.0))/(cell_a));
    return dimensiona;
}


vector<Atom*> Box::get_box_removed_atoms()
{
    return box_removed_atoms;
}


long long int Box::create_box_by_cell(Cell * first_cell)
{

    vector<Position*> first_cell_position=first_cell->get_cell_positions(); //yo creo que esto no va a funcionar y debería retornar puntero

    long long int first_cell_positions_size=first_cell_position.size();
    long long int atom_id=0;
    long long int cell_id=0;

    mt19937 &mt = RandomGenerator::Instance().get();
    uniform_real_distribution<double> dist01(0.0, 1.0);
    double randomnumber=0.0; //=dist01(mt);

    long long int select_from_prob;
    double sum_prob=0.0;


    double posx=0.0;
    double posy=0.0;
    double posz=0.0;
    double addx=0.0;
    double addy=0.0;
    double addz=0.0;

    string selected_type;

    for ( long long int i=0; i<dimension_x ;i++)
    {
        for ( long long int j=0; j<dimension_y ;j++)
        {
            for ( long long int k=0; k<dimension_z ;k++)
            {

                //create the cell and put it into the box
                Cell * cellbox =new Cell(cell_id);
                add_cell_to_box(cellbox);

                addx=from_uvw_to_x(i,j,k);
                addy=from_uvw_to_y(j,k);
                addz=from_uvw_to_z(k);

                //for each position
                for(long long int l=0;l<first_cell_positions_size;l++)
                {
                    vector<double> probs_first_cell_position=first_cell_position.at(l)->get_probs();

                    if (probs_first_cell_position.size()!=1)
                    {
                        randomnumber=dist01(mt);
                        sum_prob=0.0;

                        for ( long long int m=0;m<(long long int)probs_first_cell_position.size();m++)
                        {
                            sum_prob=sum_prob+probs_first_cell_position.at(m);

                            if (randomnumber<sum_prob)
                            {
                                select_from_prob=m;
                                break;
                            }
                        }
                    }
                    else
                    {
                        select_from_prob=0;
                    }

                    posx=first_cell_position.at(l)->get_x();
                    posy=first_cell_position.at(l)->get_y();
                    posz=first_cell_position.at(l)->get_z();
                    selected_type=first_cell_position.at(l)->get_type(select_from_prob);

                     Atom * atom = new Atom(atom_id+l,posx+addx,posy+addy,posz+addz,selected_type,false);
                     add_atom(atom);
                     cellbox->add_atom(atom);

                }
                atom_id=atom_id+first_cell_positions_size;
                cell_id++;

            }
        }
    }
    return 0;
}

long long int Box::create_box_by_cell(Cell * first_cell,long long int seed)
{

    vector<Position*> first_cell_position=first_cell->get_cell_positions();

    long long int first_cell_positions_size=first_cell_position.size();
    long long int atom_id=0;
    long long int cell_id=0;

    mt19937 &mt = RandomGenerator::Instance().get();
    mt.seed(seed);
    uniform_real_distribution<double> dist01(0.0, 1.0);
    double randomnumber=0.0;

    long long int select_from_prob;
    double sum_prob=0.0;

    double posx=0.0;
    double posy=0.0;
    double posz=0.0;
    double addx=0.0;
    double addy=0.0;
    double addz=0.0;

    string selected_type;

    for ( long long int i=0; i<dimension_x ;i++)
    {
        for ( long long int j=0; j<dimension_y ;j++)
        {
            for ( long long int k=0; k<dimension_z ;k++)
            {

                //create cell and put it into the box

                Cell * cellbox =new Cell(cell_id);
                add_cell_to_box(cellbox);

                addx=from_uvw_to_x(i,j,k);
                addy=from_uvw_to_y(j,k);
                addz=from_uvw_to_z(k);

                //for each position
                for(long long int l=0;l<first_cell_positions_size;l++)
                {

                    vector<double> probs_first_cell_position=first_cell_position.at(l)->get_probs();

                    if (probs_first_cell_position.size()!=1)
                    {
                        randomnumber=dist01(mt);
                        sum_prob=0.0;

                        for ( long long int m=0;m<(long long int)probs_first_cell_position.size();m++)
                        {
                            sum_prob=sum_prob+probs_first_cell_position.at(m);

                            if (randomnumber<sum_prob)
                            {
                                select_from_prob=m;
                                break;
                            }
                        }
                    }
                    else
                    {
                        select_from_prob=0;
                    }

                    posx=first_cell_position.at(l)->get_x();
                    posy=first_cell_position.at(l)->get_y();
                    posz=first_cell_position.at(l)->get_z();
                    selected_type=first_cell_position.at(l)->get_type(select_from_prob);



                     Atom * atom = new Atom(atom_id+l,posx+addx,posy+addy,posz+addz,selected_type,false);
                     add_atom(atom);
                     cellbox->add_atom(atom);

                }
                atom_id=atom_id+first_cell_positions_size;
                cell_id++;

            }
        }
    }
    return 0;
}

long long int Box::set_cell_neighbours_26(bool periodicity_x,bool periodicity_y,bool periodicity_z)
{
    //26 neigbours: 9 bot, 8 around, y 9 top

     long long int box_cells_size=box_cells.size();

    #pragma omp parallel for
    for (  long long int i=0; i < box_cells_size ;i++)
    {
        Cell * cell1=box_cells.at(i);
        long long int id_cell= cell1->get_cell_id();

        long long int factorxa=1;
        long long int factorxb=1;
        long long int factorya=1;
        long long int factoryb=1;
        long long int factorza=1;
        long long int factorzb=1;

        if (id_cell%dimension_z==dimension_z-1 )
        {
            factorza=-(dimension_z-1);
            if (periodicity_z) {cell1->set_cornerz(1);}
        }

        if (id_cell%dimension_z==0 )
        {
            factorzb=-(dimension_z-1);
            if (periodicity_z) {cell1->set_cornerz(-1);}
        }
        if ( (id_cell/dimension_z)%dimension_y==dimension_y-1 )
        {
            factorya=-(dimension_y-1);
            if (periodicity_y) {cell1->set_cornery(1);}
        }

        if ( (id_cell/dimension_z)%dimension_y==0 )
        {
            factoryb=-(dimension_y-1);
            if (periodicity_y) {cell1->set_cornery(-1);}
        }

        if ((id_cell/(dimension_y*dimension_z))%dimension_x==dimension_x-1 )
        {
            factorxa=-(dimension_x-1);
            if (periodicity_x) {cell1->set_cornerx(1);}
        }

        if ((id_cell/(dimension_y*dimension_z))%dimension_x==0 )
        {
            factorxb=-(dimension_x-1);
            if (periodicity_x) {cell1->set_cornerx(-1);}
        }

        if (factorza==1)                      cell1->add_cell_neighbour(box_cells.at(i+1*factorza)); //top
        if (factorza!=1 && periodicity_z)     cell1->add_cell_neighbour(box_cells.at(i+1*factorza)); //top

        if (factorzb==1) cell1->add_cell_neighbour(box_cells.at(i-1*factorzb));  //bot
        if (factorzb!=1 && periodicity_z) cell1->add_cell_neighbour(box_cells.at(i-1*factorzb));

        //j=-1
        if ((factorzb!=1 &&  !periodicity_z)  ||  (factoryb!=1 && !periodicity_y)){}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb-dimension_z*factoryb));}

        if((factorzb!=1 &&  !periodicity_z) || (factorya!=1 && !periodicity_y))   {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb+dimension_z*factorya));}

        if ((factorxb!=1 && !periodicity_x) || (factorzb!=1 &&  !periodicity_z))     {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb-dimension_z*dimension_y*factorxb));}

        if ((factorxa!=1 && !periodicity_x)|| (factorzb!=1 && !periodicity_z))  {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb+dimension_z*dimension_y*factorxa));}

        if ((factorya!=1 &&  !periodicity_y) || (factorzb!=1 && !periodicity_z) || ( factorxa!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb+dimension_z*factorya+dimension_z*dimension_y*factorxa));}

        if ((factorya!=1  &&  !periodicity_y ) || (factorzb!=1 &&  !periodicity_z) || (factorxb!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb+dimension_z*factorya-dimension_z*dimension_y*factorxb));}

        if ((factoryb!=1 && !periodicity_y) || (factorzb!=1 &&  !periodicity_z) || (factorxb!=1 && !periodicity_x )) {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb-dimension_z*factoryb-dimension_z*dimension_y*factorxb));}

        if ((factorzb!=1 && !periodicity_z)|| (factoryb!=1 && !periodicity_y) || (factorxa!=1 && !periodicity_x))  {}
        else {cell1->add_cell_neighbour(box_cells.at(i-1*factorzb-dimension_z*factoryb+dimension_z*dimension_y*factorxa));}

        //j=0
        if (factoryb==1)                  cell1->add_cell_neighbour(box_cells.at(i-dimension_z*factoryb));
        if (factoryb!=1 && periodicity_y) cell1->add_cell_neighbour(box_cells.at(i-dimension_z*factoryb));

        if (factorya==1)                  cell1->add_cell_neighbour(box_cells.at(i+dimension_z*factorya));
        if (factorya!=1 && periodicity_y) cell1->add_cell_neighbour(box_cells.at(i+dimension_z*factorya));

        if (factorxb==1)                  cell1->add_cell_neighbour(box_cells.at(i-dimension_z*dimension_y*factorxb));
        if (factorxb!=1 && periodicity_x) cell1->add_cell_neighbour(box_cells.at(i-dimension_z*dimension_y*factorxb));

        if (factorxa==1)                  cell1->add_cell_neighbour(box_cells.at(i+dimension_z*dimension_y*factorxa));
        if (factorxa!=1 && periodicity_x) cell1->add_cell_neighbour(box_cells.at(i+dimension_z*dimension_y*factorxa));

        if ((factorya!=1 && !periodicity_y) || (factorxa!=1 && !periodicity_x))  {}
        else {cell1->add_cell_neighbour(box_cells.at(i+dimension_z*factorya+dimension_z*dimension_y*factorxa));}

        if ((factorya!=1 && !periodicity_y) || (factorxb!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+dimension_z*factorya-dimension_z*dimension_y*factorxb));}

        if ((factoryb!=1 && !periodicity_y) || (factorxb!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i-dimension_z*factoryb-dimension_z*dimension_y*factorxb));}

        if ((factoryb!=1 && !periodicity_y) || (factorxa!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i-dimension_z*factoryb+dimension_z*dimension_y*factorxa));}

        //j=1

        if ((factorza!=1 && !periodicity_z) || (factoryb!=1 && !periodicity_y)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza-dimension_z*factoryb));}

        if ((factorza!=1 && !periodicity_z) ||  (factorya!=1 && !periodicity_y))   {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza+dimension_z*factorya));}

        if ((factorza!=1 && !periodicity_z) || (factorxb!=1 && !periodicity_x) )        {}
        else { cell1->add_cell_neighbour(box_cells.at(i+1*factorza-dimension_z*dimension_y*factorxb)); }

        if ((factorza!=1 && !periodicity_z) || (factorxa!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza+dimension_z*dimension_y*factorxa));}

        if ((factorya!=1 && !periodicity_y) || (factorxa!=1 && !periodicity_x) || (factorza!=1 && !periodicity_z))   {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza+dimension_z*factorya+dimension_z*dimension_y*factorxa));}

        if ((factorza!=1 && !periodicity_z) || (factorya!=1 && !periodicity_y ) || (factorxb!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza+dimension_z*factorya-dimension_z*dimension_y*factorxb));}

        if ((factorza!=1 && !periodicity_z)|| (factoryb!=1 && !periodicity_y) || ( factorxb!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza-dimension_z*factoryb-dimension_z*dimension_y*factorxb));}

        if ((factorza!=1 && !periodicity_z) || (factoryb!=1 && !periodicity_y ) || (factorxa!=1 && !periodicity_x)) {}
        else {cell1->add_cell_neighbour(box_cells.at(i+1*factorza-dimension_z*factoryb+dimension_z*dimension_y*factorxa));}
    }

    return 0;
}

long long int Box::add_cell_to_box(Cell * cell_)
{
    box_cells.push_back(cell_);
    return box_cells.size();
}

Cell* Box::add_cell_to_box(long long int id_cell)
{
    Cell * cellbox =new Cell(id_cell);
    add_cell_to_box(cellbox);
    return cellbox;
}

double Box::from_uvw_to_x(double u, double v, double w)
{
    return Util::from_uvw_to_x(u,v,w,cell_a,cell_b,cell_c,angle_beta,angle_gamma);
}
double Box::from_uvw_to_y(double v, double w)
{
    return Util::from_uvw_to_y(v,w, cell_b, cell_c, angle_alpha, angle_beta, angle_gamma);
}
double Box::from_uvw_to_z(double w)
{
    return Util::from_uvw_to_z(w,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}

double Box::from_xyz_to_u(double x, double y, double z)
{
    return Util::from_xyz_to_u(x,y,z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}
double Box::from_xyz_to_v(double y, double z)
{
    return Util::from_xyz_to_v(y,z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}
double Box::from_xyz_to_w(double z)
{
    return Util::from_xyz_to_w(z,cell_a,cell_b,cell_c,angle_alpha,angle_beta,angle_gamma);
}

vector<Cell*> Box::get_box_cells()
{
    return box_cells;
}

long long int Box::set_atom_neighbourhood(string targettype_, vector<string> neigh_type,vector<double> distances_)
{
    long long int size_box_cell= box_cells.size();
    long long int size_neighs = neigh_type.size();



    #pragma omp parallel for
    for (long long int i = 0; i<size_box_cell;i++)
    {
        double distance_to_neigh=0.0;
        double distance_to_aux_x=0.0;
        double distance_to_aux_y=0.0;
        double distance_to_aux_z=0.0;
        double distance_to_aux_xy=0.0;
        double distance_to_aux_yz=0.0;
        double distance_to_aux_xz=0.0;
        double distance_to_aux_xyz=0.0;

        long long int size_cell_atom;
        long long int size_cell_neigh_cell;
        long long int size_cell_neigh_cell2;

        Cell * cell1;
        Cell * cell2;

        Atom * atom1;
        Atom * atom2;

        double valuex=0.0;
        double valuey=0.0;
        double valuez=0.0;

        double dim_a=0.0;
        double dim_b=0.0;
        double dim_c=0.0;

        long long int cornerx=0;
        long long int cornery=0;
        long long int cornerz=0;
        long long int cornerxy=0;
        long long int cornerxz=0;
        long long int corneryz=0;
        long long int cornerxyz=0;

        //////////////////
        //////////////////

        cell1= box_cells.at(i);

        size_cell_atom=cell1->get_cell_atoms().size();

        size_cell_neigh_cell=cell1->get_cell_neighbour().size();

        cornerx=cell1->get_cornerx();
        cornery=cell1->get_cornery();
        cornerz=cell1->get_cornerz();

        if (cornerx!=0 && cornery!=0) cornerxy=1;
        if (cornerx!=0 && cornerz!=0) cornerxz=1;
        if (cornery!=0 && cornerz!=0) corneryz=1;
        if (cornerx!=0 && cornery!=0 && cornerz!=0) cornerxyz=1;

        for (long long int j=0; j<size_cell_atom;j++)
        {
            atom1=cell1->get_cell_atoms().at(j);

            Atom atom_aux_x=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_y=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_z=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xy=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_yz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xyz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());

            dim_a = from_xyz_to_u(atom1->get_x(), atom1->get_y(), atom1->get_z());
            dim_b= from_xyz_to_v(atom1->get_y(), atom1->get_z());
            dim_c= from_xyz_to_w(atom1->get_z());

            //Auxiliar atoms positions

            if (cornerx != 0 )
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b,dim_c);
                atom_aux_x.set_x(valuex);
            }
            if (cornery != 0 )
            {
                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery),dim_c);
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery),dim_c);
                atom_aux_y.set_x(valuex);
                atom_aux_y.set_y(valuey);
            }
            if (cornerz != 0 )
            {
                valuex=from_uvw_to_x(dim_a,dim_b,dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                atom_aux_z.set_x(valuex); atom_aux_z.set_y(valuey); atom_aux_z.set_z(valuez);
            }
            if (cornerxy!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b-((double)dimension_y*cornery),dim_c);
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery),dim_c);
                 atom_aux_xy.set_x(valuex); atom_aux_xy.set_y(valuey);
            }
            if (cornerxz!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b,dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                 atom_aux_xz.set_x(valuex); atom_aux_xz.set_y(valuey); atom_aux_xz.set_z(valuez);
            }
            if (corneryz!=0)
            {
                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery),dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery), dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                 atom_aux_yz.set_x(valuex); atom_aux_yz.set_y(valuey); atom_aux_yz.set_z(valuez);
            }
            if (cornerxyz!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b-((double)dimension_y*cornery),dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery), dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                atom_aux_xyz.set_x(valuex); atom_aux_xyz.set_y(valuey); atom_aux_xyz.set_z(valuez);
            }

            if (targettype_.compare(atom1->get_atom_type())==0)
            {
                //we have the atom, now we have to see neighbors atoms and the same cell. And closer cells
                //own cell

                for (long long int k=0;k<size_cell_atom;k++)
                {
                    atom2=cell1->get_cell_atoms().at(k);
                    distance_to_neigh=atom1->get_distance(atom2);
                    //id distance and type is the same. we save it into neighbors
                    for (long long int s=0;s<size_neighs;s++)
                    {
                        if (Util::isEqualDistances(distance_to_neigh,distances_.at(s)) && atom2->get_atom_type().compare(neigh_type.at(s))==0)
                        {
                            atom1->add_neighbour(atom2);
                            atom1->add_distance_to_neighbours(distances_.at(s));
                            atom1->add_linked_type(NO_LINKED_TYPE);
                            atom1->add_void_vector_in_linked();
                            break;
                        }
                    }
                }

                //rest of cells
                //we have to take into account if a cell is in a corner
                //For that we can create auxiliar atoms with displaced possition, and if there are close, put them into the original ones
                //7 ayxiliar atoms

                for (long long int r=0; r<size_cell_neigh_cell; r++)
                {
                    cell2=cell1->get_cell_neighbour().at(r);
                    size_cell_neigh_cell2=cell2->get_cell_atoms().size();

                    //atoms inside the neighbor cell
                    for (long long int h=0;h<size_cell_neigh_cell2;h++)
                    {
                        atom2=cell2->get_cell_atoms().at(h);
                        distance_to_neigh=atom1->get_distance(atom2);
                        distance_to_aux_x=atom_aux_x.get_distance(atom2);
                        distance_to_aux_y=atom_aux_y.get_distance(atom2);
                        distance_to_aux_z=atom_aux_z.get_distance(atom2);
                        distance_to_aux_xy=atom_aux_xy.get_distance(atom2);
                        distance_to_aux_xz=atom_aux_xz.get_distance(atom2);
                        distance_to_aux_yz=atom_aux_yz.get_distance(atom2);
                        distance_to_aux_xyz=atom_aux_xyz.get_distance(atom2);

                        //que pasa si hay varios a la vez que cumplen la condicion? solo estaremos metiendo 1..de eso se trata..
                        for (long long int s=0;s<size_neighs;s++)
                        {
                            if ((Util::isEqualDistances(distance_to_neigh,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_x,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_y,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_z,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xy,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xz,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_yz,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xyz,distances_.at(s))
                                 ) && atom2->get_atom_type().compare(neigh_type.at(s))==0)
                            {
                                atom1->add_neighbour(atom2);
                                atom1->add_distance_to_neighbours(distances_.at(s));
                                atom1->add_linked_type(NO_LINKED_TYPE);
                                atom1->add_void_vector_in_linked();
                                break;
                            }
                        }

                    }

                }

            }
                //atom1->set_neigh_record(); //build the neighbor record ESTO DA MAL SI LO LLAMAMOS DOS O MAS VECES.. ESTAMOS HACIENDO EL RECORD MAS DE 1 VEZ
        }
        cornerx=0;
        cornery=0;
        cornerz=0;
        cornerxy=0;
        cornerxz=0;
        corneryz=0;
        cornerxyz=0;
    }
    return 0;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

long long int Box::set_atom_neighbourhood_linked(string targettype_, vector<string> neigh_type,vector<double> distances_, vector<Linked_neighbour> linked_neighbours)
{
    long long int size_box_cell= box_cells.size();
    long long int size_neighs = neigh_type.size();

   //long long int  size_linked_neighbours= linked_neighbours.size(); deberia ser el mismo tamaño que el anterior


    #pragma omp parallel for  //igual deberia quitar esto si quiero que los linked esten sincronizados con los vecinos, NO HACE FALTA, VA A ESTAR EN ORDEN
    for (long long int i = 0; i<size_box_cell;i++)
    {
        double distance_to_neigh=0.0;
        double distance_origin_to_link=0.0;
        double distance_target_to_link=0.0;

        long long int size_cell_atom;
        long long int size_cell_neigh_cell;
        long long int size_cell_neigh_cell2;

        Cell * cell1;
        Cell * cell2;

        Atom * atom1; // el propio
        Atom * atom2; // el vecino
        Atom * atom3; // el linked al vecino


        double valuex=0.0;
        double valuey=0.0;
        double valuez=0.0;

        double dim_a=0.0;
        double dim_b=0.0;
        double dim_c=0.0;

        //////////////////
        //////////////////

        cell1= box_cells.at(i);

        size_cell_atom=cell1->get_cell_atoms().size();

        size_cell_neigh_cell=cell1->get_cell_neighbour().size();


        for (long long int j=0; j<size_cell_atom;j++)
        {
            atom1=cell1->get_cell_atoms().at(j);

            if (targettype_.compare(atom1->get_atom_type())==0)
            {
                //we have the atom, now we have to see neighbors atoms and the same cell. And closer cells
                //own cell
                //cout<<"HOLA10"<<endl;

                for (long long int s=0;s<size_neighs;s++)
                {

                    for (long long int k=0;k<size_cell_atom;k++)
                    {
                        atom2=cell1->get_cell_atoms().at(k);

                        distance_to_neigh=atom1->get_distance(atom2);
                        //id distance and type is the same. we save it into neighbors
                        //cout<<"HOLA20"<<endl;

                        //cout<<"HOLA25"<<endl;

                            if (Util::isEqualDistances(distance_to_neigh,distances_.at(s)) && atom2->get_atom_type().compare(neigh_type.at(s))==0)
                            {
                                atom1->add_neighbour(atom2);
                                atom1->add_distance_to_neighbours(distances_.at(s));
                                atom1->add_linked_type(linked_neighbours.at(s).get_type());
                                atom1->add_void_vector_in_linked();


                                //habra que annadir aqui los linked asociados al vecinos.. y no solo en la mima celda, sino mirando las contiguas.. :/
                                if (linked_neighbours.at(s).get_type() !=NO_LINKED_TYPE)
                                {
                                    //cout<<"HOLA50"<<endl;
                                    vector <string> type_linked=linked_neighbours.at(s).get_type_linked(); //tipos de los vacinos que son linked
                                    vector <double> distances_origin_to_link=linked_neighbours.at(s).get_distance_origin_linked();
                                    vector <double> distances_target_to_link=linked_neighbours.at(s).get_distance_target_linked();

                                    //CASO EN EL QUE TANTO EL VECINO COMO EL LINK ESTAN EN LA MISMA CELDA
                                    for (long long int v=0;v<size_cell_atom;v++)
                                    {
                                        atom3=cell1->get_cell_atoms().at(v);

                                        distance_origin_to_link=atom1->get_distance(atom3);
                                        distance_target_to_link=atom2->get_distance(atom3);


                                        for (long long int m=0;m<(long long int)type_linked.size();m++)
                                        {
                                            if (atom3->get_atom_type().compare(type_linked.at(m))==0 && Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                            //le añadimos a lista de linked, y le ponemos el id del linked_neighbours que corresponda
                                            {
                                                atom1->add_linked_to(atom3,atom1->get_size_neighbour()-1);  //CREO QUE NO  DEBERIA METERSE EN a s, casi seguro que no CHECKEAR
                                                break;
                                            }
                                        }
                                    }

                                    //cout<<"HOLA100"<<endl;
                                    //CASO EN EL QUE VECINO ESTA EN LA MISMA CELDA, PERO EL LINKED EN LA CONTIGUA
                                    for (long long int r=0;r<size_cell_neigh_cell;r++)
                                    {
                                        cell2=cell1->get_cell_neighbour().at(r);
                                        size_cell_neigh_cell2=cell2->get_cell_atoms().size();

                                        //atoms inside the neighbor cell
                                        for (long long int h=0;h<size_cell_neigh_cell2;h++)
                                        {
                                            atom3=cell2->get_cell_atoms().at(h);

                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////


                                            long long int cornerx3=cell2->get_cornerx();
                                            long long int cornery3=cell2->get_cornery();
                                            long long int cornerz3=cell2->get_cornerz();
                                            long long int cornerxy3=0;
                                            long long int cornerxz3=0;
                                            long long int corneryz3=0;
                                            long long int cornerxyz3=0;

                                            if (cornerx3!=0 && cornery3!=0) cornerxy3=1;
                                            if (cornerx3!=0 && cornerz3!=0) cornerxz3=1;
                                            if (cornery3!=0 && cornerz3!=0) corneryz3=1;
                                            if (cornerx3!=0 && cornery3!=0 && cornerz3!=0) cornerxyz3=1;


                                            Atom atom3_aux_x=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_y=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_z=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xy=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_yz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xyz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());

                                            dim_a = from_xyz_to_u(atom3->get_x(), atom3->get_y(), atom3->get_z());
                                            dim_b= from_xyz_to_v(atom3->get_y(), atom3->get_z());
                                            dim_c= from_xyz_to_w(atom3->get_z());

                                            if (cornerx3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b,dim_c);
                                                atom3_aux_x.set_x(valuex);
                                            }
                                            if (cornery3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery3),dim_c);
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3),dim_c);
                                                atom3_aux_y.set_x(valuex);
                                                atom3_aux_y.set_y(valuey);
                                            }
                                            if (cornerz3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                atom3_aux_z.set_x(valuex); atom3_aux_z.set_y(valuey); atom3_aux_z.set_z(valuez);
                                            }
                                            if (cornerxy3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b-((double)dimension_y*cornery3),dim_c);
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3),dim_c);
                                                 atom3_aux_xy.set_x(valuex); atom3_aux_xy.set_y(valuey);
                                            }
                                            if (cornerxz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                 atom3_aux_xz.set_x(valuex); atom3_aux_xz.set_y(valuey); atom3_aux_xz.set_z(valuez);
                                            }
                                            if (corneryz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery3),dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3), dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                 atom3_aux_yz.set_x(valuex); atom3_aux_yz.set_y(valuey); atom3_aux_yz.set_z(valuez);
                                            }
                                            if (cornerxyz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b-((double)dimension_y*cornery3),dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3), dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                atom3_aux_xyz.set_x(valuex); atom3_aux_xyz.set_y(valuey); atom3_aux_xyz.set_z(valuez);
                                            }

                                            distance_origin_to_link=atom3->get_distance(atom1);
                                            distance_target_to_link=atom3->get_distance(atom2);

                                            double distance_origin_to_linkauxiliar_x=atom3_aux_x.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_x=atom3_aux_x.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_y=atom3_aux_y.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_y=atom3_aux_y.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_z=atom3_aux_z.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_z=atom3_aux_z.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_xy=atom3_aux_xy.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_xy=atom3_aux_xy.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_xz=atom3_aux_xz.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_xz=atom3_aux_xz.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_yz=atom3_aux_yz.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_yz=atom3_aux_yz.get_distance(atom2);
                                            double distance_origin_to_linkauxiliar_xyz=atom3_aux_xyz.get_distance(atom1);
                                            double distance_target_to_linkauxiliar_xyz=atom3_aux_xyz.get_distance(atom2);



                                            for (long long int m=0;m<(long long int)type_linked.size();m++)
                                            {
                                                if (((Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_x, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_x, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_y, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_y, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_z, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_z, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xy, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xy, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xz, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_yz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_yz, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xyz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xyz, distances_target_to_link.at(m)))
                                                     )   && atom3->get_atom_type().compare(type_linked.at(m))==0)
                                                 {
                                                        atom1->add_linked_to(atom3,atom1->get_size_neighbour()-1);
                                                        break;
                                                 }
                                            }
                                        }

                                    }
                                    //cout<<"HOLA200"<<endl;
                                    atom1->remove_neighbours_with_out_linkeds();

                                }

                                //cout<<"HOLA250"<<endl;

                            }


                        //cout<<"HOLA300"<<endl;
                    }
                    //cout<<"HOLA350"<<endl;

                    //rest of cells
                    //we have to take into account if a cell is in a corner
                    //For that we can create auxiliar atoms with displaced possition, and if there are close, put them into the original ones
                    //7 ayxiliar atoms

                    for (long long int r=0; r<size_cell_neigh_cell; r++)
                    {
                        cell2=cell1->get_cell_neighbour().at(r);
                        size_cell_neigh_cell2=cell2->get_cell_atoms().size();

                        //atoms inside the neighbor cell
                        for (long long int h=0;h<size_cell_neigh_cell2;h++)
                        {
                            atom2=cell2->get_cell_atoms().at(h);

                            //Aqui nos falta el if del tipo de evento

                            ////////////////////////AUXILIARES DE ATOM2///////////////
                            ////////////////////////AUXILIARES DE ATOM2///////////////
                            ////////////////////////AUXILIARES DE ATOM2///////////////
                            ////////////////////////AUXILIARES DE ATOM2///////////////
                            long long int cornerx2=cell2->get_cornerx();
                            long long int cornery2=cell2->get_cornery();
                            long long int cornerz2=cell2->get_cornerz();
                            long long int cornerxy2=0;
                            long long int cornerxz2=0;
                            long long int corneryz2=0;
                            long long int cornerxyz2=0;

                            if (cornerx2!=0 && cornery2!=0) cornerxy2=1;
                            if (cornerx2!=0 && cornerz2!=0) cornerxz2=1;
                            if (cornery2!=0 && cornerz2!=0) corneryz2=1;
                            if (cornerx2!=0 && cornery2!=0 && cornerz2!=0) cornerxyz2=1;


                            Atom atom2_aux_x=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_y=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_z=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_xy=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_xz=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_yz=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());
                            Atom atom2_aux_xyz=Atom(atom2->get_id(),atom2->get_x(),atom2->get_y(),atom2->get_z(),atom2->get_atom_type(),atom2->get_insurface());

                            dim_a = from_xyz_to_u(atom2->get_x(), atom2->get_y(), atom2->get_z());
                            dim_b= from_xyz_to_v(atom2->get_y(), atom2->get_z());
                            dim_c= from_xyz_to_w(atom2->get_z());

                            if (cornerx2 != 0 )
                            {
                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx2),dim_b,dim_c);
                                atom2_aux_x.set_x(valuex);
                            }
                            if (cornery2 != 0 )
                            {
                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery2),dim_c);
                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery2),dim_c);
                                atom2_aux_y.set_x(valuex);
                                atom2_aux_y.set_y(valuey);
                            }
                            if (cornerz2 != 0 )
                            {
                                valuex=from_uvw_to_x(dim_a,dim_b,dim_c-((double)dimension_z*cornerz2));
                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz2));
                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz2));
                                atom2_aux_z.set_x(valuex); atom2_aux_z.set_y(valuey); atom2_aux_z.set_z(valuez);
                            }
                            if (cornerxy2!=0)
                            {
                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx2),dim_b-((double)dimension_y*cornery2),dim_c);
                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery2),dim_c);
                                 atom2_aux_xy.set_x(valuex); atom2_aux_xy.set_y(valuey);
                            }
                            if (cornerxz2!=0)
                            {
                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx2),dim_b,dim_c-((double)dimension_z*cornerz2));
                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz2));
                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz2));
                                 atom2_aux_xz.set_x(valuex); atom2_aux_xz.set_y(valuey); atom2_aux_xz.set_z(valuez);
                            }
                            if (corneryz2!=0)
                            {
                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery2),dim_c-((double)dimension_z*cornerz2));
                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery2), dim_c-((double)dimension_z*cornerz2));
                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz2));
                                 atom2_aux_yz.set_x(valuex); atom2_aux_yz.set_y(valuey); atom2_aux_yz.set_z(valuez);
                            }
                            if (cornerxyz2!=0)
                            {
                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx2),dim_b-((double)dimension_y*cornery2),dim_c-((double)dimension_z*cornerz2));
                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery2), dim_c-((double)dimension_z*cornerz2));
                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz2));
                                atom2_aux_xyz.set_x(valuex); atom2_aux_xyz.set_y(valuey); atom2_aux_xyz.set_z(valuez);
                            }

                            double distance_to_neigh=atom2->get_distance(atom1);
                            double distance_to_aux_x=atom2_aux_x.get_distance(atom1);
                            double distance_to_aux_y=atom2_aux_y.get_distance(atom1);
                            double distance_to_aux_z=atom2_aux_z.get_distance(atom1);
                            double distance_to_aux_xy=atom2_aux_xy.get_distance(atom1);
                            double distance_to_aux_xz=atom2_aux_xz.get_distance(atom1);
                            double distance_to_aux_yz=atom2_aux_yz.get_distance(atom1);
                            double distance_to_aux_xyz=atom2_aux_xyz.get_distance(atom1);



                            if ((Util::isEqualDistances(distance_to_neigh,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_x,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_y,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_z,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xy,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xz,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_yz,distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xyz,distances_.at(s))
                                 ) && atom2->get_atom_type().compare(neigh_type.at(s))==0)
                            {
                                //tengo que quitar los vecinos que no esten linkados



                                atom1->add_neighbour(atom2);
                                atom1->add_distance_to_neighbours(distances_.at(s));
                                atom1->add_linked_type(linked_neighbours.at(s).get_type());
                                atom1->add_void_vector_in_linked();


                                if (linked_neighbours.at(s).get_type() !=NO_LINKED_TYPE)
                                {
                                    vector <string> type_linked=linked_neighbours.at(s).get_type_linked(); //tipos de los vacinos que son linked
                                    vector <double> distances_origin_to_link=linked_neighbours.at(s).get_distance_origin_linked();
                                    vector <double> distances_target_to_link=linked_neighbours.at(s).get_distance_target_linked();

                                    //CASO EN EL QUE EL VECINO ESTA EN LA CELDA CONTIGUA Y EL LINK EN LA MISMA
                                    for (long long int v=0;v<size_cell_atom;v++)
                                    {
                                        atom3=cell1->get_cell_atoms().at(v);




                                        distance_origin_to_link=atom1->get_distance(atom3);
                                        distance_target_to_link=atom2->get_distance(atom3);


                                        double distance_target_to_linkauxiliar_x=atom2_aux_x.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_y=atom2_aux_y.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_z=atom2_aux_z.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_xy=atom2_aux_xy.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_xz=atom2_aux_xz.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_yz=atom2_aux_yz.get_distance(atom3);
                                        double distance_target_to_linkauxiliar_xyz=atom2_aux_xyz.get_distance(atom3);

                                        //LA DISTANCIA TARGET_LINK VA A NECESITAR AUXILIARES, LA OTRA NO

                                        //CREO QUE NO ESTA PILLANDO ESTO
                                        for (long long int m=0;m<(long long int)type_linked.size();m++)
                                        {
                                            if (   ( (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_x, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_y, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_z, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xy, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xz, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_yz, distances_target_to_link.at(m)))
                                                    || (Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_linkauxiliar_xyz, distances_target_to_link.at(m))))
                                                    &&  atom3->get_atom_type().compare(type_linked.at(m))==0)
                                            //le añadimos a lista de linked, y le ponemos el id del linked_neighbours que corresponda
                                            {
                                                    atom1->add_linked_to(atom3,atom1->get_size_neighbour()-1);
                                                    break;
                                            }
                                        }
                                    }

                                    //CASO EN EL QUE AMBOS ESTAN EN LA CELDA CONTIGUA ¿puede estar cada uno en una celda contigua distinta?
                                    for (long long int r=0;r<size_cell_neigh_cell;r++)
                                    {
                                       Cell * cell4=cell1->get_cell_neighbour().at(r);   //creo que aqui estamos modificando cell4
                                       long long int size_cell_neigh_cell4=cell4->get_cell_atoms().size();

                                        //atoms inside the neighbor cell
                                        for (long long int h=0;h<size_cell_neigh_cell4;h++)
                                        {
                                            atom3=cell4->get_cell_atoms().at(h);

                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////
                                            ////////////////////////AUXILIARES DE ATOM3///////////////

                                            long long int cornerx3=cell4->get_cornerx();
                                            long long int cornery3=cell4->get_cornery();
                                            long long int cornerz3=cell4->get_cornerz();
                                            long long int cornerxy3=0;
                                            long long int cornerxz3=0;
                                            long long int corneryz3=0;
                                            long long int cornerxyz3=0;

                                            if (cornerx3!=0 && cornery3!=0) cornerxy3=1;
                                            if (cornerx3!=0 && cornerz3!=0) cornerxz3=1;
                                            if (cornery3!=0 && cornerz3!=0) corneryz3=1;
                                            if (cornerx3!=0 && cornery3!=0 && cornerz3!=0) cornerxyz3=1;


                                            Atom atom3_aux_x=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_y=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_z=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xy=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_yz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());
                                            Atom atom3_aux_xyz=Atom(atom3->get_id(),atom3->get_x(),atom3->get_y(),atom3->get_z(),atom3->get_atom_type(),atom3->get_insurface());

                                            dim_a = from_xyz_to_u(atom3->get_x(), atom3->get_y(), atom3->get_z());
                                            dim_b= from_xyz_to_v(atom3->get_y(), atom3->get_z());
                                            dim_c= from_xyz_to_w(atom3->get_z());

                                            if (cornerx3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b,dim_c);
                                                atom3_aux_x.set_x(valuex);
                                            }
                                            if (cornery3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery3),dim_c);
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3),dim_c);
                                                atom3_aux_y.set_x(valuex);
                                                atom3_aux_y.set_y(valuey);
                                            }
                                            if (cornerz3 != 0 )
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                atom3_aux_z.set_x(valuex); atom3_aux_z.set_y(valuey); atom3_aux_z.set_z(valuez);
                                            }
                                            if (cornerxy3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b-((double)dimension_y*cornery3),dim_c);
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3),dim_c);
                                                 atom3_aux_xy.set_x(valuex); atom3_aux_xy.set_y(valuey);
                                            }
                                            if (cornerxz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                 atom3_aux_xz.set_x(valuex); atom3_aux_xz.set_y(valuey); atom3_aux_xz.set_z(valuez);
                                            }
                                            if (corneryz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery3),dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3), dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                 atom3_aux_yz.set_x(valuex); atom3_aux_yz.set_y(valuey); atom3_aux_yz.set_z(valuez);
                                            }
                                            if (cornerxyz3!=0)
                                            {
                                                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx3),dim_b-((double)dimension_y*cornery3),dim_c-((double)dimension_z*cornerz3));
                                                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery3), dim_c-((double)dimension_z*cornerz3));
                                                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz3));
                                                atom3_aux_xyz.set_x(valuex); atom3_aux_xyz.set_y(valuey); atom3_aux_xyz.set_z(valuez);
                                            }
                                            //NO HACE FALTA AUXILIAR TARGET_TO_LINK
                                            //HACEN FALTA AUXILIAR ORIGEN TO LINK

                                            //CREO QUE TENGO QUE METER LOS AUXILIARES DE ATOM 2 PQ PUEDE SER QUE CADA VECINO ESTE EN CELDAS DISTINTAS
                                            //Y JUSTO COINCIDA CON QUE EL LINK ESTA EN LA PARTE PERIODICA

                                           // if(atom1->get_id()==1621  && atom2->get_id()==5  && atom3->get_id()==2)
                                            //{


                                            distance_origin_to_link=atom1->get_distance(atom3);
                                            distance_target_to_link=atom2->get_distance(atom3);

                                            double distance_origin_to_linkauxiliar_x=atom3_aux_x.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_y=atom3_aux_y.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_z=atom3_aux_z.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_xy=atom3_aux_xy.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_xz=atom3_aux_xz.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_yz=atom3_aux_yz.get_distance(atom1);
                                            double distance_origin_to_linkauxiliar_xyz=atom3_aux_xyz.get_distance(atom1);

                                            double distance_target_to_linkauxiliar_x=atom2_aux_x.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_y=atom2_aux_y.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_z=atom2_aux_z.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_xy=atom2_aux_xy.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_xz=atom2_aux_xz.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_yz=atom2_aux_yz.get_distance(atom3);
                                            double distance_target_to_linkauxiliar_xyz=atom2_aux_xyz.get_distance(atom3);




                                            for (long long int m=0;m<(long long int)type_linked.size();m++)
                                            {
                                                if (atom3->get_atom_type().compare(type_linked.at(m))==0)
                                                {
                                                    if (   Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_x, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_y, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_z, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_xy, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_xz, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_yz, distances_origin_to_link.at(m))
                                                        || Util::isEqualDistances(distance_origin_to_linkauxiliar_xyz, distances_origin_to_link.at(m))
                                                        )
                                                    {
                                                        if (   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_x, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_y, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_z, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_xy, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_xz, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_yz, distances_target_to_link.at(m))
                                                            || Util::isEqualDistances(distance_target_to_linkauxiliar_xyz, distances_target_to_link.at(m))
                                                            )
                                                        {
                                                            atom1->add_linked_to(atom3,atom1->get_size_neighbour()-1);
                                                            break;
                                                        }

                                                    }
                                                }
                                                /**
                                                if (((Util::isEqualDistances(distance_origin_to_link, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_x, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_y, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_z, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xy, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_yz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     || (Util::isEqualDistances(distance_origin_to_linkauxiliar_xyz, distances_origin_to_link.at(m)) &&   Util::isEqualDistances(distance_target_to_link, distances_target_to_link.at(m)))
                                                     )   && atom3->get_atom_type().compare(type_linked.at(m))==0)
                                                 {







                                                        //le añadimos a lista de linked, y le ponemos el id del linked_neighbours que corresponda


                                                 }
                                                 */
                                            }
                                            //}
                                        }
                                    }
                                    atom1->remove_neighbours_with_out_linkeds();

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}








long long int Box::set_atom_affected(string targettype_, vector<string> affected_type,vector<double> affected_distances_) //similar to the previous
{

    long long int size_box_cell= box_cells.size();
    long long int size_affecteds = affected_type.size();

    #pragma omp parallel for
    for (long long int i = 0; i<size_box_cell;i++)
    {

        double distance_to_affected=0.0;
        double distance_to_aux_x=0.0;
        double distance_to_aux_y=0.0;
        double distance_to_aux_z=0.0;
        double distance_to_aux_xy=0.0;
        double distance_to_aux_yz=0.0;
        double distance_to_aux_xz=0.0;
        double distance_to_aux_xyz=0.0;

        long long int size_cell_atom;
        long long int size_cell_neigh_cell;
        long long int size_cell_neigh_cell2;
        Cell * cell1;
        Cell * cell2;
        Atom * atom1;
        Atom * atom2;

        double valuex=0.0;
        double valuey=0.0;
        double valuez=0.0;

        double dim_a=0.0;
        double dim_b=0.0;
        double dim_c=0.0;

        long long int cornerx=0;
        long long int cornery=0;
        long long int cornerz=0;
        long long int cornerxy=0;
        long long int cornerxz=0;
        long long int corneryz=0;
        long long int cornerxyz=0;

        //////////////////////////////
        //////////////////////////////

        cell1= box_cells.at(i);

        size_cell_atom=cell1->get_cell_atoms().size();

        size_cell_neigh_cell=cell1->get_cell_neighbour().size();

        cornerx=cell1->get_cornerx();
        cornery=cell1->get_cornery();
        cornerz=cell1->get_cornerz();

        if (cornerx!=0 && cornery!=0) cornerxy=1;
        if (cornerx!=0 && cornerz!=0) cornerxz=1;
        if (cornery!=0 && cornerz!=0) corneryz=1;
        if (cornerx!=0 && cornery!=0 && cornerz!=0) cornerxyz=1;


        for (long long int j=0; j<size_cell_atom;j++)
        {
            atom1=cell1->get_cell_atoms().at(j);

            Atom atom_aux_x=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_y=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_z=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xy=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_yz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());
            Atom atom_aux_xyz=Atom(atom1->get_id(),atom1->get_x(),atom1->get_y(),atom1->get_z(),atom1->get_atom_type(),atom1->get_insurface());

            dim_a = from_xyz_to_u(atom1->get_x(), atom1->get_y(), atom1->get_z());
            dim_b= from_xyz_to_v(atom1->get_y(), atom1->get_z());
            dim_c= from_xyz_to_w(atom1->get_z());

            if (cornerx != 0 )
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b,dim_c);
                atom_aux_x.set_x(valuex);
            }
            if (cornery != 0 )
            {
                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery),dim_c);
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery),dim_c);
                atom_aux_y.set_x(valuex);
                atom_aux_y.set_y(valuey);
            }
            if (cornerz != 0 )
            {
                valuex=from_uvw_to_x(dim_a,dim_b,dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                atom_aux_z.set_x(valuex); atom_aux_z.set_y(valuey); atom_aux_z.set_z(valuez);
            }
            if (cornerxy!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b-((double)dimension_y*cornery),dim_c);
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery),dim_c);
                 atom_aux_xy.set_x(valuex); atom_aux_xy.set_y(valuey);
            }
            if (cornerxz!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b,dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b,dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                 atom_aux_xz.set_x(valuex); atom_aux_xz.set_y(valuey); atom_aux_xz.set_z(valuez);
            }
            if (corneryz!=0)
            {
                valuex=from_uvw_to_x(dim_a,dim_b-((double)dimension_y*cornery),dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery), dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                 atom_aux_yz.set_x(valuex); atom_aux_yz.set_y(valuey); atom_aux_yz.set_z(valuez);
            }
            if (cornerxyz!=0)
            {
                valuex=from_uvw_to_x(dim_a-((double)dimension_x*cornerx),dim_b-((double)dimension_y*cornery),dim_c-((double)dimension_z*cornerz));
                valuey=from_uvw_to_y(dim_b-((double)dimension_y*cornery), dim_c-((double)dimension_z*cornerz));
                valuez=from_uvw_to_z(dim_c-((double)dimension_z*cornerz));
                atom_aux_xyz.set_x(valuex); atom_aux_xyz.set_y(valuey); atom_aux_xyz.set_z(valuez);
            }

            if (targettype_.compare(atom1->get_atom_type())==0)
            {
                for (long long int k=0;k<size_cell_atom;k++)
                {
                    atom2=cell1->get_cell_atoms().at(k);
                    distance_to_affected=atom1->get_distance(atom2);

                    for (long long int s=0;s<size_affecteds;s++)
                    {
                        if (Util::isEqualDistances(distance_to_affected,affected_distances_.at(s)) && atom2->get_atom_type().compare(affected_type.at(s))==0)
                        {
                            atom1->add_affected(atom2);
                            break;
                        }
                    }
                }

                for (long long int r=0; r<size_cell_neigh_cell; r++)
                {
                    cell2=cell1->get_cell_neighbour().at(r);
                    size_cell_neigh_cell2=cell2->get_cell_atoms().size();

                    for (long long int h=0;h<size_cell_neigh_cell2;h++)
                    {
                        atom2=cell2->get_cell_atoms().at(h);
                        distance_to_affected=atom1->get_distance(atom2);
                        distance_to_aux_x=atom_aux_x.get_distance(atom2);
                        distance_to_aux_y=atom_aux_y.get_distance(atom2);
                        distance_to_aux_z=atom_aux_z.get_distance(atom2);
                        distance_to_aux_xy=atom_aux_xy.get_distance(atom2);
                        distance_to_aux_xz=atom_aux_xz.get_distance(atom2);
                        distance_to_aux_yz=atom_aux_yz.get_distance(atom2);
                        distance_to_aux_xyz=atom_aux_xyz.get_distance(atom2);

                        for (long long int s=0;s<size_affecteds;s++)
                        {
                            if ((Util::isEqualDistances(distance_to_affected,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_x,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_y,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_z,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xy,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xz,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_yz,affected_distances_.at(s)) ||
                                 Util::isEqualDistances(distance_to_aux_xyz,affected_distances_.at(s))
                                 ) && atom2->get_atom_type().compare(affected_type.at(s))==0)
                            {
                                atom1->add_affected(atom2);
                                break;
                            }
                        }

                    }

                }

            }
        }
        cornerx=0;
        cornery=0;
        cornerz=0;
        cornerxy=0;
        cornerxz=0;
        corneryz=0;
        cornerxyz=0;
    }
    return 0;
}


long long int Box::set_initial_surface_atoms(vector<Event_definition*>  defined_events)
{
    //all the atoms that have less than its maximun coordination, or the are insoluble, are on the surface.

    long long int boxsize= box_atoms.size();
    long long int sizedefinedevents=defined_events.size();



    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        vector <long long int> max_to_bulk;
        vector <double> distance_neighbours;
        vector <string> type_neighbours;
        vector <Record> neigh_record;
        //vector <Linked_neighbour>  linked_neighbours;
        vector<vector<Atom*>> linkeds;
        vector<long long int> type_linkeds;

        bool next=false;

        vector<Atom*>  neighbours_involvedatom= box_atoms.at(i)->get_neighbours();
        vector<double>  distances_involvedatom= box_atoms.at(i)->get_distances_to_neighbours();



        ////////////////////
        ///////////////////

        //para cada evento definido (en el que estan los maximos vecinos)
        for (long long int j = 0;j<sizedefinedevents; j++)
        {
            //si el evento definido tiene el mismo type que el atomo que estamos mirando
            if (defined_events.at(j)->get_involved_atom_type().compare(box_atoms.at(i)->get_atom_type())==0)
            {
                neigh_record=box_atoms.at(i)->get_neigh_record();
                //tamano del record de vecinos que tiene el atomo que estamos mirando
                //size_neigh_record=neigh_record.size(); no hace falta


                //vector con los valores de vecinos maximos
                max_to_bulk= defined_events.at(j)->get_max_to_bulk();
                //vector con las distancias a los vecinos
                distance_neighbours= defined_events.at(j)->get_distance_neighbours();
                //vector con los tipos de los vecinos
                type_neighbours=defined_events.at(j)->get_type_neighbours();



                //miramos dentro del record, el tamanno del record, de los types, y de las distancias son las mismas y siguen el orden

                for ( long long int k=0; k<(long long int)max_to_bulk.size();k++)
                {
                    if (max_to_bulk.at(k)!=0) //solo le estamos metiendo a vecinos si hemos dicho que tiene maximo para ser bulk
                    {
                        long long int typelink_from_event_def=defined_events.at(j)->get_list_linked_neighbour().at(k).get_type();

                        long long int truenumber=box_atoms.at(i)->get_number_complex(type_neighbours.at(k),distance_neighbours.at(k),typelink_from_event_def);


                        //cout<<" antes  " <<endl;
                        //cout<<" type_neighbours.at(k)  " <<type_neighbours.at(k)<<endl;
                        //cout<<" distance_neighbours.at(k)  " <<distance_neighbours.at(k)<<endl;

                        //cout<<"truenumber: "<<truenumber<<endl;


                        //cout<<" despues " <<endl;
                        //cout<<" type_neighbours.at(k)  " <<type_neighbours.at(k)<<endl;
                        //cout<<" distance_neighbours.at(k)  " <<distance_neighbours.at(k)<<endl;
                        //cout<<"truenumber: "<<truenumber<<endl;

                        if (max_to_bulk.at(k)>truenumber &&
                            box_atoms.at(i)->get_type()== NORMAL)
                        {

                            #pragma omp critical
                            {
                                add_surface_atom(box_atoms.at(i));
                                box_atoms.at(i)->set_insurface(true);
                            }
                            next=true;
                            break;
                        }
                    }
                }
                if (next) break;
            }
        }
    }
    return 0;
}

Cell* Box::select_cell_by_position(long long int a, long long int b, long long int c)
{
    long long int pos=c+dimension_z*b+dimension_z*dimension_y*a;
    return box_cells.at(pos);
}


//GOOD IRREGULARITIES
//------------------------------------------------------------------------
//------------------------------------------------------------------------
long long int Box::change_type_sym_cell(long long int dim_a, long long int dim_b, long long int dim_c,long long int type)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
        long long int sizeatoms=cellatoms.size();

        for (long long int i=0;i<sizeatoms;i++)
        {
            cellatoms.at(i)->set_type(type);
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::change_type_cell_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, long long int type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
                if(change_type_sym_cell(dimension_a, dimension_b+i, dimension_c+j, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_type_cell_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, long long int type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(change_type_sym_cell(dimension_a+i, dimension_b, dimension_c+j, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_type_cell_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, long long int type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(change_type_sym_cell(dimension_a+i, dimension_b+j, dimension_c, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}


long long int Box::add_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z


    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz>fromz && atomz<untilz)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }

    }
    return 0;

}

long long int Box::add_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }
    }
    return 0;

}

long long int Box::add_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomy>fromy && atomy<untily)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }
    }
    return 0;

}

long long int Box::add_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }

    }
    return 0;

}

long long int Box::add_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)//x,z postion at z=0
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomx>fromx && atomx<untilx)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }

    }
    return 0;
}

long long int Box::add_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            box_atoms.at(i)->set_type(FREEPOSITION);
        }

    }
    return 0;
}

//Defining planes
//Ax+By+Cz+D=0

long long int Box::change_type_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_, long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

        if (distance<distance_)
        {
            box_atoms.at(i)->set_type(type);
        }
    }
    return 0;

}

//defining gap

long long int Box::change_type_gap(double x_gap,double y_gap, double z_gap,double radious_gap,long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
            atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
            atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
        {
            box_atoms.at(i)->set_type(type);
        }
    }
    return 0;
}

long long int Box::change_type_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere,long long int type)
{
    //(x-x0)^2+(y-y0)^2+(z-z0)^2=r^2
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
        {
            box_atoms.at(i)->set_type(type);
        }
    }
    return 0;

}

long long int Box::change_type_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, long long int type)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
        {
            box_atoms.at(i)->set_type(type);
        }
    }
    return 0;

}

long long int Box::change_type_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,long long int type)
{
    //Ax^2+By^2+Cz^2+Dxy+Exz+Fyz+Gx+Hy+Jz+K<1
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
        {
            box_atoms.at(i)->set_type(type);
        }
    }
    return 0;

}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------


//SELECTIVE BY TYPE IRREGULARITIES SELECTIVAS
//------------------------------------------------------------------------
//------------------------------------------------------------------------
long long int Box::change_type_sym_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c,long long int type)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
         long long int sizeatoms=cellatoms.size();

        #pragma omp parallel for
        for ( long long int i=0;i<sizeatoms;i++)
        {
            if (atom_type.compare(cellatoms.at(i)->get_atom_type())==0)
            {
                cellatoms.at(i)->set_type(type);
            }
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::change_type_cell_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, long long int type)
{
    bool all=true;

    #pragma omp parallel for shared (all)
    for ( long long int i=0; i<sidey; i++)
    {
        for ( long long int j=0; j<sidez; j++)
        {
            if(change_type_sym_cell(atom_type,dimension_a, dimension_b+i, dimension_c+j, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_type_cell_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, long long int type)
{
    bool all=true;

    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(change_type_sym_cell(atom_type,dimension_a+i, dimension_b, dimension_c+j, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_type_cell_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, long long int type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(change_type_sym_cell(atom_type,dimension_a+i, dimension_b+j, dimension_c, type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}


long long int Box::add_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz>fromz && atomz<untilz)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }
    }
    return 0;

}

long long int Box::add_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }

    }
    return 0;

}

long long int Box::add_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomy>fromy && atomy<untily)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }

    }
    return 0;
}

long long int Box::add_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }

    }
    return 0;

}

long long int Box::add_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomx>fromx && atomx<untilx)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }

    }
    return 0;
}

long long int Box::add_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                box_atoms.at(i)->set_type(FREEPOSITION);
            }
        }
    }
    return 0;
}

long long int Box::change_type_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_, long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

            if (distance<distance_)
            {
                box_atoms.at(i)->set_type(type);
            }
        }
    }
    return 0;

}

long long int Box::change_type_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap,long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
                atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
                atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
            {
                box_atoms.at(i)->set_type(type);
            }
        }
    }
    return 0;

}

long long int Box::change_type_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere,long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
            {
                box_atoms.at(i)->set_type(type);
            }
        }
    }
    return 0;

}

long long int Box::change_type_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, long long int type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
            {
                box_atoms.at(i)->set_type(type);
            }
        }
    }
    return 0;

}

long long int Box::change_type_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,long long int type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for (long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
            {
                box_atoms.at(i)->set_type(type);
            }
        }
    }
    return 0;

}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//REMOVE FROM SURFACE

long long int Box::remove_from_surface_cell(long long int dim_a, long long int dim_b, long long int dim_c)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX REMOVING FROM SURFACE!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
        long long int sizeatoms=cellatoms.size();

        for (long long int i=0;i<sizeatoms;i++)
        {
            rm_surface_atom(cellatoms.at(i)->get_id());
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::remove_from_surface_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;
    //paralelizarlo puede traer problemas porque podemos estar quitando simultaneamente de surface
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
                if(remove_from_surface_cell(dimension_a, dimension_b+i, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::remove_from_surface_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(remove_from_surface_cell(dimension_a+i, dimension_b, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::remove_from_surface_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(remove_from_surface_cell(dimension_a+i, dimension_b+j, dimension_c)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}



long long int Box::remove_from_surface_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz>fromz && atomz<untilz)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }

    }
    return 0;

}

long long int Box::remove_from_surface_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

     #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomy>fromy && atomy<untily)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }

    }
    return 0;

}

long long int Box::remove_from_surface_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)//x,z postion at z=0
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomx>fromx && atomx<untilx)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }

    }
    return 0;
}

long long int Box::remove_from_surface_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

     #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }

    }
    return 0;
}




long long int Box::remove_from_surface_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

        if (distance<distance_)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}

//defining gap

long long int Box::remove_from_surface_gap(double x_gap,double y_gap, double z_gap,double radious_gap)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
            atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
            atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;
}

long long int Box::remove_from_surface_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere)
{
    //(x-x0)^2+(y-y0)^2+(z-z0)^2=r^2
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K)
{
    //Ax^2+By^2+Cz^2+Dxy+Exz+Fyz+Gx+Hy+Jz+K<1
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
        {
            #pragma omp critical
            {
                rm_surface_atom(box_atoms.at(i)->get_id());
            }
        }
    }
    return 0;

}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//ADD TO SURFACE

long long int Box::add_to_surface_cell(long long int dim_a, long long int dim_b, long long int dim_c)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX ADDING TO SURFACE!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
        long long int sizeatoms=cellatoms.size();

        for (long long int i=0;i<sizeatoms;i++)
        {
            add_surface_atom(cellatoms.at(i));
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::add_to_surface_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;
    //paralelizarlo puede traer problemas porque podemos estar quitando simultaneamente de surface
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
                if(add_to_surface_cell(dimension_a, dimension_b+i, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::add_to_surface_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(add_to_surface_cell(dimension_a+i, dimension_b, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::add_to_surface_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(add_to_surface_cell(dimension_a+i, dimension_b+j, dimension_c)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}



long long int Box::add_to_surface_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz>fromz && atomz<untilz)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }

    }
    return 0;

}

long long int Box::add_to_surface_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

     #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomy>fromy && atomy<untily)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }

    }
    return 0;

}

long long int Box::add_to_surface_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)//x,z postion at z=0
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomx>fromx && atomx<untilx)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }

    }
    return 0;
}

long long int Box::add_to_surface_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

     #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }

    }
    return 0;
}




long long int Box::add_to_surface_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

        if (distance<distance_)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}

//defining gap

long long int Box::add_to_surface_gap(double x_gap,double y_gap, double z_gap,double radious_gap)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
            atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
            atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;
}

long long int Box::add_to_surface_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere)
{
    //(x-x0)^2+(y-y0)^2+(z-z0)^2=r^2
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K)
{
    //Ax^2+By^2+Cz^2+Dxy+Exz+Fyz+Gx+Hy+Jz+K<1
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
        {
            #pragma omp critical
            {
                add_surface_atom(box_atoms.at(i));
            }
        }
    }
    return 0;

}
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------



//ADD_TO_SURFACE SELECTIVAS
//------------------------------------------------------------------------
//------------------------------------------------------------------------
long long int Box::add_to_surface_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
         long long int sizeatoms=cellatoms.size();

        for ( long long int i=0;i<sizeatoms;i++)
        {
            if (atom_type.compare(cellatoms.at(i)->get_atom_type())==0)
            {
                add_surface_atom(cellatoms.at(i));
            }
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::add_to_surface_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;

    for ( long long int i=0; i<sidey; i++)
    {
        for ( long long int j=0; j<sidez; j++)
        {
            if(add_to_surface_cell(atom_type,dimension_a, dimension_b+i, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::add_to_surface_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;

    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(add_to_surface_cell(atom_type,dimension_a+i, dimension_b, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::add_to_surface_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(add_to_surface_cell(atom_type,dimension_a+i, dimension_b+j, dimension_c)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::add_to_surface_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz>fromz && atomz<untilz)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }

    }
    return 0;

}

long long int Box::add_to_surface_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomy>fromy && atomy<untily)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }

    }
    return 0;
}

long long int Box::add_to_surface_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }

    }
    return 0;

}

long long int Box::add_to_surface_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)
{

    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomx>fromx && atomx<untilx)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }

    }
    return 0;
}

long long int Box::add_to_surface_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;
}

long long int Box::add_to_surface_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

            if (distance<distance_)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
                atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
                atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}

long long int Box::add_to_surface_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for (long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
            {
                #pragma omp critical
                {
                    add_surface_atom(box_atoms.at(i));
                }
            }
        }
    }
    return 0;

}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

long long int Box::change_element_cell(long long int dim_a, long long int dim_b, long long int dim_c, string new_atom_type)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
        long long int sizeatoms=cellatoms.size();

        for (long long int i=0;i<sizeatoms;i++)
        {
            cellatoms.at(i)->set_atom_type(new_atom_type);
        }
        return 0;
    }
    else
    {
        return 1;
    }
}

long long int Box::change_element_yzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, string new_atom_type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidey; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
                if(change_element_cell(dimension_a, dimension_b+i, dimension_c+j, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xzplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, string new_atom_type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(change_element_cell(dimension_a+i, dimension_b, dimension_c+j, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xyplane(long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, string new_atom_type)
{
    bool all=true;
    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(change_element_cell(dimension_a+i, dimension_b+j, dimension_c, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xy_heli_dislocation(double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca, string new_atom_type)//x,y postion at z=0
{
     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz>fromz && atomz<untilz)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }

    }
    return 0;
}

long long int Box::change_element_xy_heli_dislocation(double x, double y, double angle_xz, double angle_yz, double radiousdisloca, string new_atom_type)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    //possitions x and y th dislocation have at atom's z

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomz*cos(angle_xz*PI/180.0);
        double disloy=y+atomz*cos(angle_yz*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;

}

long long int Box::change_element_xz_heli_dislocation(double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca, string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomy>fromy && atomy<untily)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;
}

long long int Box::change_element_xz_heli_dislocation(double x, double z, double angle_xy, double angle_zy, double radiousdisloca, string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double dislox=x+atomy*cos(angle_xy*PI/180.0);
        double disloz=z+atomy*cos(angle_zy*PI/180.0);

        if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }

    }
    return 0;
}

long long int Box::change_element_yz_heli_dislocation(double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca, string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
            && atomx>fromx && atomx<untilx)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }

    }
    return 0;
}

long long int Box::change_element_yz_heli_dislocation(double y, double z, double angle_yx, double angle_zx, double radiousdisloca, string new_atom_type)//x,z postion at z=0
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double disloy=y+atomx*cos(angle_yx*PI/180.0);
        double disloz=z+atomx*cos(angle_zx*PI/180.0);

        if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
            && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }

    }
    return 0;
}

//Definimos planos
//Ax+By+Cz+D=0
long long int Box::change_element_plane_distance(double A_plane, double B_plane, double C_plane, double D_plane, double distance_, string new_atom_type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

        if (distance<distance_)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;
}


long long int Box::change_element_gap(double x_gap,double y_gap, double z_gap,double radious_gap, string new_atom_type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
            atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
            atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;
}

long long int Box::change_element_sphere(double x_sphere,double y_sphere, double z_sphere,double radious_sphere, string new_atom_type)
{
     //(x-x0)^2+(y-y0)^2+(z-z0)^2=r^2
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;
}

long long int Box::change_element_ellipsoid(double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, string new_atom_type)
{

    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;

}


long long int Box::change_element_general_ellipsoid(double A, double B,double C,double D, double E, double F, double G,double H,double J, double K,string new_atom_type)
{
    //Ax^2+By^2+Cz^2+Dxy+Exz+Fyz+Gx+Hy+Jz+K<1
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        double atomx=box_atoms.at(i)->get_x();
        double atomy=box_atoms.at(i)->get_y();
        double atomz=box_atoms.at(i)->get_z();

        if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
        {
            box_atoms.at(i)->set_atom_type(new_atom_type);
        }
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


long long int Box::change_element_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c, string new_atom_type)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
         long long int sizeatoms=cellatoms.size();

        #pragma omp parallel for
        for ( long long int i=0;i<sizeatoms;i++)
        {
            if (atom_type.compare(cellatoms.at(i)->get_atom_type())==0)
            {
                cellatoms.at(i)->set_atom_type(new_atom_type);
            }
        }
        return 0;
    }
    else
    {
        return 1;
    }
}

long long int Box::change_element_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez, string new_atom_type)
{
    bool all=true;

    #pragma omp parallel for shared (all)
    for ( long long int i=0; i<sidey; i++)
    {
        for ( long long int j=0; j<sidez; j++)
        {
            if(change_element_cell(atom_type,dimension_a, dimension_b+i, dimension_c+j, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez, string new_atom_type)
{
    bool all=true;

    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(change_element_cell(atom_type,dimension_a+i, dimension_b, dimension_c+j, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey, string new_atom_type)
{
    bool all=true;

    #pragma omp parallel for shared (all)
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(change_element_cell(atom_type,dimension_a+i, dimension_b+j, dimension_c, new_atom_type)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::change_element_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca,string new_atom_type)//x,y postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz>fromz && atomz<untilz)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}

long long int Box::change_element_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca,string new_atom_type)//x,y postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}

long long int Box::change_element_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca,string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomy>fromy && atomy<untily)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }

    }
    return 0;
}

long long int Box::change_element_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca,string new_atom_type)//x,z postion at z=0
{

    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }

    }
    return 0;
}

long long int Box::change_element_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca,string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomx>fromx && atomx<untilx)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }

    }
    return 0;
}

long long int Box::change_element_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca,string new_atom_type)//x,z postion at z=0
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}

//Definimos planos
//Ax+By+Cz+D=0
long long int Box::change_element_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_, string new_atom_type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

            if (distance<distance_)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;

}

long long int Box::change_element_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap,string new_atom_type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
                atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
                atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}

long long int Box::change_element_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere,string new_atom_type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}

long long int Box::change_element_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z, string new_atom_type)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;

}

long long int Box::change_element_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K, string new_atom_type)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for (long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
            {
                box_atoms.at(i)->set_atom_type(new_atom_type);
            }
        }
    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////




//REVOME_FROM_SURFACE SELECTIVAS
//------------------------------------------------------------------------
//------------------------------------------------------------------------
long long int Box::remove_from_surface_cell(string atom_type,long long int dim_a, long long int dim_b, long long int dim_c)
{
    if (dim_a>dimension_x-1 || dim_b>dimension_y-1 || dim_c>dimension_z-1
        ||dim_a<0 || dim_b <0 || dim_c < 0)
    {
        cout<<"CELL OUT OF BOX!! "<<endl;
        return 1;
    }

    Cell * targetcell=select_cell_by_position(dim_a,dim_b,dim_c);
    if (targetcell!=nullptr)
    {
        vector<Atom*> cellatoms=targetcell->get_cell_atoms();
         long long int sizeatoms=cellatoms.size();

        for ( long long int i=0;i<sizeatoms;i++)
        {
            if (atom_type.compare(cellatoms.at(i)->get_atom_type())==0)
            {
                rm_surface_atom(cellatoms.at(i)->get_id());
            }
        }
        return 0;
    }
    else
    {
        return 1;
    }

}

long long int Box::remove_from_surface_yzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidey, long long int sidez)
{
    bool all=true;

    for ( long long int i=0; i<sidey; i++)
    {
        for ( long long int j=0; j<sidez; j++)
        {
            if(remove_from_surface_cell(atom_type,dimension_a, dimension_b+i, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::remove_from_surface_xzplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidez)
{
    bool all=true;

    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidez; j++)
        {
            if(remove_from_surface_cell(atom_type,dimension_a+i, dimension_b, dimension_c+j)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::remove_from_surface_xyplane(string atom_type,long long int dimension_a, long long int dimension_b, long long int dimension_c, long long int sidex, long long int sidey)
{
    bool all=true;
    for (long long int i=0; i<sidex; i++)
    {
        for (long long int j=0; j<sidey; j++)
        {
            if(remove_from_surface_cell(atom_type,dimension_a+i, dimension_b+j, dimension_c)!=0) all=false;
        }
    }
    if(all) return 0;
    else return 1;
}

long long int Box::remove_from_surface_xy_heli_dislocation(string atom_type,double x, double y, double fromz, double untilz, double angle_xz, double angle_yz, double radiousdisloca)//x,y postion at z=0
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz>fromz && atomz<untilz)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_xy_heli_dislocation(string atom_type,double x, double y, double angle_xz, double angle_yz, double radiousdisloca)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomz*cos(angle_xz*PI/180.0);
            double disloy=y+atomz*cos(angle_yz*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }

    }
    return 0;

}

long long int Box::remove_from_surface_xz_heli_dislocation(string atom_type,double x, double z, double fromy, double untily, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomy>fromy && atomy<untily)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }

    }
    return 0;
}

long long int Box::remove_from_surface_xz_heli_dislocation(string atom_type,double x, double z, double angle_xy, double angle_zy, double radiousdisloca)
{

     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double dislox=x+atomy*cos(angle_xy*PI/180.0);
            double disloz=z+atomy*cos(angle_zy*PI/180.0);

            if (atomx<dislox+radiousdisloca && atomx>dislox-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }

    }
    return 0;

}

long long int Box::remove_from_surface_yz_heli_dislocation(string atom_type,double y, double z, double fromx, double untilx, double angle_yx, double angle_zx, double radiousdisloca)
{

    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca
                && atomx>fromx && atomx<untilx)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }

    }
    return 0;
}

long long int Box::remove_from_surface_yz_heli_dislocation(string atom_type,double y, double z, double angle_yx, double angle_zx, double radiousdisloca)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double disloy=y+atomx*cos(angle_yx*PI/180.0);
            double disloz=z+atomx*cos(angle_zx*PI/180.0);

            if (atomy<disloy+radiousdisloca && atomy>disloy-radiousdisloca
                && atomz<disloz+radiousdisloca && atomz>disloz-radiousdisloca)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;
}

long long int Box::remove_from_surface_plane_distance(string atom_type,double A_plane, double B_plane, double C_plane, double D_plane, double distance_)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            double distance=Util::distance_point_to_plane(atomx,atomy,atomz,A_plane,B_plane, C_plane, D_plane);

            if (distance<distance_)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_gap(string atom_type,double x_gap,double y_gap, double z_gap,double radious_gap)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (atomx<x_gap+radious_gap && atomx>x_gap-radious_gap &&
                atomy<y_gap+radious_gap && atomy>y_gap-radious_gap &&
                atomz<z_gap+radious_gap && atomz>z_gap-radious_gap)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_sphere(string atom_type,double x_sphere,double y_sphere, double z_sphere,double radious_sphere)
{
     long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_sphere)*(atomx-x_sphere) + (atomy-y_sphere)*(atomy-y_sphere) + (atomz-z_sphere)*(atomz-z_sphere) < radious_sphere*radious_sphere)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_ellipsoid(string atom_type,double x_ellipsoid,double y_ellipsoid, double z_ellipsoid,double radious_x, double radious_y, double radious_z)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for ( long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if((atomx-x_ellipsoid)*(atomx-x_ellipsoid)/(radious_x * radious_x)+ (atomy-y_ellipsoid)*(atomy-y_ellipsoid)/(radious_y*radious_y) + (atomz-z_ellipsoid)*(atomz-z_ellipsoid)/(radious_z*radious_z) < 1.0)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}

long long int Box::remove_from_surface_general_ellipsoid(string atom_type,double A, double B,double C,double D, double E, double F, double G,double H,double J, double K)
{
    long long int boxsize=box_atoms.size();

    #pragma omp parallel for
    for (long long int i=0;i<boxsize;i++)
    {
        if (atom_type.compare(box_atoms.at(i)->get_atom_type())==0)
        {
            double atomx=box_atoms.at(i)->get_x();
            double atomy=box_atoms.at(i)->get_y();
            double atomz=box_atoms.at(i)->get_z();

            if (A*atomx*atomx + B*atomy*atomy + C*atomz*atomz + D*atomx*atomy + E*atomx*atomz+ F*atomy*atomz+ G*atomx+ H*atomy + J * atomz+ K< 1.0)
            {
                #pragma omp critical
                {
                    rm_surface_atom(box_atoms.at(i)->get_id());
                }
            }
        }
    }
    return 0;

}
