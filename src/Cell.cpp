#include "Cell.h"

Cell::Cell(long long int cell_id_)
{
    cell_id=cell_id_;
}

Cell::~Cell()
{

}

long long int Cell::get_cell_id()
{
    return cell_id;
}

long long int Cell::add_atom(Atom * incellatom)
{
    cell_atoms.push_back(incellatom);
    return 0;
}

long long int Cell::rm_atom(long long int id)
{
     long long int value;
    volatile bool flag=false;

    #pragma omp parallel for shared(flag)
    for ( long long int i=0; i<(long long int)cell_atoms.size(); i++)
    {
        if (flag) continue;
        if (id == cell_atoms.at(i)->get_id())
        {
            flag=true;
            value=i;
        }
    }
    if (flag) {cell_atoms.erase(cell_atoms.begin()+value);  return 0;}
    return 1;
}



long long int Cell::add_position(double x_, double y_, double z_,string type_, double prob_)
{
    long long int pos_id=get_cell_positions().size()+1;
    Position * auxposition= new Position(x_,y_,z_,pos_id,type_, prob_);
    cell_positions.push_back(auxposition);
    return 0;
}

 long long int Cell::compare_and_add(double x_, double y_, double z_,string type_, double prob_)
 {
     long long int cell_positions_size=cell_positions.size();
     long long int posid=get_cell_positions().size()+1;

     Position pos1=Position(x_,y_,z_,posid,type_,prob_);
     Position * pointerpos1=&pos1;
     for (long long int i=0 ; i< cell_positions_size;i++)
     {
         if (cell_positions.at(i)->compare_position(pointerpos1))
         {
             cell_positions.at(i)->add_type_prob(type_,prob_);

             if(cell_positions.at(i)->get_sum_probs()>1.000001) {cout<<"probability in position "<< i+1 <<"is greater than 1.0 "<< endl; exit(1);}
             return 1;
         }
     }

     add_position(x_,y_,z_,type_,prob_);
     return 0;
 }


long long int Cell::rm_position(long long int pos_id_)
{
     long long int cell_positions_size=cell_positions.size();

     long long int value;
    volatile bool flag=false;

    #pragma omp parallel for shared(flag)
    for ( long long int i=0; i<cell_positions_size; i++)
    {
        if(flag) continue;
        if (pos_id_ == cell_positions.at(i)->get_position_id())
        {
            Position * todelete = cell_positions.at(i);
            cell_removed_positions.push_back(todelete);
            value=i;
        }
    }
    if (flag) {cell_positions.erase(cell_positions.begin()+value); return 0;}
    return 1;
}

Position * Cell::get_pos_in(double x_,double y_, double z_)
{
     long long int cell_positions_size=cell_positions.size();

     for (long long int i=0 ; i< cell_positions_size;i++)
     {
         if (Util::isEqual(cell_positions.at(i)->get_x(), x_) &&
             Util::isEqual(cell_positions.at(i)->get_y(), y_) &&
             Util::isEqual(cell_positions.at(i)->get_z(), z_))
         {
             return cell_positions.at(i);
         }
     }
    return nullptr;
}


vector<Atom*> Cell::get_cell_atoms()
{
    return cell_atoms;
}

vector<Position*> Cell::get_cell_positions()
{
    return cell_positions;
}
vector<Position*> Cell::get_cell_removed_positions()
{
    return cell_removed_positions;
}

void Cell::set_cornerx(long long int cornerx_)
{
    cornerx=cornerx_;
}
long long int Cell::get_cornerx()
{
    return cornerx;
}

void Cell::set_cornery(long long int cornery_)
{
    cornery=cornery_;
}
long long int Cell::get_cornery()
{
    return cornery;
}

void Cell::set_cornerz(long long int cornerz_)
{
    cornerz=cornerz_;
}
long long int Cell::get_cornerz()
{
    return cornerz;
}

long long int Cell::add_cell_neighbour(Cell* cell)
{
    cell_neighbours.push_back(cell);
    return cell_neighbours.size();
}

vector<Cell*> Cell::get_cell_neighbour()
{
    return cell_neighbours;
}
