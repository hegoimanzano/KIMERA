#include "Atom.h"

Atom::Atom( long long int id_, double x_, double y_, double z_,long long int type_, bool insurface_)
{
    id = new  long long int;
    x = new double;
    y = new double;
    z = new double;
    atom_type = new string;
    type = new long long int;
    insurface = new bool;
    mass=new double;
    *id=id_;
    *x=x_;
    *y=y_;
    *z=z_;
    *type=type_;
    *insurface=insurface_;
    *mass=1.0;

}

Atom::Atom( long long int id_, double x_, double y_, double z_, string atom_type_, bool insurface_)
{
    id = new  long long int;
    x = new double;
    y = new double;
    z = new double;
    type = new long long int;
    atom_type = new string;
    insurface = new bool;
    mass=new double;
    *type=NORMAL;
    *id=id_;
    *x=x_;
    *y=y_;
    *z=z_;
    *atom_type=atom_type_;
    *insurface=insurface_;
    *mass=1.0;

}

Atom::~Atom()
{
    delete x;
    delete y;
    delete z;
    delete type;
    delete atom_type;
    delete id;
    delete insurface;
    delete mass;
}



void Atom::set_id( long long int id_)
{
    *id=id_;
}
 long long int Atom::get_id()
{
    return *id;
}

void Atom::set_mass(double mass_)
{
    *mass=mass_;
}

double Atom::get_mass()
{
    return *mass;
}


void Atom::set_x(double x_)
{
    *x=x_;
}
double Atom::get_x()
{
    return *x;
}


void Atom::set_y(double y_)
{
    *y=y_;
}
double Atom::get_y()
{
    return *y;
}

void Atom::set_z(double z_)
{
    *z=z_;
}
double Atom::get_z()
{
    return *z;
}

void Atom::set_type(long long int type_)
{
    *type=type_;
}
long long int Atom::get_type()
{
    return *type;
}

void Atom::set_atom_type(string atom_type_)
{
    *atom_type=atom_type_;
}
string Atom::get_atom_type()
{
    return *atom_type;
}

void Atom::set_insurface(bool insurface_)
{
    *insurface=insurface_;
}
bool Atom::get_insurface()
{
    return *insurface;
}

double Atom::get_distance(Atom *otheratom)
{
    double otheratomx=otheratom->get_x();
    double otheratomy=otheratom->get_y();
    double otheratomz=otheratom->get_z();

    double thisatomx=get_x();
    double thisatomy=get_y();
    double thisatomz=get_z();

    double distance=sqrt((otheratomx-thisatomx)*(otheratomx-thisatomx)+
                         (otheratomy-thisatomy)*(otheratomy-thisatomy)+
                         (otheratomz-thisatomz)*(otheratomz-thisatomz));


    return distance;
}



long long int Atom::add_neighbour(Atom *neighbour)
{
    neighbours.push_back(neighbour);
    return 0;
}

Atom* Atom::get_neighbour(long long int pos_)
{
    return neighbours.at(pos_);
}

long long int Atom::get_size_neighbour()
{
    return neighbours.size();
}

long long int Atom::rm_neighbour( long long int id)
{
    for ( long long int i=0; i<(long long int)neighbours.size(); i++)
    {
        if (id == neighbours.at(i)->get_id())
        {
            neighbours.erase(neighbours.begin()+i);
            return 0;
        }
    }
    return 1;
}

long long int Atom::rm_allneighbour()
{
    neighbours.clear();
    return 0;
}

long long int Atom::get_id_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_id();
}

bool Atom::get_insurface_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_insurface();
}

double Atom::get_x_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_x();
}

double Atom::get_y_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_y();
}

double Atom::get_z_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_z();
}
long long int Atom::get_type_neighbour(long long int pos_)
{
    return neighbours.at(pos_)->get_type();
}

vector<Atom*> Atom::get_neighbours()
{
    return neighbours;
}

vector<double> Atom::get_distances_to_neighbours()
{
    return distances_to_neighbours;
}
long long int Atom::add_distance_to_neighbours(double distance_)
{
    distances_to_neighbours.push_back(distance_);
    return distances_to_neighbours.size();
}
double Atom::get_distance_to_neighbours_by_pos(long long int pos_)
{
    return distances_to_neighbours.at(pos_);
}

long long int Atom::set_neigh_record() //lo estamos llamando dos veces.. que pasa si lo llamamos dos veces??!! pues que estamos metiendo el doble de vecinos.. \clap
{
    long long int sizeneighbours=neighbours.size();
    bool any=false;

    for (long long int i=0; i < sizeneighbours; i++)
    {
        long long int sizeneighrecord=neigh_record.size();
        for (long long int j=0; j< sizeneighrecord;j++)
        {

            if (neighbours.at(i)->get_atom_type().compare(neigh_record.at(j).get_type())==0 &&
                Util::isEqualDistances(distances_to_neighbours.at(i),neigh_record.at(j).get_distance()))
            {
                neigh_record.at(j).increment_number_by_1();
                any=true;
                break;
            }
        }
        if (!any)
        {
            //create a new record
            Record record1=Record(neighbours.at(i)->get_atom_type(),distances_to_neighbours.at(i));
            //add
            neigh_record.push_back(record1);
        }
        any=false;
    }
    return neigh_record.size();
}

vector<Record> Atom::get_neigh_record()
{
    return neigh_record;
}

long long int Atom::reduce_record_count_by_1(string type_, double distance_)
{
    long long int count;
    long long int sizeneighrecord=neigh_record.size();

    for (long long int i=0; i< sizeneighrecord;i++)
    {
        if (type_.compare(neigh_record.at(i).get_type())==0 &&
            Util::isEqualDistances(distance_,neigh_record.at(i).get_distance()))
        {
            count = neigh_record.at(i).reduce_number_by_1();
            if (count <0)
            {
                return 2;       // if size smaller than 0

            }
            return 0;
        }
    }
    return 1; //not finded
}

long long int Atom::rm_neighbour_and_record(long long int id_neigh)
{
    double distance_to_neigh=0.0;
    string type_neigh;

    for ( long long int i=0; i<(long long int)neighbours.size(); i++)
    {
        if (id_neigh == neighbours.at(i)->get_id())
        {
            distance_to_neigh=get_distances_to_neighbours().at(i);
            type_neigh=neighbours.at(i)->get_atom_type();
            reduce_record_count_by_1( type_neigh , distance_to_neigh);
            distances_to_neighbours.erase(distances_to_neighbours.begin()+i);
            neighbours.erase(neighbours.begin()+i);
            if (linked.size()>0)
            {
                linked.erase(linked.begin()+i);
                linked_type.erase(linked_type.begin()+i);
            }

            return 0;
        }
    }
    return 1;
}

long long int Atom::rm_neighbour_and_record_both_directions(long long int id_neigh)    //id of the atom which is dissolved
{
    double distance_to_neigh=0.0;
    string type_neigh;
    string this_type;
    long long int thisid=get_id();

    for ( long long int i=0; i<(long long int)neighbours.size(); i++)
    {
        if (id_neigh == neighbours.at(i)->get_id())
        {
            distance_to_neigh=get_distances_to_neighbours().at(i);
            type_neigh=neighbours.at(i)->get_atom_type();
            reduce_record_count_by_1( type_neigh , distance_to_neigh);
            distances_to_neighbours.erase(distances_to_neighbours.begin()+i);

            neighbours.at(i)->rm_neighbour_and_record(thisid);

            neighbours.erase(neighbours.begin()+i);
            if (linked.size()>0)
            {
                linked.erase(linked.begin()+i);
                linked_type.erase(linked_type.begin()+i);
            }
            return 0;
        }
    }
    return 1;
}

long long int Atom::get_number(string type_, double distance_)
{
    long long int size_neigh_record=neigh_record.size();

    for (long long int i=0;i<size_neigh_record;i++)
    {
        if (neigh_record.at(i).get_type().compare(type_)==0 && Util::isEqualDistances(distance_,neigh_record.at(i).get_distance()))
        {
            return neigh_record.at(i).get_number();
        }

    }
    cout<<"error in Atom::get_number(), it is possible that the definition of the events is not correct, otherwise, please contact the authors"<<endl;
    return 0;
}

vector<Atom*> Atom::get_affected()
{
    return affected;
}

long long int Atom::add_affected(Atom *affected_)
{
    affected.push_back(affected_);
    return affected.size();
}

long long int Atom::rm_affected( long long int id)
{
    for ( long long int i=0; i<(long long int)affected.size(); i++)
    {
        if (id == affected.at(i)->get_id())
        {
            affected.erase(affected.begin()+i);
            return 0;
        }
    }
    return 1;
}

vector<vector<Atom*>> Atom::get_linked()
{
    return linked;
}

void Atom::add_void_vector_in_linked()
{
    vector<Atom*> emptyvector;
    linked.push_back(emptyvector);
}

void Atom::remove_neighbours_with_out_linkeds()
{
/**
    for(long long int joder=0;joder<linked_type.size();joder++)
    {
        cout<<"neighbours.at(joder)->get_type():  "<<distances_to_neighbours.at(joder);
        cout<<endl;
        cout<<"linked.at(joder).size():  "<<linked.at(joder).size();
        cout<<endl;
        cout<<"linked_type.at(joder):  "<<linked_type.at(joder);
        cout<<endl;
        cout<<endl;
    }
    cout<<endl<<endl<<endl;;
    */

    long long int sizelinkedtype=linked_type.size();

    if (sizelinkedtype !=0)
    {
        if (linked.at(sizelinkedtype-1).size()==0)
        {
            distances_to_neighbours.erase(distances_to_neighbours.begin()+(sizelinkedtype-1));
            neighbours.erase(neighbours.begin()+(sizelinkedtype-1));
            linked.erase(linked.begin()+(sizelinkedtype-1));
            linked_type.erase(linked_type.begin()+(sizelinkedtype-1));
        }
    }
}

long long int  Atom::add_linked_to(Atom * linked_, long long int pos_)
{
    linked.at(pos_).push_back(linked_);
    return linked.at(pos_).size();
}

long long int  Atom::rm_linked(long long int pos_)
{
    linked.erase(linked.begin()+pos_);
    return linked.size();
}

bool Atom::is_bulk(vector<long long int> max_to_bulk, vector<string>types_, vector<double>distances_)
{
    long long int sizetypes=types_.size();
    long long int sizerecord=neigh_record.size();

    for (long long int i=0; i<sizetypes ;i++)
    {
        if (max_to_bulk.at(i)!=0)
        {
            for (long long int j=0;j<sizerecord;j++)
            {
                if (types_.at(i).compare(neigh_record.at(j).get_type())==0 && Util::isEqual(distances_.at(i),neigh_record.at(j).get_distance())
                     && get_number(types_.at(i),distances_.at(i))==max_to_bulk.at(i))
                {
                    return true;
                }
            }
        }
    }
    return false;
}


bool Atom::is_bulk_with_links(vector<long long int> max_to_bulk, vector<string> types_, vector<double> distances_, vector<Linked_neighbour> link_vector)
{
    long long int sizetypes=types_.size();

    for (long long int i=0; i<sizetypes ;i++)
    {
        if (max_to_bulk.at(i)!=0)
        {
            long long int typelink_from_event_def=link_vector.at(i).get_type();

            long long int truenumber= get_number_complex(types_.at(i),distances_.at(i),typelink_from_event_def);

            if (truenumber<max_to_bulk.at(i)) return false;
        }
    }
    return true;
}


long long int Atom::get_number_complex( string type_event_definition, double distance_event_definition, long long int typelink_from_event_def)
{

    long long int truenumber=0;


    //#pragma omp parallel for shared (truenumber)  //este aumenta mucho el tiempo de sim
    for(long long int s=0; s<(long long int)linked.size() ;s++)
    {
        long long int type = linked_type.at(s);
        double distance_to_neight=distances_to_neighbours.at(s);
        string type_neight=neighbours.at(s)->get_atom_type();

        if (type==NO_LINKED_TYPE && typelink_from_event_def==NO_LINKED_TYPE  && type_event_definition.compare(type_neight)==0  && Util::isEqualDistances(distance_event_definition,distance_to_neight) )
        {
            truenumber++;
        }
        if (type == LINKED_TYPE_NORMAL && typelink_from_event_def==LINKED_TYPE_NORMAL && type_event_definition.compare(type_neight)==0  && Util::isEqualDistances(distance_event_definition,distance_to_neight) )
        {
            //comprobamos que los linked que tiene el mismo id correspondiente, estan (todos ellos) tenemos que seleccionar por distancia y tipo
            vector<Atom*> inlinkeds= linked.at(s);
            bool enter=true;
            for (long long int l=0;l<(long long int)inlinkeds.size() ;l++)
            {
                if  (inlinkeds.at(l)->get_type()!=NORMAL)
                {
                    enter=false; break;
                }
            }
            if(enter) truenumber++;
        }
        if (type == LINKED_TYPE_DISSOLVED && typelink_from_event_def==LINKED_TYPE_DISSOLVED && type_event_definition.compare(type_neight)==0  && Util::isEqualDistances(distance_event_definition,distance_to_neight) )
        {
            //comprobamos que los linked que tiene el mismo id correspondiente, estan (todos ellos)
            vector<Atom*> inlinkeds= linked.at(s);
            bool enter=true;
            for (long long int l=0;l<(long long int)inlinkeds.size() ;l++)
            {
                if  (inlinkeds.at(l)->get_type()!=DISSOLVED && inlinkeds.at(l)->get_type()!=REMOVED)
                {
                    enter=false; break;
                }
            }
            if(enter) truenumber++;
        }
    }
    return truenumber;
}
