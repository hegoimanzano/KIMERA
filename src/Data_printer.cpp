#include "Data_printer.h"

Data_printer::Data_printer(string filename_)
{
        filename=filename_;

        previousremoved=0.0;
        previouspreviousremoved=0.0;

        previousdispersionx=0.0;
        previouspreviousdispersionx=0.0;

        previousdispersiony=0.0;
        previouspreviousdispersiony=0.0;

        previousdispersionz=0.0;
        previouspreviousdispersionz=0.0;

        previousalpha=0.0;
        previouspreviousalpha=0.0;

        previousgyradious=0.0;
        previouspreviousgyradious=0.0;

        previoustime=0.0;
        previousprevioustime=0.0;

        controldata=0.0;
        previouscontroldata=0.0;

        controlbox=0.0;
        previouscontrolbox=0.0;

        controllayer=0.0;
        previouscontrollayer=0.0;

        controlsurface=0.0;
        previouscontrolsurface=0.0;

        controlmean=0.0;
        previouscontrolmean=0.0;

        controldatastep=0;
        previouscontroldatastep=0;

        controlboxstep=0;
        previouscontrolboxstep=0;

        controllayerstep=0;
        previouscontrollayerstep=0;

        controlsurfacestep=0;
        previouscontrolsurfacestep=0;

        controlmeanstep=0;
        previouscontrolmeanstep=0;
}

Data_printer::~Data_printer()
{
    //dtor
}



long long int Data_printer::print_state(Box * box, Control *control) //Para guardar el estado... TODO Imprimiria la box y el control
{                                                             //Luego tendria que implementar el reader de estado
    return 0;                                               //tendria que meterle en el main e imprimirlo en caso de que no termine
}                                                           //la simulacion



long long int Data_printer::print_initial_box(Box * box, string name, bool affected_, bool links_)
{
    string fileinitialboxname=filename+"."+name+"kimerabox";

    ofstream fileinibox;

    fileinibox.open(fileinitialboxname,ios::app);

    long long int sizecells=box->get_box_cells().size();
    long long int size_atoms_in_cell;
    vector <Atom*> atoms_in_cell;

    vector <Cell*> cells=  box->get_box_cells();

    for (long long int i=0;i<sizecells;i++)
    {
        fileinibox<<"CELL "<<cells.at(i)->get_cell_id()<<endl;

        size_atoms_in_cell=cells.at(i)->get_cell_atoms().size();

        atoms_in_cell=cells.at(i)->get_cell_atoms();

         #pragma omp parallel for
        for (long long int j=0;j<size_atoms_in_cell;j++)
        {
                //if (atoms_in_cell.at(j)->get_type()==NORMAL || atoms_in_cell.at(j)->get_type()==INSOLUBLE)  //esto ahora que hay linkeds deberia no estar
                //{
                    #pragma omp critical
                    {
                       fileinibox<<atoms_in_cell.at(j)->get_id()<<"   "<< atoms_in_cell.at(j)->get_atom_type() <<
                        "   "<< atoms_in_cell.at(j)->get_type() << "   "<< atoms_in_cell.at(j)->get_x() <<
                        "   "<< atoms_in_cell.at(j)->get_y() <<"   "<< atoms_in_cell.at(j)->get_z()<<
                        "   "<< atoms_in_cell.at(j)->get_insurface()<<endl ;
                    }
                //}
        }
    }

    fileinibox<<endl<<"NEIGHBORS "<<endl;

    for (long long int i=0;i<sizecells;i++)
    {
        size_atoms_in_cell=cells.at(i)->get_cell_atoms().size();

        atoms_in_cell=cells.at(i)->get_cell_atoms();

        #pragma omp parallel for
        for (long long int j=0;j<size_atoms_in_cell;j++)
        {
            vector <Atom*> neigh_from_atom=atoms_in_cell.at(j)->get_neighbours();

            if(neigh_from_atom.size()!=0)
            {
                if (atoms_in_cell.at(j)->get_type()==NORMAL || atoms_in_cell.at(j)->get_type()==INSOLUBLE)
                {
                    #pragma omp critical
                    {
                        fileinibox<<atoms_in_cell.at(j)->get_id()<<"   ";

                        for ( long long int k=0; k<(long long int)neigh_from_atom.size();k++)
                        {
                            fileinibox<<neigh_from_atom.at(k)->get_id()<<"   ";
                            fileinibox<<atoms_in_cell.at(j)->get_distances_to_neighbours().at(k)<<"   ";
                        }
                        fileinibox<<endl;
                    }
                }
            }
        }
    }

    if(affected_)
    {

        fileinibox<<endl<<"AFFECTED "<<endl;

        for (long long int i=0;i<sizecells;i++)
        {
            size_atoms_in_cell=cells.at(i)->get_cell_atoms().size();

            atoms_in_cell=cells.at(i)->get_cell_atoms();

            #pragma omp parallel for
            for (long long int j=0;j<size_atoms_in_cell;j++)
            {
                vector <Atom*> affected_from_atom=atoms_in_cell.at(j)->get_affected();

                if (affected_from_atom.size()!=0)
                {
                    if (atoms_in_cell.at(j)->get_type()==NORMAL || atoms_in_cell.at(j)->get_type()==INSOLUBLE)
                    {
                        #pragma omp critical
                        {
                            fileinibox<<atoms_in_cell.at(j)->get_id()<<"   ";

                            for ( long long int k=0; k<(long long int)affected_from_atom.size();k++)
                            {
                                if(affected_from_atom.at(k)->get_type()==NORMAL || affected_from_atom.at(k)->get_type()==INSOLUBLE)
                                {
                                   fileinibox<<affected_from_atom.at(k)->get_id()<<"   ";
                                }
                            }
                            fileinibox<<endl;
                        }
                    }
                }
            }
        }
    }

    fileinibox<<endl<<"LINKED "<<endl;

    for (long long int i=0;i<sizecells;i++)
    {
        size_atoms_in_cell=cells.at(i)->get_cell_atoms().size();

        atoms_in_cell=cells.at(i)->get_cell_atoms();

        #pragma omp parallel for
        for (long long int j=0;j<size_atoms_in_cell;j++)
        {

            vector <vector<Atom*>> linked_from_atom=atoms_in_cell.at(j)->get_linked();
            vector <long long int> type_linked_from_atom=atoms_in_cell.at(j)->get_linked_type();

            if (linked_from_atom.size()!=0)
            {
                if (atoms_in_cell.at(j)->get_type()==NORMAL || atoms_in_cell.at(j)->get_type()==INSOLUBLE)
                {
                    #pragma omp critical
                    {
                        fileinibox<<atoms_in_cell.at(j)->get_id();


                        for ( long long int k=0; k<(long long int)type_linked_from_atom.size();k++)
                        {
                            fileinibox<<"  L  "<<type_linked_from_atom.at(k);
                            for (long long int s=0;s<(long long int)linked_from_atom.at(k).size();s++)
                            {
                                fileinibox<<"   "<<linked_from_atom.at(k).at(s)->get_id();
                            }

                         }
                         fileinibox<<endl;
                    }
                }
            }
        }
    }


    return 0;
}



long long int Data_printer::print_all_data(Box * box, Control *control,Tracker * tracker)
{

    double time= control->get_time();

    double deltat=1.0e-30;
    double targettime=control->get_targettime();
    double boxsteps=(double)control->get_boxsteps();

    if (control->get_boxframe())
    {
        if (control->get_step()==0){print_complete_system(box,control);}
        if (time-controlbox>targettime/boxsteps && time-previouscontrolbox>deltat)
        {
            print_complete_system(box,control);
            previouscontrolbox=time;
            controlbox+=targettime/boxsteps;
        }
    }

    double datastep=(double)control->get_datasteps();

    if (control->get_dataanalisys())
    {
        if (control->get_step()==0){print_data(box,control,tracker);}
        if(time-controldata>targettime/datastep && time-previouscontroldata>deltat)
        {
            print_data(box,control,tracker);
            previouscontrolbox=time;
            controldata+=targettime/datastep;
        }
    }

    double layerstep=(double)control->get_layersteps();

    if (control->get_layerxanalysis() || control->get_layeryanalysis()|| control->get_layerzanalysis())
    {
        if (control->get_step()==0 && control->get_layerxanalysis()){print_xlayer(box,control);}
        if (control->get_step()==0 && control->get_layeryanalysis()){print_ylayer(box,control);}
        if (control->get_step()==0 && control->get_layerzanalysis()){print_zlayer(box,control);}

        if (time-controllayer>targettime/layerstep && time-previouscontrollayer>deltat)
        {
            if(control->get_layerxanalysis()) {print_xlayer(box,control);}
            if(control->get_layeryanalysis()) {print_ylayer(box,control);}
            if(control->get_layerzanalysis()) {print_zlayer(box,control);}
            previouscontrollayer=time;
            controllayer+=targettime/layerstep;
        }
    }

    double surfacestep=(double)control->get_surfacesteps();

    if(control->get_surfaceframe())
    {
        if (control->get_step()==0) {print_surface(box,control);}

        if(time-controlsurface>targettime/surfacestep && time- previouscontrolsurface>deltat)
        {
            print_surface(box,control);
            previouscontrolsurface=time;
            controlsurface+=targettime/surfacestep;
        }
    }

    double meanstep= (double)control->get_meansteps();

    if (control->get_meandiscoord())
    {
        if (control->get_step()==0) {print_mean_dissolved_coordination(control,tracker);}

        if(time-controlmean>targettime/meanstep && time - previouscontrolmean> deltat)
        {
            print_mean_dissolved_coordination(control,tracker);
            previouscontrolmean=time;
            controlmean+=targettime/meanstep;
        }
    }

    return 0;
}

long long int Data_printer::print_all_data_step(Box * box, Control *control,Tracker * tracker)
{

    long long int step= control->get_step();

    //double deltat=1.0e-30; //la diferencia deber 1 paso
    long long int targetsteps=control->get_targetsteps();
    long long int boxsteps=control->get_boxsteps();

    if (control->get_boxframe())
    {
        if (step==0){print_complete_system(box,control);}
        if (step-controlboxstep>targetsteps/boxsteps && step-previouscontrolboxstep>1)
        {
            print_complete_system(box,control);
            previouscontrolboxstep=step;
            controlboxstep+=targetsteps/boxsteps;
        }
    }

    long long int datastep=control->get_datasteps();

    if (control->get_dataanalisys())
    {
        if (step==0){print_data(box,control,tracker);}
        if(step-controldatastep>targetsteps/datastep && step-previouscontroldatastep>1)
        {
            print_data(box,control,tracker);
            previouscontrolboxstep=step;
            controldatastep+=targetsteps/datastep;
        }
    }

    long long int layerstep=control->get_layersteps();

    if (control->get_layerxanalysis() || control->get_layeryanalysis()|| control->get_layerzanalysis())
    {
        if (step==0 && control->get_layerxanalysis()){print_xlayer(box,control);}
        if (step==0 && control->get_layeryanalysis()){print_ylayer(box,control);}
        if (step==0 && control->get_layerzanalysis()){print_zlayer(box,control);}

        if (step-controllayerstep>targetsteps/layerstep && step-previouscontrollayerstep>1)
        {
            if(control->get_layerxanalysis()) {print_xlayer(box,control);}
            if(control->get_layeryanalysis()) {print_ylayer(box,control);}
            if(control->get_layerzanalysis()) {print_zlayer(box,control);}
            previouscontrollayerstep=step;
            controllayerstep+=targetsteps/layerstep;
        }
    }

    long long int surfacestep=control->get_surfacesteps();

    if(control->get_surfaceframe())
    {
        if (step==0) {print_surface(box,control);}

        if(step-controlsurfacestep>targetsteps/surfacestep && step- previouscontrolsurfacestep>1)
        {
            print_surface(box,control);
            previouscontrolsurfacestep=step;
            controlsurfacestep+=targetsteps/surfacestep;
        }
    }

    long long int meanstep= control->get_meansteps();

    if (control->get_meandiscoord())
    {
        if (step==0) {print_mean_dissolved_coordination(control,tracker);}

        if(step-controlmeanstep>targetsteps/meanstep && step - previouscontrolmeanstep> 1)
        {
            print_mean_dissolved_coordination(control,tracker);
            previouscontrolmeanstep=step;
            controlmeanstep+=targetsteps/meanstep;
        }
    }

    return 0;
}

long long int Data_printer::print_data(Box * box, Control *control, Tracker * tracker)
{

    long long int sizebox=box->get_box_atoms().size();
    long long int size_typesinboxini;

    double time=control->get_time();
    long long int removed_atoms=box->get_box_removed_atoms().size();
    long long int total_atoms=removed_atoms+sizebox;
    double atoms_fraction=(double)removed_atoms/ (double)total_atoms;
    long long int events_accepted=control->get_eventsaccepted();
    double surface_dispersion_x;
    double surface_dispersion_y;
    double surface_dispersion_z;
    double gyradious;
    double alpha;

    double derivadaremoved=0.0;
    double derivadafraction=0.0;
    double derivadadispersionx=0.0;
    double derivadadispersiony=0.0;
    double derivadadispersionz=0.0;
    double derivadaalpha=0.0;
    double derivadagyradious=0.0;

    string filedataname=filename+".data";

    ofstream filedata;

    filedata.open(filedataname,ios::app);

    if (control->get_step()==0)
    {
        //detect the atoms type inside the box

        if( set_initial_type_number(box) !=0) {cout<<"There is an error in the fraction data"<<endl;}
        size_typesinboxini=typesinbox_initial.size();

        filedata<<"time \t\t ";

        for (long long int i=0;i<size_typesinboxini;i++)
        {
            filedata<<"removed_atoms_"<<typesinbox_initial.at(i)<<"\t\t";
            filedata<<"fraction_atoms_"<<typesinbox_initial.at(i)<<"\t\t";
        }

        filedata<< "removed_atoms_total\t\tatoms_fraction_total\t\tevents_accepted";
        if (control->get_plane_y() && control->get_plane_z()) filedata<<"\t\tsurface_dispersion_a";
        if (control->get_plane_x() && control->get_plane_z()) filedata<<"\t\tsurface_dispersion_b";
        if (control->get_plane_x() && control->get_plane_y()) filedata<<"\t\tsurface_dispersion_c";
        if (control->get_grain()) filedata<<"\t\tgyradious\t\talpha";
        filedata<<"\t\tdremoved_atoms/dt\t\tdatoms_fraction/dt";
        if (control->get_plane_y() && control->get_plane_z()) filedata<<"\t\tdsurface_dispersion_a/dt";
        if (control->get_plane_x() && control->get_plane_z()) filedata<<"\t\tdsurface_dispersion_b/dt";
        if (control->get_plane_x() && control->get_plane_y()) filedata<<"\t\tdsurface_dispersion_c/dt";
        if (control->get_grain())
        {
            filedata<<"\t\tdgyradious/dt\t\tdalpha/dt";
             gyradious=get_grain_gyradious(box);
             control->set_gyradiousini(gyradious);
             //alpha=get_grain_alpha(box,control);
             alpha=get_grain_alpha(box,control,gyradious);
        }

        filedata<<endl;
    }

    filedata<<time<<"\t\t";
    long long int number;
    size_typesinboxini=typesinbox_initial.size();
    for (long long int i=0;i<size_typesinboxini;i++)
    {
        number=tracker->get_number(typesinbox_initial.at(i));
        filedata<<number<<"\t\t";
        filedata<<(double)number/(double)numberinbox_initial.at(i)<<"\t\t";
    }
    filedata<<removed_atoms<<"\t\t"<<atoms_fraction<< "\t\t"<< events_accepted << "\t\t"  ;
    if (control->get_plane_y() && control->get_plane_z()) {surface_dispersion_x=get_surface_dispersion_x(box); filedata<< surface_dispersion_x<<"\t\t";}
    if (control->get_plane_x() && control->get_plane_z()) {surface_dispersion_y=get_surface_dispersion_y(box); filedata<< surface_dispersion_y<<"\t\t";}
    if (control->get_plane_x() && control->get_plane_y()) {surface_dispersion_z=get_surface_dispersion_z(box); filedata<< surface_dispersion_z<<"\t\t";}
    if (control->get_grain())
    {
        gyradious=get_grain_gyradious(box);
        //alpha=get_grain_alpha(box,control);
        alpha=get_grain_alpha(box,control,gyradious);
        filedata<<gyradious<<"\t\t"<<alpha;
    }

    if (control->get_step()<2) filedata<<endl;

    if (control->get_step()>=2)
    {
        //derivatives
        derivadaremoved=((double)removed_atoms-(double)previouspreviousremoved)/(time-previousprevioustime);
        derivadafraction=((double)removed_atoms-(double)previouspreviousremoved)/(total_atoms*(time-previousprevioustime));
        if (control->get_plane_y() && control->get_plane_z())
        {
            derivadadispersionx=(surface_dispersion_x-previouspreviousdispersionx)/(time-previousprevioustime);
        }
        if (control->get_plane_x() && control->get_plane_z())
        {
            derivadadispersiony=(surface_dispersion_y-previouspreviousdispersiony)/(time-previousprevioustime);
        }
        if (control->get_plane_x() && control->get_plane_y())
        {
            derivadadispersionz=(surface_dispersion_z-previouspreviousdispersionz)/(time-previousprevioustime);
        }
        if (control->get_grain())
        {
            derivadaalpha=(alpha-previouspreviousalpha)/(time-previousprevioustime);
            derivadagyradious=(gyradious-previouspreviousgyradious)/(time-previousprevioustime);
        }


        filedata<<"\t\t"<<derivadaremoved<<"\t\t"<<derivadafraction;
        if (control->get_plane_y() && control->get_plane_z()) filedata<<"\t\t"<<derivadadispersionx;
        if (control->get_plane_x() && control->get_plane_z()) filedata<<"\t\t"<<derivadadispersiony;
        if (control->get_plane_x() && control->get_plane_y()) filedata<<"\t\t"<<derivadadispersionz;
        if (control->get_grain()) filedata<<"\t\t"<<derivadagyradious<<"\t\t"<<derivadaalpha;
        filedata<<endl;
    }

    previouspreviousremoved=previousremoved;
    previouspreviousalpha=previousalpha;
    previouspreviousgyradious=previousgyradious;
    previouspreviousdispersionx=previousdispersionx;
    previouspreviousdispersiony=previousdispersiony;
    previouspreviousdispersionz=previousdispersionz;
    previousprevioustime=previoustime;

    previousremoved=removed_atoms;
    previousalpha=alpha;
    previousgyradious=gyradious;
    previousdispersionx=surface_dispersion_x;
    previousdispersiony=surface_dispersion_y;
    previousdispersionz=surface_dispersion_z;
    previousprevioustime=previoustime;

    filedata.close();

    return 0;
}

long long int Data_printer::print_complete_system(Box * box, Control *control)
{

    ofstream filecompletesystem;

    string filecompletesystemname;

    filecompletesystemname=filename+".box";

    filecompletesystem.open(filecompletesystemname,ios::app);

    filecompletesystem<<"ITEM: TIMESTEP" <<endl;
    filecompletesystem<<control->get_step() <<"      " <<setprecision(12)<<control->get_time() <<endl;

    filecompletesystem<<"ITEM: NUMBER OF ATOMS" <<endl;
    filecompletesystem<<box->get_box_atoms().size()<<endl;

    filecompletesystem<<"ITEM: BOX BOUNDS xy xz yz"<<endl;
    filecompletesystem<< control->get_box_boundxl0() <<"  "<< control->get_box_boundxl()<<"  "<< control->get_tilt_factors_xy() <<endl;
    filecompletesystem<< control->get_box_boundyl0() <<"  "<< control->get_box_boundyl()<<"  "<< control->get_tilt_factors_xz() <<endl;
    filecompletesystem<< control->get_box_boundzl0() <<"  "<< control->get_box_boundzl()<<"  "<< control->get_tilt_factors_yz() <<endl;

    filecompletesystem<<"ITEM: ATOMS id type i2 coord x y z"<< endl;
    long long int limit=box->get_box_atoms().size();

    #pragma omp parallel for
    for (long long int i=0; i<limit;i++)
    {
        Atom *selectedatom=box->get_box_atoms().at(i);

        #pragma omp critical
        {
            filecompletesystem<<selectedatom->get_id()<<" "<<selectedatom->get_atom_type()<<" "<<selectedatom->get_type()<<" "
            <<selectedatom->get_size_neighbour()<<" "<<selectedatom->get_x()<<" "<<selectedatom->get_y()
            <<" "<<selectedatom->get_z()<<endl;
        }
    }

    filecompletesystem.close();

    return 0;
}

long long int Data_printer::print_surface(Box * box, Control *control)
{
    ofstream filesurface;

    string filesurfacename;

    filesurfacename=filename+".surface";

    filesurface.open(filesurfacename,ios::app);

    filesurface<<"ITEM: TIMESTEP" <<endl;
    filesurface<<control->get_step() <<"      " <<setprecision(12)<<control->get_time() <<endl;

    filesurface<<"ITEM: NUMBER OF ATOMS" <<endl;
    filesurface<<box->get_surface_size()<<endl;

    filesurface<<"ITEM: BOX BOUNDS xy xz yz"<<endl;
    filesurface<< control->get_box_boundxl0() <<"  "<< control->get_box_boundxl()<<"  "<< control->get_tilt_factors_xy() <<endl;
    filesurface<< control->get_box_boundyl0() <<"  "<< control->get_box_boundyl()<<"  "<< control->get_tilt_factors_xz() <<endl;
    filesurface<< control->get_box_boundzl0() <<"  "<< control->get_box_boundzl()<<"  "<< control->get_tilt_factors_yz() <<endl;

    filesurface<<"ITEM: ATOMS id type coord x y z"<< endl;
    long long int limit=box->get_surface_size();
    #pragma omp parallel for
    for (long long int i=0; i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        #pragma omp critical
        {
        filesurface<<selectedatom->get_id()<<" "<<selectedatom->get_atom_type()<<" "
        <<selectedatom->get_size_neighbour()<<" "<<selectedatom->get_x()<<" "<<selectedatom->get_y()
        <<" "<<selectedatom->get_z()<<endl;
        }
    }

    filesurface.close();

    return 0;
}

long long int Data_printer::print_coord(Box * box, Control *control)  //Este creo que no tendria sentido
{
    /**
    ofstream filecoord;

    char filecoordname[24];

    strcpy(filecoordname,filename);

    strncat(filecoordname,".coord",24);

    filecoord.open(filecoordname,ios::app);

    //primero determinar cual es la coordinacion inicial maxima (se guarda en control)
    if (control->get_step()==0)
    {
         long long int maxcoordini=0;
         long long int coord=0;
        long long int limit=box->get_box_atoms().size();

        for (long long int i=0;i<limit;i++)
        {
            Atom *selectedatom=box->get_box_atoms().at(i);
            coord=selectedatom->get_size_neighbour();
            if (coord>maxcoordini) maxcoordini=coord;
        }
        control->set_maxcoordini(maxcoordini);

        //Aprovechamos y metemos la cabecera

        filecoord<<"time         ";

        for (long long int i=0;i<=maxcoordini;i++)
        {
            filecoord<<"coordination"<<i<<"   ";
        }
        filecoord<<endl;
    }

    long long int coordination[control->get_maxcoordini()+1];

    for (long long int i=0;i<=control->get_maxcoordini();i++)
    {
        coordination[i]=0;
    }

    long long int limit= box->get_box_atoms().size();

     long long int coord=0;

    for (long long int i=0; i<limit;i++)
    {
        Atom *selectedatom=box->get_box_atoms().at(i);
        coord=selectedatom->get_size_neighbour();
        coordination[coord]++;
    }

    //ponemos tambien los de coordinacion 0 (atomos disueltos)

    coordination[0]=box->get_box_removed_atoms().size();

    filecoord<<control->get_time()<<"            ";

    for (long long int i=0;i<=control->get_maxcoordini();i++)
    {
        filecoord<<coordination[i]<<"           ";
    }

    filecoord<<endl;

    filecoord.close();
*/
return 0;
}

long long int Data_printer::print_xlayer(Box * box, Control *control) // amount of atoms in each layer x case
{
    ofstream filelayerx;

    string filelayername;

    filelayername=filename+".alayer";

    filelayerx.open(filelayername,ios::app);



    long long int totalxlayers=control->get_dimension_x();

    vector<long long int> layer;


    for (long long int i=0;i<totalxlayers;i++)
    {
        layer.push_back(0);
    }


    //long long int sizecell=0;


    for (long long int i =0 ; i< totalxlayers;i++)
    {
        long long int counter=0;
        long long int limit=box->get_dimension_y();
        #pragma omp parallel for
        for (long long int j=0; j<limit;j++)
        {
            for (long long int k=0; k<box->get_dimension_z();k++)
            {
                Cell * targetcell= box->select_cell_by_position(i,j,k);
                //we see the amount of atoms in the cell
                long long int sizecell=targetcell->get_cell_atoms().size();

                for (long long int h=0;h<sizecell;h++)
                {
                    if (targetcell->get_cell_atoms().at(h)->get_type()==NORMAL || targetcell->get_cell_atoms().at(h)->get_type()==INSOLUBLE)    //vamos a suponer que todos los atomos los hemos definido con type 1
                    {
                        #pragma omp critical
                        {
                            counter++;
                        }
                    }
                }
            }

        }
        layer.at(i)=counter;
    }

    if (control->get_step()==0) //prlong long int the header
    {

        filelayerx<<"time         ";

        for (long long int i=1;i<=totalxlayers;i++)
        {
            filelayerx<<"layer_a"<<i<<"   ";
        }

        filelayerx<<endl;

    }

    filelayerx<<control->get_time();

    for (long long int i=0;i<totalxlayers;i++)
    {
        filelayerx<<"     "<<layer.at(i);
    }
    filelayerx<<endl;

    filelayerx.close();


 return 0;
}

long long int Data_printer::print_ylayer(Box * box, Control *control) // amount of atoms in each layer y case
{

    ofstream filelayery;

    string filelayername;

    filelayername=filename+".blayer";

    filelayery.open(filelayername,ios::app);


    long long int totalylayers=control->get_dimension_y();

    vector<long long int> layer;

    for (long long int i=0;i<totalylayers;i++)
    {
        layer.push_back(0);
    }

    long long int sizecell=0;


    for (long long int i =0 ; i< totalylayers;i++)
    {
        long long int counter=0;
        long long int limit=box->get_dimension_x();
        #pragma omp parallel for
        for (long long int j=0; j<limit;j++)
        {
            for (long long int k=0; k<box->get_dimension_z();k++)
            {
                Cell * targetcell= box->select_cell_by_position(j,i,k);

                sizecell=targetcell->get_cell_atoms().size();

                for (long long int h=0;h<sizecell;h++)
                {
                    if (targetcell->get_cell_atoms().at(h)->get_type()==NORMAL || targetcell->get_cell_atoms().at(h)->get_type()==INSOLUBLE)
                    {
                        #pragma omp critical
                        {
                            counter++;
                        }
                    }

                }

            }
        }
        layer.at(i)=counter;
    }

    if (control->get_step()==0) //prlong long int the header
    {

        filelayery<<"time         ";

        for (long long int i=1;i<=totalylayers;i++)
        {
            filelayery<<"layer_b"<<i<<"   ";
        }

        filelayery<<endl;

    }

    filelayery<<control->get_time();

    for (long long int i=0;i<totalylayers;i++)
    {
        filelayery<<"     "<<layer.at(i);
    }
    filelayery<<endl;

    filelayery.close();


 return 0;
}

long long int Data_printer::print_zlayer(Box * box, Control *control) // amount of atoms in each layer z case
{
    ofstream filelayerz;

    string filelayername;

    filelayername=filename+".clayer";

    filelayerz.open(filelayername,ios::app);

    long long int totalzlayers=control->get_dimension_z();

    vector<long long int> layer;


    for (long long int i=0;i<totalzlayers;i++)
    {
        layer.push_back(0);
    }

    long long int sizecell=0;


    for (long long int i =0 ; i< totalzlayers;i++)
    {
        long long int counter=0;
        long long int limit=box->get_dimension_x();
        //elegimos la celda
        #pragma omp parallel for
        for (long long int j=0; j<limit;j++)
        {
            for (long long int k=0; k<box->get_dimension_y();k++)
            {
                Cell * targetcell= box->select_cell_by_position(j,k,i);

                //miramos la cantidad de atomos dentro de la celda

                sizecell=targetcell->get_cell_atoms().size();

                for (long long int h=0;h<sizecell;h++)
                {
                    if (targetcell->get_cell_atoms().at(h)->get_type()==NORMAL || targetcell->get_cell_atoms().at(h)->get_type()==INSOLUBLE)    //vamos a suponer que todos los atomos los hemos definido con type 1
                    {
                        #pragma omp critical
                        {
                            counter++;
                        }
                    }

                }
            }
        }
        layer.at(i)=counter;
    }

    if (control->get_step()==0) //prlong long int the header
    {

        filelayerz<<"time         ";

        for (long long int i=1;i<=totalzlayers;i++)
        {
            filelayerz<<"layer_c"<<i<<"   ";
        }

        filelayerz<<endl;

    }

    filelayerz<<control->get_time();

    for (long long int i=0;i<totalzlayers;i++)
    {
        filelayerz<<"     "<<layer.at(i);
    }
    filelayerz<<endl;

    filelayerz.close();


 return 0;
}

double Data_printer::get_surface_dispersion_x(Box * box)
{

    double meanhigh=0.0;
    double dispersionh=0.0;

    long long int limit=box->get_surface_size();

    //#pragma omp parallel for shared (meanhigh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        meanhigh+=selectedatom->get_x();
    }

    meanhigh=meanhigh/(double)limit;

    //calculate height dispersion

    //#pragma omp parallel for shared (dispersionh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        double x=selectedatom->get_x();
        dispersionh+=(x-meanhigh)*(x-meanhigh);
    }

    dispersionh=dispersionh/(double)limit;


    return dispersionh;

}
double Data_printer::get_surface_dispersion_y(Box * box)
{

    double meanhigh=0.0;
    double dispersionh=0.0;

    long long int limit=box->get_surface_size();

    //#pragma omp parallel for shared (meanhigh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        meanhigh+=selectedatom->get_y();
    }

    meanhigh=meanhigh/(double)limit;

    //calculate height dispersion
    //#pragma omp parallel for shared (dispersionh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        double y=selectedatom->get_y();
        dispersionh+=(y-meanhigh)*(y-meanhigh);
    }

    dispersionh=dispersionh/(double)limit;


    return dispersionh;

}
double Data_printer::get_surface_dispersion_z(Box * box)
{

    double meanhigh=0.0;
    double dispersionh=0.0;

    long long int limit=box->get_surface_size();

    //#pragma omp parallel for shared (meanhigh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        meanhigh+=selectedatom->get_z();
    }

    meanhigh=meanhigh/(double)limit;

    //calculate height dispersion

    //#pragma omp parallel for shared (dispersionh)
    for(long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_surface_atoms().at(i);
        double z=selectedatom->get_z();
        dispersionh+=(z-meanhigh)*(z-meanhigh);
    }

    dispersionh=dispersionh/(double)limit;


    return dispersionh;

}

double Data_printer::get_grain_gyradious(Box * box)
{
    //first get the mass center, later the momentum
    double summassx=0.0;
    double summassy=0.0;
    double summassz=0.0;

    double totalmass=0.0;

    long long int limit=box->get_box_atoms().size();

    //#pragma omp parallel for shared (summassx,summassy,summassz,totalmass)
    for (long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_box_atoms().at(i);

        double selectedatommass=selectedatom->get_mass(); //necesitamos comando para meter la masa

        summassx+=selectedatom->get_x()*selectedatommass;
        summassy+=selectedatom->get_y()*selectedatommass;
        summassz+=selectedatom->get_z()*selectedatommass;

        totalmass+=selectedatommass;
    }

    double rcmx=summassx/totalmass;
    double rcmy=summassy/totalmass;
    double rcmz=summassz/totalmass;

    double inertiamoment=0.0;

    //#pragma omp parallel for shared (inertiamoment)
    for (long long int i=0;i<limit;i++)
    {
        Atom *selectedatom=box->get_box_atoms().at(i);

        double selectedatomx = selectedatom->get_x();
        double selectedatomy = selectedatom->get_y();
        double selectedatomz = selectedatom->get_z();
        double selectedatommass=selectedatom->get_mass();

        inertiamoment+=((selectedatomx-rcmx)*(selectedatomx-rcmx)
                      +(selectedatomy-rcmy)*(selectedatomy-rcmy)
                      +(selectedatomz-rcmz)*(selectedatomz-rcmz))*selectedatommass;
    }

    double gyradious=sqrt(inertiamoment/totalmass);

    return gyradious;

}


double Data_printer::get_grain_alpha(Box * box, Control *control)
{

    double gyradious=get_grain_gyradious(box);
    double alpha= 1.0 - pow(gyradious,3.0)/pow(control->get_gyradiousini(),3.0);
    return alpha;

}

double Data_printer::get_grain_alpha(Box * box, Control *control, double gyradious)
{
    double alpha= 1.0 - pow(gyradious,3.0)/pow(control->get_gyradiousini(),3.0);
    return alpha;
}

long long int Data_printer::set_initial_type_number(Box * box) //se puede paralelizar pero solo se llama 1 vez
{
     long long int sizebox=box->get_box_atoms().size();
    bool invector=false;
    long long int posinvector=0;

    string targettype;

    for ( long long int i=0;i<sizebox;i++)
    {
        targettype=box->get_box_atoms().at(i)->get_atom_type();
        for ( long long int j=0; j<(long long int)typesinbox_initial.size(); j++)
        {
            if (targettype.compare(typesinbox_initial.at(j))==0) {invector=true; posinvector=j; break;}
        }

        if(invector)
        {
            numberinbox_initial.at(posinvector)++;
        }
        else
        {
            typesinbox_initial.push_back(targettype);
            numberinbox_initial.push_back(1);
        }

        invector=false;

    }

    long long int totalnumber=0;
    long long int numbersize = numberinbox_initial.size();
    for (long long int i=0;i<numbersize;i++)
    {
        totalnumber+=numberinbox_initial.at(i);
    }

    return sizebox-totalnumber; //0 if box size and the method ones
}

long long int Data_printer::get_number_inbox_initial(string targettype_) //pasamos de paralelizar, no se llama
{
    long long int size_typesinbox_initial=typesinbox_initial.size();
    for (long long int i=0; i <  size_typesinbox_initial;i++)
    {
        if (targettype_.compare(typesinbox_initial.at(i))==0)
        {
            return numberinbox_initial.at(i);
        }
    }
    return 0;

}

long long int Data_printer::print_mean_dissolved_coordination( Control *control, Tracker * tracker)
{
    string filedataname=filename+".meandiscoord";

    ofstream filemeandiscoord;

    filemeandiscoord.open(filedataname,ios::app);

    vector<Event_definition*> eventdeflist=control->get_list_event_definition();

    long long int sizeeventdeflist= eventdeflist.size();
    long long int sizeeventdef=0;

    if (control->get_step()==0)
    {
        filemeandiscoord<<"time\t\t";

        for (long long int i=0;i<sizeeventdeflist;i++)
        {
            Event_definition * eventdef=eventdeflist.at(i);

            sizeeventdef = eventdef->get_type_neighbours().size();

            for (long long int j=0; j<sizeeventdef; j++)
            {
                filemeandiscoord<<eventdef->get_involved_atom_type()<<"-"<<eventdef->get_type_neighbours_by_pos(j)<<"-"
                <<eventdef->get_distance_neighbours_by_pos(j)<<"\t\t";
            }

        }
        filemeandiscoord<<endl;
    }


    filemeandiscoord<<control->get_time()<<"\t\t";

    for (long long int i=0;i<sizeeventdeflist;i++)
    {
        Event_definition * eventdef=eventdeflist.at(i);

        sizeeventdef = eventdef->get_type_neighbours().size();

        for (long long int j=0; j<sizeeventdef; j++)
        {
            filemeandiscoord<<tracker->get_mean(eventdef->get_involved_atom_type(),
            eventdef->get_type_neighbours_by_pos(j),eventdef->get_distance_neighbours_by_pos(j))<<"\t\t";
        }
    }
    filemeandiscoord<<endl;


    return 0;
}
