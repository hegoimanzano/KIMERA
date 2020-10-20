#include "Data_reader.h"

Data_reader::Data_reader()
{
    //ctor
}

Data_reader::~Data_reader()
{

}

long long int Data_reader::read_create_xyz_file(string * filename_,  Cell * cell)
{

    ifstream xyzfile(filename_->c_str());

    string keyword;
    string value;


    if (xyzfile.is_open())
    {
        while(!xyzfile.eof())
        {
            string number;
            string type_s;
            string x_s;
            string y_s;
            string z_s;

            getline(xyzfile, number);

            long long int numberpos=stoll(number); //this line contains the number of following possitions

            getline(xyzfile, number);  //this line is for comments

            for (long long int i=0;i<numberpos;i++)
            {
                xyzfile >> type_s;
                xyzfile >> x_s;
                xyzfile >> y_s;
                xyzfile >> z_s;

                {
                    cell->add_position(stof(x_s),stof(y_s),stof(z_s), type_s, 1.0);
                }

            }
            break;
        }
    }
    return 0;
}

long long int Data_reader::read_create_xyz_file_without_corners(string * filename_,  Cell * cell, Control * control)
{ //remove readed positions that are in the corners

    ifstream xyzfile(filename_->c_str());

    string keyword;
    string value;


    if (xyzfile.is_open())
    {
        while(!xyzfile.eof())
        {
            string number;
            string type_s;
            string x_s;
            string y_s;
            string z_s;

            getline(xyzfile, number);

            long long int numberpos=stoll(number); //this line contains the number of following possitions

            getline(xyzfile, number);  //this line is for comments

            for (long long int i=0;i<numberpos;i++)
            {
                xyzfile >> type_s;
                xyzfile >> x_s;
                xyzfile >> y_s;
                xyzfile >> z_s;

                double x_value=stof(x_s);
                double y_value=stof(y_s);
                double z_value=stof(z_s);


                //mirar a ver si la posicion esta en las esquinas y en tal caso no annadirla

                if (control->from_xyz_to_u(x_value, y_value, z_value)> 0.9999
                    || control->from_xyz_to_v( y_value, z_value)> 0.9999
                    || control->from_xyz_to_w(z_value)> 0.9999 )
                {
                    //tu gozo en un pozo

                    cout<<"Position at " << x_value<< "  "<< y_value << "  " << z_value << "  discarded because it is in the corner"<<endl;
                }
                else
                {
                    cell->add_position(x_value,y_value,z_value, type_s, 1.0);
                }
            }
            break;
        }
    }
    return 0;
}


long long int Data_reader::read_create_cif_file(string * filename_, Control * control, Cell * cell)  //TODO.. impossible
{
    // check if is truly a cif file
    long long int namesize = filename_->length();
    string extension = filename_->substr(namesize-4,4);

    cout<<extension<<endl;

    if (extension.compare(".cif")!=0  && extension.compare(".CIF")!=0)
    {
        //printf("\nNot valid format,input file has to be a .cif\n");
        cout<<"\nNot valid format,input file has to be a .cif\n "<<endl;
        exit(1);
    }

    ifstream ciffile(filename_->c_str());

    string keyword;
    string value;

    //Starting the reading

    // Previous check -> Find all neccessary keywords
    bool flag_cell_length_a=false;
    bool flag_cell_length_b=false;
    bool flag_cell_length_c=false;
    bool flag_cell_angle_alpha=false;
    bool flag_cell_angle_beta=false;
    bool flag_cell_angle_gamma=false;
    bool flag_space_group_symop_operation_xyz=false;
    bool flag_atom_site_label=false;
    bool flag_atom_site_fract_x=false;
    bool flag_atom_site_fract_y=false;
    bool flag_atom_site_fract_z=false;
    bool flag_atom_site_occupancy=false;


    if (ciffile.is_open())
    {
        while(!ciffile.eof())
        {
            ciffile >> keyword;

            if (keyword.compare("_cell_length_a")==0) flag_cell_length_a=true;
            if (keyword.compare("_cell_length_b")==0) flag_cell_length_b=true;
            if (keyword.compare("_cell_length_c")==0) flag_cell_length_c=true;
            if (keyword.compare("_cell_angle_alpha")==0) flag_cell_angle_alpha=true;
            if (keyword.compare("_cell_angle_beta")==0) flag_cell_angle_beta=true;
            if (keyword.compare("_cell_angle_gamma")==0) flag_cell_angle_gamma=true;
            if (keyword.compare("_space_group_symop_operation_xyz")==0) flag_space_group_symop_operation_xyz=true;
            if (keyword.compare("_atom_site_label")==0) flag_atom_site_label=true;
            if (keyword.compare("_atom_site_fract_x")==0) flag_atom_site_fract_x=true;
            if (keyword.compare("_atom_site_fract_y")==0) flag_atom_site_fract_y=true;
            if (keyword.compare("_atom_site_fract_z")==0) flag_atom_site_fract_z=true;
            if (keyword.compare("_atom_site_occupancy")==0) flag_atom_site_occupancy=true;
        }

        if (!flag_cell_length_a) cout<<"missing _cell_length_a"<<endl;
        if (!flag_cell_length_b) cout<<"missing _cell_length_b"<<endl;
        if (!flag_cell_length_c) cout<<"missing _cell_length_c"<<endl;
        if (!flag_cell_angle_alpha) cout<<"missing _cell_angle_alpha"<<endl;
        if (!flag_cell_angle_beta) cout<<"missing _cell_angle_beta"<<endl;
        if (!flag_cell_angle_gamma) cout<<"missing _cell_angle_gamma"<<endl;
        if (!flag_space_group_symop_operation_xyz) cout<<"missing _space_group_symop_operation_xyz"<<endl;
        if (!flag_atom_site_label) cout<<"missing _atom_site_label"<<endl;
        if (!flag_atom_site_fract_x) cout<<"missing _atom_site_fract_x"<<endl;
        if (!flag_atom_site_fract_y) cout<<"missing _atom_site_fract_y"<<endl;
        if (!flag_atom_site_fract_z) cout<<"missing _atom_site_fract_z"<<endl;

        if( !(flag_cell_length_a && flag_cell_length_b && flag_cell_length_c &&
              flag_cell_angle_alpha && flag_cell_angle_beta && flag_cell_angle_gamma &&
              flag_space_group_symop_operation_xyz &&
              flag_atom_site_label && flag_atom_site_fract_x && flag_atom_site_fract_y && flag_atom_site_fract_z) )
        {
            ciffile.close();
            exit(1);
        }

        cout<<"OK... All required data is in the cif file."<<endl;
    }

    ciffile.clear();
    ciffile.seekg(0);

    vector<Sym_equation*> list_sym_equations;



    if (ciffile.is_open())
    {
        while(!ciffile.eof())
        {
            ciffile >> keyword;

            if (keyword.compare("_cell_length_a")==0)
            {
                ciffile >> value;
                double cell_length_a= stod(value);
                control->set_cell_a(cell_length_a);
            }

            if (keyword.compare("_cell_length_b")==0)
            {
                ciffile >> value;
                double cell_length_b= stod(value);
                control->set_cell_b(cell_length_b);
            }

            if (keyword.compare("_cell_length_c")==0)
            {
                ciffile >> value;
                double cell_length_c= stod(value);
                control->set_cell_c(cell_length_c);
            }

            if (keyword.compare("_cell_angle_alpha")==0)
            {
                ciffile >> value;
                double _cell_angle_alpha= stod(value);
                control->set_angle_alpha(_cell_angle_alpha);
            }
            if (keyword.compare("_cell_angle_beta")==0)
            {
                ciffile >> value;
                double _cell_angle_beta= stod(value);
                control->set_angle_beta(_cell_angle_beta);
            }
            if (keyword.compare("_cell_angle_gamma")==0)
            {
                ciffile >> value;
                double _cell_angle_gamma= stod(value);
                control->set_angle_gamma(_cell_angle_gamma);
            }

            if (keyword.compare("_space_group_symop_operation_xyz")==0)    //CAMBIAR ESTO?
            {
                //crear una instancia de sym equation por cada ecuacion
                ciffile >> value;
                while (value.compare("loop_")!=0)
                {
                    Sym_equation * sym_equation=new Sym_equation(value);
                    list_sym_equations.push_back(sym_equation);
                    ciffile >> value;
                }
            }

            if (keyword.compare("_atom_site_fract_z")==0)
            {
                long long int extra=0;
                ciffile >> value;

                //avanzamos una posicion mas si esta la ocupacion
                if (flag_atom_site_occupancy)
                {
                    ciffile >> value;
                }

                //por si hay mas cosas de las que queremos
                while (value.compare(0,1,"_")==0)
                {
                    extra++;
                    ciffile >> value;
                }

                string type;
                double x_;
                double y_;
                double z_;
                double occupancy;
                long long int id_ = control->get_totalnumberpositions();

                // mientras no llegue al final, o la segunda palabra que pille no sea un numero..



                //leer posiciones mientras que la palabra siguiente sea distinta a loop o que no sea el final del documento
                while (value.compare("loop_")!=0 && !ciffile.eof())
                {

                    type=value;

                    cout<<value.compare("loop_")<<endl;
                    cout<<value<<endl;

                    ciffile >> value;
                    x_=stod(value);

                    ciffile >> value;
                    y_=stod(value);

                    ciffile >> value;
                    z_=stod(value);

                    if (flag_atom_site_occupancy)
                    {
                        ciffile >> value;
                        occupancy=stod(value);
                    }
                    else
                    {
                        occupancy=1.0;
                    }

                    //Pasamos de los extras

                    for (long long int i=0;i<extra;i++)
                    {
                        ciffile >> value;
                    }


                    //Miramos que la posición que vamos a meter no exista ya

                    vector<Position*> list_positions = cell->get_cell_positions();
                    Position positionl = Position(x_,y_,z_,id_, type, occupancy);

                    bool exist=false;
                    for ( long long int i = 0; i<(long long int)list_positions.size();i++)
                    {
                        if (positionl.compare_position(list_positions.at(i)))
                        {
                            exist=true;

                            //metemos a la posicion que existe el type_prob
                            list_positions.at(i)->add_type_prob(type, occupancy);
                            break;
                        }

                    }

                    //Si no existe la metemos
                    if (!exist) cell->add_position(x_,y_,z_, type, occupancy);


                    ciffile >> value;
                }

            }

        }
        ciffile.close();
    }
//METER TMB LAS SIMETRICAS







    return 0; //Everything OK
}


long long int Data_reader::read_input_file_dis_events(string * filename_, Control * control)
{
    ifstream inputfile(filename_->c_str());
    bool readed=false;

    string keyword;
    string value;

    string value2;
    string value3;
    string value4;
    string value5;
    string value6;
    string value7;
    bool any=false;

    bool inffd=false;
    bool inffp=false;
    bool indeltaG=false;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if(!readed) inputfile >> keyword;
            readed=false;

            if (keyword.compare("DEFINE_DISSOLUTION_EVENT")==0)
            {
                string type;
                inputfile >> type;

                Event_definition * event_definition1 = new Event_definition();
                event_definition1->set_involved_atom_type(type);

                //Nuevo
                string list_length;
                string linkedtype;
                string neighbourtype;
                double neighbourdistance;
                double neighbourEd;
                double neighbourEp;
                long long int max_neighbour_to_bulk;
                //Nuevo
                long long int longitud_lista;

                string affectedtype;
                double affecteddistance;

                double ffd;
                double ffp;
                double deltaG;


                while(!inputfile.eof())
                {
                    if(!readed) inputfile >> keyword;
                    readed=false;
                    if (keyword.compare("NEIGHBOUR")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();
                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;
                        linked_neighbour1.set_type(NO_LINKED_TYPE);  //esto quiere decir que no es linked
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        inputfile>>neighbourtype;
                        inputfile>>value;
                        inputfile>>value2;
                        inputfile>>value3;
                        inputfile>>value4;

                        try{
                            neighbourdistance=stod(value);
                            neighbourEd=stod(value2);
                            neighbourEp=stod(value3);
                            max_neighbour_to_bulk=stoll(value4);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }
                        if (neighbourdistance<0.0 || neighbourEd<0.0 || neighbourEp<0.0 || max_neighbour_to_bulk<0)
                        {
                            cout<<endl<<"Not valid number in NEIGHBOUR"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        cout <<endl<< "Added NEIGHBOUR of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED " <<   neighbourEd  << " KT units and precipitation energy EP  " <<   neighbourEp
                        << " KT units.  The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;


                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }
                    //Nuevo
                    if (keyword.compare("NEIGHBOUR_DIRECT_LIST")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();

                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;
                        double valorEd;
                        double valorEp;

                        linked_neighbour1.set_type(NO_LINKED_TYPE);  //esto quiere decir que no es linked
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        inputfile>>neighbourtype;
                        inputfile>>value;
                        inputfile>>list_length;

                        //comprobamos que es el comando LIST_LENGTH

                        if (list_length.compare("LIST_LENGTH")!=0)
                        {
                            cout << "Invalid argument in NEIGHBOUR_DIRECT_LIST: LIST_LENGTH is missing " << endl;
                            inputfile.close();
                            exit(1);
                        }

                        inputfile>>value5;

                        try{
                            neighbourdistance=stod(value);
                            longitud_lista=stoll(value5);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        //comprobar que la longitud_lista es entero positiva

                        if (neighbourdistance<0.0 || longitud_lista<0 )
                        {
                            cout<<endl<<"Not valid number in NEIGHBOUR_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        //los n siguientes numeros, son los valores de Ed, y los n siguientes los de Ep

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            inputfile>>value6;
                            inputfile>>value7;

                            try{
                                valorEd =stod(value6);
                                valorEp =stod(value7);
                            }

                            catch (const invalid_argument& ia) {
                                cout << "Invalid argument in the list in NEIGHBOUR_DIRECT_LIST: " << ia.what() << endl;
                                inputfile.close();
                                exit(1);
                            }

                            Ed_direct.push_back(valorEd);
                            Ep_direct.push_back(valorEp);

                        }

                        //ahora me falta el max_to_bulk


                        inputfile>>value4;

                        try{
                            max_neighbour_to_bulk=stoll(value4);

                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in neighbours_to_bulk in NEIGHBOUR_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        if (max_neighbour_to_bulk<0)
                        {
                            cout<<endl<<"Not valid number in neighbours_to_bulk in NEIGHBOUR_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);


                        cout <<endl<< "Added NEIGHBOUR_DIRECT_LIST of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED and precipitation energy EP "<<endl;

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            cout<<"Number of Neighbours: "<< i+1 << "  ED:  " << Ed_direct.at(i) << "  EP:  "<< Ep_direct.at(i) << "     KT units   " <<endl;
                        }

                        cout<< " The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;

                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }




                    if (keyword.compare("NEIGHBOUR_LINKED_DIRECT_LIST")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();
                        linked_neighbour1.set_type(LINKED_TYPE_NORMAL);  //esto quiere decir que el linkeado es si el atomo es normal.
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;
                        double valorEd;
                        double valorEp;


                        inputfile>>neighbourtype;
                        inputfile>>value;

                        try{
                            neighbourdistance=stod(value);
                        }
                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        inputfile>>value;

                        while(!inputfile.eof())
                        {
                            if (value.compare("LINK")==0)
                            {
                                double distance_ori_lnk=0.0;
                                double distance_target_lnk=0.0;

                                inputfile>>linkedtype;
                                inputfile>>value2;
                                inputfile>>value3;

                                try{
                                    distance_ori_lnk=stod(value2);
                                    distance_target_lnk=stod(value3);
                                }
                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in LINK in NEIGHBOUR_LINKED_DIRECT_LIST: " << ia.what() << endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                linked_neighbour1.add_linked_neighbour(linkedtype,distance_ori_lnk,distance_target_lnk);

                                inputfile>>value;
                            }
                            else{break;}
                        }


                        //comprobamos que es el comando LIST_LENGTH

                        if (value.compare("LIST_LENGTH")!=0)
                        {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DIRECT_LIST: LIST_LENGTH is missing " << endl;
                            inputfile.close();
                            exit(1);
                        }

                        inputfile>>value5;

                        try{
                            longitud_lista=stoll(value5);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        //comprobar que la longitud_lista es entero positiva

                        if (longitud_lista<0 )
                        {
                            cout<<endl<<"Not valid number in NEIGHBOUR_LINKED_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        //los n siguientes numeros, son los valores de Ed, y los n siguientes los de Ep

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            inputfile>>value6;
                            inputfile>>value7;

                            try{
                                valorEd =stod(value6);
                                valorEp =stod(value7);
                            }

                            catch (const invalid_argument& ia) {
                                cout << "Invalid argument in the list in NEIGHBOUR_LINKED_DIRECT_LIST: " << ia.what() << endl;
                                inputfile.close();
                                exit(1);
                            }

                            Ed_direct.push_back(valorEd);
                            Ep_direct.push_back(valorEp);

                        }


                        //ahora va el max to bulk
                        inputfile>>value4;

                        try{
                            max_neighbour_to_bulk=stoll(value4);

                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in neighbours_to_bulk in NEIGHBOUR_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        if (max_neighbour_to_bulk<0)
                        {
                            cout<<endl<<"Not valid number in neighbours_to_bulk in NEIGHBOUR_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);



                        cout <<endl<< "Added NEIGHBOUR_LINKED_DIRECT_LIST of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED and precipitation energy EP "<<endl;

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            cout<<"Number of Neighbours: "<< i+1 << "  ED:  " << Ed_direct.at(i) << "  EP:  "<< Ep_direct.at(i) << "     KT units   " <<endl;
                        }

                        cout<< " The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;

                        cout << "Following linked atoms must remain not dissolved: "<<endl;
                        long long int linkedsize=linked_neighbour1.get_distance_origin_linked().size();

                        for (long long int i=0;i<linkedsize;i++)
                        {
                            cout<< linked_neighbour1.get_type_linked().at(i) << " at a distance from the origin atom of "
                            << linked_neighbour1.get_distance_origin_linked().at(i) << " and at a distance from the neighbour atom of "
                            << linked_neighbour1.get_distance_target_linked().at(i) << endl;
                        }

                        control->set_linked(true);  //creo que es un vestigio
                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }


                    if (keyword.compare("NEIGHBOUR_LINKED")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();
                        linked_neighbour1.set_type(LINKED_TYPE_NORMAL);  //esto quiere decir que el linkeado es si el atomo es normal.
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;

                        inputfile>>neighbourtype;
                        inputfile>>value;

                        try{
                            neighbourdistance=stod(value);
                        }
                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        inputfile>>value;

                        while(!inputfile.eof())
                        {
                            if (value.compare("LINK")==0)
                            {
                                double distance_ori_lnk=0.0;
                                double distance_target_lnk=0.0;

                                inputfile>>linkedtype;
                                inputfile>>value2;
                                inputfile>>value3;

                                try{
                                    distance_ori_lnk=stod(value2);
                                    distance_target_lnk=stod(value3);
                                }
                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in LINK in NEIGHBOUR_LINKED: " << ia.what() << endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                linked_neighbour1.add_linked_neighbour(linkedtype,distance_ori_lnk,distance_target_lnk);

                                inputfile>>value;
                            }
                            else{break;}
                        }

                        //value tiene el valor de Ed
                        inputfile>>value2;
                        inputfile>>value3;

                        try{
                            neighbourEd=stod(value);
                            neighbourEp=stod(value2);
                            max_neighbour_to_bulk=stoll(value3);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in in NEIGHBOUR_LINKED: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);


                        cout <<endl<< "Added LINKED_NEIGHBOUR of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED " <<   neighbourEd  << " KT units and precipitation energy EP  " <<   neighbourEp
                        << " KT units.  The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;

                        cout << "Following linked atoms must remain not dissolved: "<<endl;
                        long long int linkedsize=linked_neighbour1.get_distance_origin_linked().size();

                        for (long long int i=0;i<linkedsize;i++)
                        {
                            cout<< linked_neighbour1.get_type_linked().at(i) << " at a distance from the origin atom of "
                            << linked_neighbour1.get_distance_origin_linked().at(i) << " and at a distance from the neighbour atom of "
                            << linked_neighbour1.get_distance_target_linked().at(i) << endl;
                        }

                        control->set_linked(true);
                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }


                    if (keyword.compare("NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();
                        linked_neighbour1.set_type(LINKED_TYPE_DISSOLVED);  //esto quiere decir que el linkeado es si el atomo es normal.
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;
                        double valorEd;
                        double valorEp;


                        inputfile>>neighbourtype;
                        inputfile>>value;

                        try{
                            neighbourdistance=stod(value);
                        }
                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        inputfile>>value;

                        while(!inputfile.eof())
                        {
                            if (value.compare("LINK")==0)
                            {
                                double distance_ori_lnk=0.0;
                                double distance_target_lnk=0.0;

                                inputfile>>linkedtype;
                                inputfile>>value2;
                                inputfile>>value3;

                                try{
                                    distance_ori_lnk=stod(value2);
                                    distance_target_lnk=stod(value3);
                                }
                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in LINK in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST: " << ia.what() << endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                linked_neighbour1.add_linked_neighbour(linkedtype,distance_ori_lnk,distance_target_lnk);

                                inputfile>>value;
                            }
                            else{break;}
                        }


                        //comprobamos que es el comando LIST_LENGTH

                        if (value.compare("LIST_LENGTH")!=0)
                        {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DIRECT_LIST: LIST_LENGTH is missing " << endl;
                            inputfile.close();
                            exit(1);
                        }

                        inputfile>>value5;

                        try{
                            longitud_lista=stoll(value5);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        //comprobar que la longitud_lista es entero positiva

                        if (longitud_lista<0 )
                        {
                            cout<<endl<<"Not valid number in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        //los n siguientes numeros, son los valores de Ed, y los n siguientes los de Ep

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            inputfile>>value6;
                            inputfile>>value7;

                            try{
                                valorEd =stod(value6);
                                valorEp =stod(value7);
                            }

                            catch (const invalid_argument& ia) {
                                cout << "Invalid argument in the list in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST: " << ia.what() << endl;
                                inputfile.close();
                                exit(1);
                            }

                            Ed_direct.push_back(valorEd);
                            Ep_direct.push_back(valorEp);

                        }

                        //ahora va el max to bulk
                        inputfile>>value4;

                        try{
                            max_neighbour_to_bulk=stoll(value4);

                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in neighbours_to_bulk in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }


                        if (max_neighbour_to_bulk<0)
                        {
                            cout<<endl<<"Not valid number in neighbours_to_bulk in NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST"<<endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);


                        cout <<endl<< "Added NEIGHBOUR_LINKED_DISSOLVED_DIRECT_LIST of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED and precipitation energy EP "<<endl;

                        for (long long int i=0;i<longitud_lista;i++)
                        {
                            cout<<"Number of Neighbours: "<< i+1 << "  ED:  " << Ed_direct.at(i) << "  EP:  "<< Ep_direct.at(i) << "     KT units   " <<endl;
                        }

                        cout<< " The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;

                        cout << "Following linked atoms must remain dissolved: "<<endl;
                        long long int linkedsize=linked_neighbour1.get_distance_origin_linked().size();

                        for (long long int i=0;i<linkedsize;i++)
                        {
                            cout<< linked_neighbour1.get_type_linked().at(i) << " at a distance from the origin atom of "
                            << linked_neighbour1.get_distance_origin_linked().at(i) << " and at a distance from the neighbour atom of "
                            << linked_neighbour1.get_distance_target_linked().at(i) << endl;
                        }

                        control->set_linked(true);  //creo que es un vestigio
                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }


                    if (keyword.compare("NEIGHBOUR_LINKED_DISSOLVED")==0)
                    {
                        Linked_neighbour linked_neighbour1=Linked_neighbour();
                        linked_neighbour1.set_type(LINKED_TYPE_DISSOLVED);  //esto quiere decir que el linkeado es si el atomo es disuelto.
                        control ->set_id_linked(control->get_id_linked()+1);
                        linked_neighbour1.set_id(control->get_id_linked());

                        //Nuevo
                        vector<double> Ed_direct;
                        vector<double> Ep_direct;

                        inputfile>>neighbourtype;
                        inputfile>>value;

                        try{
                            neighbourdistance=stod(value);
                        }
                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOUR_LINKED_DISSOLVED: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        inputfile>>value;

                        while(!inputfile.eof())
                        {
                            if (value.compare("LINK")==0)
                            {
                                double distance_ori_lnk=0.0;
                                double distance_target_lnk=0.0;

                                inputfile>>linkedtype;
                                inputfile>>value2;
                                inputfile>>value3;

                                try{
                                    distance_ori_lnk=stod(value2);
                                    distance_target_lnk=stod(value3);
                                }
                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in LINK in NEIGHBOUR_LINKED_DISSOLVED: " << ia.what() << endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                linked_neighbour1.add_linked_neighbour(linkedtype,distance_ori_lnk,distance_target_lnk);

                                inputfile>>value;

                            }
                            else{break;}
                        }

                        //value tiene el valor de Ed
                        inputfile>>value2;
                        inputfile>>value3;

                        try{
                            neighbourEd=stod(value);
                            neighbourEp=stod(value2);
                            max_neighbour_to_bulk=stoll(value3);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in in NEIGHBOUR_LINKED_DISSOLVED: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_neighbours(neighbourtype);
                        event_definition1->add_distance_neighbours(neighbourdistance);
                        event_definition1->add_Ed_neighbours(neighbourEd);
                        event_definition1->add_Ep_neighbours(neighbourEp);
                        event_definition1->add_max_to_bulk(max_neighbour_to_bulk);

                        //Nuevo
                        event_definition1->add_Ed_neighbours_direct(Ed_direct);
                        event_definition1->add_Ep_neighbours_direct(Ep_direct);


                        event_definition1->add_to_list_linked_neighbour(linked_neighbour1);


                        cout <<endl<< "Added LINKED_NEIGHBOUR_DISSOLVED of  "<< neighbourtype <<"  to an atom of  " << type << "  in distance   " <<  neighbourdistance
                        << " with dissolution energy ED " <<   neighbourEd  << " KT units and precipitation energy EP  " <<   neighbourEp
                        << " KT units.  The amount of neighbours which define a bulk atom are "   <<max_neighbour_to_bulk<< " (0->not define)" << endl;

                        cout << "Following linked atoms must remain dissolved: "<<endl;
                        long long int linkedsize=linked_neighbour1.get_distance_origin_linked().size();

                        for (long long int i=0;i<linkedsize;i++)
                        {
                            cout<< linked_neighbour1.get_type_linked().at(i) << " at a distance from the origin atom of "
                            << linked_neighbour1.get_distance_origin_linked().at(i) << " and at a distance from the neighbour atom of "
                            << linked_neighbour1.get_distance_target_linked().at(i) << endl;
                        }

                        control->set_linked(true);
                        any=true;
                        inputfile >> keyword;
                        readed=true;
                    }


                    if (keyword.compare("AFFECTED")==0)
                    {
                        inputfile>>affectedtype;
                        inputfile>>value;

                        try{

                            affecteddistance=stod(value);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in AFFECTED: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }

                        event_definition1->add_type_affected(affectedtype);
                        event_definition1->add_distance_affected(affecteddistance);

                        cout <<endl<< "Added AFFECTED of  "<< affectedtype <<"  to an atom of  " << type  << " in distance "<< affecteddistance << endl;

                        control->set_affected(true);
                        inputfile >> keyword;
                        readed=true;

                    }

                    if (keyword.compare("FFD")==0) // initialized latter if there is no FFD
                    {
                        inputfile>>value;

                        try{
                            ffd=stod(value);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in FFD: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }
                        inffd=true;
                        event_definition1->set_ffd(ffd);

                        cout <<endl<< "DISSOLUTION_FUNDAMENTAL_FREQUENCY for this event:   "<< ffd << endl;

                        inputfile >> keyword;
                        readed=true;

                    }

                    if (keyword.compare("FFP")==0) // initialized latter if there is no FFP
                    {
                        inputfile>>value;

                        try{
                            ffp=stod(value);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in FFP: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }
                        inffp=true;
                        event_definition1->set_ffp(ffp);
                        cout <<endl<< "PRECIPITATION_FUNDAMENTAL_FREQUENCY for this event:   "<< ffp << endl;
                        inputfile >> keyword;
                        readed=true;

                    }

                    if (keyword.compare("DG*")==0) // initialized latter if there is no DG*
                    {
                        inputfile>>value;

                        try{
                            deltaG=stod(value);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in DG*: " << ia.what() << endl;
                            inputfile.close();
                            exit(1);
                        }
                        indeltaG=true;
                        event_definition1->set_deltaG(deltaG);
                        cout <<endl<< "DELTA_G* for this event:   "<< deltaG << endl;
                        inputfile >> keyword;
                        readed=true;

                    }

                    if (keyword.compare("DEFINE_DISSOLUTION_EVENT")==0){break;}
                }

            if (!inffd) event_definition1->set_ffd(control->get_fd());
            if (!inffp) event_definition1->set_ffp(control->get_fr());
            if (!indeltaG) event_definition1->set_deltaG(control->get_deltag());
            control->add_list_event_definition(event_definition1);
            inffd=false;
            inffp=false;
            indeltaG=false;
            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    if (!any) cout<< "You should use at least one DEFINE_DISSOLUTION_EVENT" <<endl;
    return 0;


}

long long int Data_reader::read_input_file(string * filename_, Control * control)
{
    bool ffdbool=false;
    bool ffrbool=false;
    bool targettimebool=false;
    bool stimatetimebool=false;
    bool targetstepbool=false;
    bool parallelize=false;
    bool seedsim=false;
    bool readed=false;
    bool boolxyzread=false;
    bool boolkimeraread=false;

    ifstream inputfile(filename_->c_str());

    string keyword;
    string value;

    string value2;
    string value3;
    string value4;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if(!readed) inputfile >> keyword;

            readed=false;


            if (keyword.compare("DIMENSION_A")==0)
            {
                long long int dimension_x;
                inputfile >> value;


                try {
                    dimension_x=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DIMENSION_A: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dimension_x<1)
                {
                    cout<<endl<<"Not valid number in DIMENSION_A"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_dimension_x(dimension_x);
                 control->set_tilt_bounds_factors();
                 cout<< "DIMENSION_A:  " << dimension_x<<endl;

                 continue;
            }

            if (keyword.compare("DIMENSION_B")==0)
            {
                long long int dimension_y;
                inputfile >> value;


                try {
                    dimension_y=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DIMENSION_B: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dimension_y<1)
                {
                    cout<<endl<<"Not valid number in DIMENSION_B"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_dimension_y(dimension_y);
                 control->set_tilt_bounds_factors();
                 cout<< "DIMENSION_B:  " << dimension_y<<endl;

                 continue;
            }

            if (keyword.compare("DIMENSION_C")==0)
            {
                long long int dimension_z;
                inputfile >> value;

                try {
                    dimension_z=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DIMENSION_C: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dimension_z<1)
                {
                    cout<<endl<<"Not valid number in DIMENSION_C"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_dimension_z(dimension_z);
                 control->set_tilt_bounds_factors();
                 cout<< "DIMENSION_C:  " << dimension_z<<endl;

                 continue;
            }

            if (keyword.compare("CELL_A")==0)
            {
                double cell_a;
                inputfile >> value;

                try {
                    cell_a=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_A: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (cell_a<0)
                {
                    cout<<endl<<"Not valid number in CELL_A"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_cell_a(cell_a);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_A:  " << cell_a<<endl;

                 continue;
            }

            if (keyword.compare("CELL_B")==0)
            {
                double cell_b;
                inputfile >> value;

                try {
                    cell_b=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_B: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (cell_b<0)
                {
                    cout<<endl<<"Not valid number in CELL_B"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_cell_b(cell_b);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_B:  " << cell_b<<endl;

                 continue;
            }

            if (keyword.compare("CELL_C")==0)
            {
                double cell_c;
                inputfile >> value;

                try {
                    cell_c=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_C: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (cell_c<0)
                {
                    cout<<endl<<"Not valid number in CELL_C"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_cell_c(cell_c);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_C:  " << cell_c<<endl;

                 continue;
            }

            if (keyword.compare("CELL_ALPHA")==0)
            {
                double cell_alpha;
                inputfile >> value;

                try {
                    cell_alpha=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_ALPHA: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_angle_alpha(cell_alpha);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_ALPHA:  " << cell_alpha<<endl;

                 continue;
            }

            if (keyword.compare("CELL_BETA")==0)
            {
                double cell_beta;
                inputfile >> value;

                try {
                    cell_beta=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_BETA: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_angle_beta(cell_beta);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_BETA:  " << cell_beta<<endl;

                 continue;
            }

            if (keyword.compare("CELL_GAMMA")==0)
            {
                double cell_gamma;
                inputfile >> value;

                try {
                    cell_gamma=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_GAMMA: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_angle_gamma(cell_gamma);
                 control->set_tilt_bounds_factors();
                 cout<< "CELL_GAMMA:  " << cell_gamma <<endl;

                 continue;
            }

            if (keyword.compare("TARGET_TIME")==0)
            {

                inputfile >> value;
                double dis_time;

                try {
                    dis_time=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in TARGET_TIME: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dis_time<0.0)
                {
                    cout<<endl<<"Not valid number in TARGET_TIME"<<endl;
                    inputfile.close();
                    exit(1);
                }
                 targettimebool=true;
                 control->set_targettime(dis_time);
                 cout<< "TARGET_TIME:  " << dis_time <<endl;

                 continue;
            }

            if (keyword.compare("TARGET_STEP")==0)
            {

                inputfile >> value;
                double target_step;

                try {
                    target_step=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in TARGET_STEP: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (target_step<1)
                {
                    cout<<endl<<"Not valid number in TARGET_STEP"<<endl;
                    inputfile.close();
                    exit(1);
                }
                 targetstepbool=true;
                 control->set_targetsteps(target_step);
                 control->set_print_by_steps(true);
                 cout<< "TARGET_STEP:  " << target_step <<endl;

                 continue;
            }

            if (keyword.compare("ESTIMATE_TIME")==0)
            {
                 stimatetimebool=true;
                 control->set_estimatetime(true);
                 cout<< "ESTIMATE_TIME: true "<<endl;

                 continue;
            }

            if (keyword.compare("PARALLELIZE_SIMULATION")==0)
            {
                inputfile >> value;
                long long int numbercores;

                try {
                    numbercores=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in PARALLELIZE_SIMULATION: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }
                if (numbercores<1)
                {
                    cout << "The number of cores in PARALLELIZE_SIMULATION must greater than 1: " << endl;
                    cout<< "PARALLELIZE_SIMULATION:  " <<  numbercores <<endl;
                    exit(1);
                }


                parallelize=true;
                if (numbercores==1){parallelize=false;}

                control->set_parallelize(true);
                control->set_parallelizenumber(numbercores);

                cout<< "PARALLELIZE_SIMULATION: "<< numbercores <<endl;

                continue;
            }

            /**
            if (keyword.compare("TEMPERATURE")==0)
            {

                inputfile >> value;
                double temperature;

                try {
                    temperature=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in TEMPERATURE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (temperature<0.0000001)
                {
                    cout<<endl<<"Not valid number in TEMPERATURE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_temperature(temperature);
                 cout<< "TEMPERATURE:  " << temperature <<endl;

                 continue;
            }
            */

            if (keyword.compare("DELTA_G*")==0)
            {
                double deltaG;
                inputfile >> value;

                try {
                    deltaG=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DELTA_G*: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                 control->set_deltag(deltaG);
                 cout<< "DELTA_G*:  " << deltaG <<endl;

                 continue;
            }

            if (keyword.compare("READ_POSITIONS_FROM_XYZ_FILE")==0)
            {
                inputfile >> value;

                boolxyzread=true;
                control->set_xyzfile(true);
                control->set_pathtoxyzfile(value);
                cout<< "READ_POSITIONS_FROM_XYZ_FILE:  " << control->get_pathtoxyzfile() <<endl;

                continue;
            }

            if (keyword.compare("READ_SYSTEM_FROM_KIMERA_FILE")==0)
            {
                inputfile >> value;

                boolkimeraread=true;
                control->set_kimerafile(true);
                control->set_pathtokimerafile(value);
                cout<< "READ_SYSTEM_FROM_KIMERA_FILE:  " << control->get_pathtokimerafile() <<endl;

                continue;
            }


            if (keyword.compare("SEED_BOX")==0)
            {
                long long int seedbox;
                inputfile >> value;

                try {
                    seedbox=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in SEED_BOX: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }
                 control->set_seedboxbool(true);
                 control->set_seedbox(seedbox);
                 cout<< "SEED_BOX:  " << seedbox <<endl;

                 continue;
            }

            if (keyword.compare("LINEAL_SEARCH")==0)
            {
                 control->set_linealsearch(true);
                 cout<< "LINEAL_SEARCH: true " <<endl;

                 continue;
            }

            if (keyword.compare("INITIAL_KIMERA_STATE")==0)
            {
                 control->set_printinitialkimerafile(true);
                 cout<< "INITIAL_KIMERA_STATE: true " <<endl;

                 continue;
            }

            if (keyword.compare("FINAL_KIMERA_STATE")==0)
            {
                 control->set_printfinalkimerafile(true);
                 cout<< "FINAL_KIMERA_STATE: true " <<endl;

                 continue;
            }

            if (keyword.compare("SEED_SIMULATION")==0)
            {
                long long int seedsimnumber;
                inputfile >> value;

                try {
                    seedsimnumber=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in SEED_SIMULATION: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }
                 seedsim=true;
                 control->set_seedsimbool(true);
                 control->set_seedsim(seedsimnumber);
                 cout<< "SEED_SIMULATION:  " << seedsimnumber <<endl;

                 continue;
            }

            if (keyword.compare("DISSOLUTION_FUNDAMENTAL_FREQUENCY")==0)
            {
                double ffd;
                inputfile >> value;

                try {
                    ffd=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DISSOLUTION_FUNDAMENTAL_FREQUENCY: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (ffd<0.0)
                {
                    cout<<endl<<"Not valid number in DISSOLUTION_FUNDAMENTAL_FREQUENCY"<<endl;
                    inputfile.close();
                    exit(1);
                }
                ffdbool=true;
                 control->set_fd(ffd);
                 cout<< "DISSOLUTION_FUNDAMENTAL_FREQUENCY:  " << ffd <<endl;

                 continue;
            }

            if (keyword.compare("PRECIPITATION_FUNDAMENTAL_FREQUENCY")==0)
            {
                inputfile >> value;
                double ffr;

                try {
                    ffr=stod(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in PRECIPITATION_FUNDAMENTAL_FREQUENCY: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (ffr<0.0)
                {
                    cout<<endl<<"Not valid number in PRECIPITATION_FUNDAMENTAL_FREQUENCY"<<endl;
                    inputfile.close();
                    exit(1);
                }
                 ffrbool=true;
                 control->set_fr(ffr);
                 cout<< "PRECIPITATION_FUNDAMENTAL_FREQUENCY:  " << ffr <<endl;

                 continue;
            }

            if (keyword.compare("ACCURACY")==0)
            {
                inputfile >> value;
                double accuracy;

                try {
                    accuracy=stof(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ACCURACY: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (accuracy<0.0)
                {
                    cout<<endl<<"Not valid number in ACCURACY"<<endl;
                    inputfile.close();
                    exit(1);
                }
                 control->set_unstuckaccuracy(accuracy);
                 cout<< "ACCURACY:  " << accuracy <<endl;

                 continue;
            }


            if (keyword.compare("DISTANCE_ACCURACY")==0)
            {
                inputfile >> value;
                double distance_accuracy;

                try {
                    distance_accuracy=stof(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DISTANCE_ACCURACY: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance_accuracy<0.0)
                {
                    cout<<endl<<"Not valid number in DISTANCE_ACCURACY"<<endl;
                    inputfile.close();
                    exit(1);
                }
                 double previous_distance_accuracy=control->get_distance_accuracy();
                 control->change_distance_accuracy(distance_accuracy);
                 cout<< "Watch out!! you have changed DISTANCE_ACCURACY, which may have a great impact in your definition of events.  "
                 <<endl << " Previous DISTANCE_ACCURACY of  "<< previous_distance_accuracy<<" has been changed to  "<<
                 control->get_distance_accuracy() <<endl;

                 continue;
            }

            if (keyword.compare("MEAN_DISSOLVED_ANALYSIS")==0)
            {
                inputfile >> value;
                long long int mean_frames;

                try {
                    mean_frames=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in MEAN_DISSOLVED_ANALYSIS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (mean_frames<1)
                {
                    cout<<endl<<"Not valid number in MEAN_DISSOLVED_ANALYSIS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                control->set_meandiscoord(true);
                control->set_meansteps(mean_frames);
                cout<< "MEAN_DISSOLVED_ANALYSIS:  " << mean_frames <<endl;

                continue;
            }


            if (keyword.compare("DATA_ANALYSIS")==0)
            {
                inputfile >> value;
                long long int data_frames;

                try {
                    data_frames=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DATA_ANALYSIS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (data_frames<1)
                {
                    cout<<endl<<"Not valid number in DATA_ANALYSIS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                control->set_dataanalisys(true);
                control->set_datastepts(data_frames);
                cout<< "DATA_ANALYSIS:  " << data_frames <<endl;

                continue;
            }

            if (keyword.compare("BOX_FRAMES")==0)
            {
                inputfile >> value;
                long long int box_frames;

                try {
                    box_frames=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in BOX_FRAMES: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (box_frames<1)
                {
                    cout<<endl<<"Not valid number in BOX_FRAMES"<<endl;
                    inputfile.close();
                    exit(1);
                }
                control->set_boxframe(true);
                control->set_boxsteps(box_frames);
                cout<< "BOX_FRAMES:  " << box_frames <<endl;

                continue;
            }

            if (keyword.compare("SURFACE_FRAMES")==0)
            {
                inputfile >> value;
                long long int surface_frames;

                try {
                    surface_frames=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in SURFACE_FRAMES: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (surface_frames<1)
                {
                    cout<<endl<<"Not valid number in SURFACE_FRAMES"<<endl;
                    inputfile.close();
                    exit(1);
                }
                control->set_surfaceframe(true);
                control->set_surfacesteps(surface_frames);
                cout<< "SURFACE_FRAMES:  " << surface_frames <<endl;

                continue;
            }

            if (keyword.compare("LAYER_ANALYSIS")==0)
            {
                long long int layer_frames;

                inputfile >> value;

                if (value.compare("A")==0 )
                {
                    control->set_layerxanalysis(true);
                    inputfile >> value;
                }
                if (value.compare("B")==0)
                {
                    control->set_layeryanalysis(true);
                    inputfile >> value;
                }
                if (value.compare("C")==0)
                {
                    control->set_layerzanalysis(true);
                    inputfile >> value;
                }


                try {
                    layer_frames=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in LAYER_ANALYSIS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (layer_frames<1)
                {
                    cout<<endl<<"Not valid number in LAYER_ANALYSIS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                control->set_layersteps(layer_frames);
                cout<< "LAYER_ANALYSIS:  " << layer_frames <<endl;

                continue;
            }

            if (keyword.compare("WORK_NAME")==0)  //control->set_plane_x()
            {
                inputfile >> value;

                control->set_workname(value);
                cout<< "WORK_NAME: " << value <<endl;

                continue;
            }

            if (keyword.compare("PERIODICITY")==0)  //control->set_plane_x()
            {
                inputfile >> keyword;

                if (keyword.compare("A")==0 )
                {
                    control->set_plane_x(true);
                    control->set_grain(false);
                    inputfile >> keyword;
                    readed=true;
                }
                if (keyword.compare("B")==0)
                {
                    control->set_plane_y(true);
                    control->set_grain(false);
                    inputfile >> keyword;
                    readed=true;
                }
                if (keyword.compare("C")==0)
                {
                    control->set_plane_z(true);
                    control->set_grain(false);
                    readed=false;
                }
                string periodicitya="false";
                string periodicityb="false";
                string periodicityc="false";

                if (control->get_plane_x()) {periodicitya="true";}
                if (control->get_plane_y()) {periodicityb="true";}
                if (control->get_plane_z()) {periodicityc="true";}

                cout<< "PERIODICITY IN A: " << periodicitya << " B : "<< periodicityb << " C : "<< periodicityc <<endl;

                continue;
            }

        }
        inputfile.close();
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    if ((targettimebool && stimatetimebool) || (targettimebool && targetstepbool) || (stimatetimebool && targetstepbool)) {cout<< "Wacth out!! ESTIMATE_TIME, TARGET_TIME and TARGET_STEP are incompatible" << endl; exit(1);}
    if (!targettimebool && !stimatetimebool && !targetstepbool)
    {
        cout<< "Wacth out!! Neither ESTIMATE_TIME, TARGET_TIME or TARGET_STEP is used,  Considered TARGET_STEP is: "<< control->get_targetsteps() << endl;
        control->set_print_by_steps(true);
    }
    if (seedsim && parallelize) {cout<< "Wacth out!! PARALLELIZE_SIMULATION and SEED_SIMULATION are incompatible" << endl; exit(1); }
    if (boolkimeraread && boolxyzread) {cout<< "Wacth out!! READ_POSITIONS_FROM_XYZ_FILE and READ_SYSTEM_FROM_KIMERA_FILE are incompatible" << endl; exit(1);}
    if (!ffdbool){control->set_fd_from_T();} //nothing to say, supposed ffd from fundamental frequency expresion with T=300K
    if (!ffrbool){control->set_fr(control->get_fd());}//nothing to say, , supposed ffp from fundamental frequency expresion with T=300K

    return 0;
}



long long int Data_reader::read_input_file_cell(string * filename_, Cell * cell)
{

    ifstream inputfile(filename_->c_str());

    string keyword;
    string value;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            inputfile >> keyword;

            if (keyword.compare("POSITION")==0)
            {
            string type_s;
            string x_s;
            string y_s;
            string z_s;
            string occupancy;

            double xvalue;
            double yvalue;
            double zvalue;
            double occupancyvalue;


            inputfile >> type_s;
            inputfile >> x_s;
            inputfile >> y_s;
            inputfile >> z_s;
            inputfile >> occupancy;

            try {
                xvalue=stof(x_s);
                yvalue=stof(y_s);
                zvalue=stof(z_s);
                occupancyvalue=stof(occupancy);
            }
            catch (const invalid_argument& ia) {
                cout << "Invalid argument in POSITION: " << ia.what() << endl;
                inputfile.close();
                exit(1);
            }

            if (occupancyvalue<0.0)
            {
                cout<<endl<<"Not valid number in POSITION"<<endl;
                inputfile.close();
                exit(1);
            }

            cell->compare_and_add(xvalue,yvalue,zvalue, type_s, occupancyvalue);

            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}


long long int Data_reader::read_input_file_topography(string * filename_, Box * box)
{
    bool readed=false;

   ifstream inputfile(filename_->c_str());


    string keyword;
    string value;

    string value2;
    string value3;
    string value4;
    string value5;
    string value6;
    string value7;
    string value8;
    string value9;
    string value10;
    string value11;



    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if (!readed) {inputfile >> keyword;}
            readed=false;

            if (keyword.compare("ADD_XY_DISLOCATION")==0)
            {
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XY_DISLOCATION: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XY_DISLOCATION"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in ADD_XY_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in ADD_XY_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }


                if (fromto)
                {
                    box->add_xy_heli_dislocation(x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_xy_heli_dislocation(x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                continue;
            }

            if (keyword.compare("ADD_XZ_DISLOCATION")==0)
            {
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XZ_DISLOCATION: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XZ_DISLOCATION"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in ADD_XZ_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in ADD_XZ_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_xz_heli_dislocation(x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_xz_heli_dislocation(x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                continue;
            }

            if (keyword.compare("ADD_YZ_DISLOCATION")==0)
            {
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_YZ_DISLOCATION: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_YZ_DISLOCATION"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in ADD_YZ_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in ADD_YZ_DISLOCATION: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_yz_heli_dislocation(y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_yz_heli_dislocation(y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                continue;
            }

            if (keyword.compare("DEFINE_AB_INSOLUBLE_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_INSOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_INSOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(dim_a, dim_b, dim_c, side_a, side_b, INSOLUBLE);
                cout<< "DEFINE_AB_INSOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AC_INSOLUBLE_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_INSOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_INSOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(dim_a, dim_b, dim_c, side_a, side_c, INSOLUBLE);
                cout<< "DEFINE_AC_INSOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_BC_INSOLUBLE_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;


                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_INSOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_INSOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(dim_a, dim_b, dim_c, side_b, side_c, INSOLUBLE);
                cout<< "DEFINE_BC_INSOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }


////////////////////////////
            if (keyword.compare("DEFINE_AB_SOLUBLE_CELLS")==0)  //DUPLICADO
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_SOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_SOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(dim_a, dim_b, dim_c, side_a, side_b, NORMAL);
                cout<< "DEFINE_AB_SOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AC_SOLUBLE_CELLS")==0)  //DUPLICADO
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_SOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_SOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(dim_a, dim_b, dim_c, side_a, side_c, NORMAL);
                cout<< "DEFINE_AC_SOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_BC_SOLUBLE_CELLS")==0)  //DUPLICADO
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;


                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_SOLUBLE_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_SOLUBLE_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(dim_a, dim_b, dim_c, side_b, side_c, NORMAL);
                cout<< "DEFINE_BC_SOLUBLE_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

 /////////////////////////////////////////////
 /**
            if (keyword.compare("DEFINE_AB_DISLOCATION_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_DISLOCATION_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_DISLOCATION_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(dim_a, dim_b, dim_c, side_a, side_b,INDISLOCATION);
                cout<< "DEFINE_AB_DISLOCATION_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;
            }

            if (keyword.compare("DEFINE_AC_DISLOCATION_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_DISLOCATION_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_DISLOCATION_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(dim_a, dim_b, dim_c, side_a, side_c, INDISLOCATION);
                cout<< "DEFINE_AC_DISLOCATION_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;
            }

            if (keyword.compare("DEFINE_BC_DISLOCATION_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;


                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_DISLOCATION_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_DISLOCATION_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(dim_a, dim_b, dim_c, side_b, side_c, INDISLOCATION);
                cout<< "DEFINE_BC_DISLOCATION_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;
            }

*/







            //////////////////
            if (keyword.compare("REMOVE_AB_PLANE_BY_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AB_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AB_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(dim_a, dim_b, dim_c, side_a, side_b, FREEPOSITION);
                cout<< "REMOVE_AB_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("REMOVE_AC_PLANE_BY_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AC_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AC_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(dim_a, dim_b, dim_c, side_a, side_c, FREEPOSITION);
                cout<< "REMOVE_AC_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_BC_PLANE_BY_CELLS")==0)
            {

                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_BC_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_BC_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(dim_a, dim_b, dim_c, side_b, side_c, FREEPOSITION);
                cout<< "REMOVE_BC_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_AB_PLANE_BY_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AB_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in ADD_AB_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(dim_a, dim_b, dim_c, side_a, side_b, NORMAL);
                cout<< "ADD_AB_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("ADD_AC_PLANE_BY_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AC_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_AC_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(dim_a, dim_b, dim_c, side_a, side_c, NORMAL);
                cout<< "ADD_AC_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_BC_PLANE_BY_CELLS")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_BC_PLANE_BY_CELLS: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_BC_PLANE_BY_CELLS"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(dim_a, dim_b, dim_c, side_b, side_c, NORMAL);
                cout<< "ADD_BC_PLANE_BY_CELLS at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_CELL")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CELL: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_CELL"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sym_cell(dim_a, dim_b, dim_c, NORMAL);
                cout<< "ADD_CELL at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<< endl;

                continue;
            }

            if (keyword.compare("REMOVE_CELL")==0)
            {

                long long int dim_a;
                long long int dim_b;
                long long int dim_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CELL: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CELL"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sym_cell(dim_a, dim_b, dim_c, FREEPOSITION);
                cout<< "REMOVE_CELL at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<< endl;

                continue;
            }

            if (keyword.compare("ADD_CUBE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CUBE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_CUBE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(x_value,y_value,z_value,radius, NORMAL);
                cout<< "ADD_CUBE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_CUBE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CUBE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CUBE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(x_value,y_value,z_value,radius,FREEPOSITION);
                cout<< "REMOVE_CUBE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_CUBE")==0)
            {

                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_CUBE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_CUBE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(x_value,y_value,z_value,radius,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_CUBE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_SPHERE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_SPHERE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_SPHERE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(x_value,y_value,z_value,radius,FREEPOSITION);
                cout<< "REMOVE_SPHERE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_SPHERE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_SPHERE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_SPHERE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(x_value,y_value,z_value,radius,NORMAL);
                cout<< "ADD_SPHERE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_SPHERE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_SPHERE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_SPHERE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(x_value,y_value,z_value,radius,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_SPHERE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_ELLIPSOID")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_ELLIPSOID"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz,NORMAL);
                cout<< "ADD_ELLIPSOID at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_ELLIPSOID")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_ELLIPSOID"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz,FREEPOSITION);
                cout<< "REMOVE_ELLIPSOID at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_ELLIPSOID")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);


                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_ELLIPSOID"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_ELLIPSOID at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("ADD_GENERAL_ELLIPSOID")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;


                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_GENERAL_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_type_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,NORMAL);
                cout<< "ADD_GENERAL_ELLIPSOID A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_GENERAL_ELLIPSOID")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_GENERAL_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_type_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,FREEPOSITION);
                cout<< "REMOVE_GENERAL_ELLIPSOID A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_GENERAL_ELLIPSOID")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_GENERAL_ELLIPSOID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_type_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_GENERAL_ELLIPSOID A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_SOLUBLE_PLANE")==0)  //DUPLICADO
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_SOLUBLE_PLANE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_SOLUBLE_PLANE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(A_value,B_value,C_value,D_value,distance,NORMAL);
                cout<< "DEFINE_SOLUBLE_PLANE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_PLANE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_PLANE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_PLANE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(A_value,B_value,C_value,D_value,distance,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_PLANE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_PLANE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_PLANE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_PLANE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(A_value,B_value,C_value,D_value,distance,FREEPOSITION);
                cout<< "REMOVE_PLANE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("ADD_PLANE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_PLANE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_PLANE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(A_value,B_value,C_value,D_value,distance,NORMAL);
                cout<< "ADD_PLANE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

        //TYPE SELECTIVE TOPOGRAPHY
            if (keyword.compare("ADD_XY_DISLOCATION_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XY_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XY_DISLOCATION_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in ADD_XY_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in ADD_XY_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_xy_heli_dislocation(atom_type,x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_xy_heli_dislocation(atom_type,x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_XZ_DISLOCATION_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XZ_DISLOCATION_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in ADD_XZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in ADD_XZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_xz_heli_dislocation(atom_type,x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_xz_heli_dislocation(atom_type, x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_YZ_DISLOCATION_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_YZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_YZ_DISLOCATION_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in ADD_YZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in ADD_YZ_DISLOCATION_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_yz_heli_dislocation(atom_type,y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_yz_heli_dislocation(atom_type,y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_TYPE to "<<atom_type << " at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("DEFINE_AB_INSOLUBLE_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_INSOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_INSOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b, INSOLUBLE);
                cout<< "DEFINE_AB_INSOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AC_INSOLUBLE_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_INSOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_INSOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c, INSOLUBLE);
                cout<< "DEFINE_AC_INSOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_BC_INSOLUBLE_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_INSOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_INSOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(atom_type, dim_a, dim_b, dim_c, side_b, side_c, INSOLUBLE);
                cout<< "DEFINE_BC_INSOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AB_SOLUBLE_CELLS_TO_TYPE")==0)  //DUPLICADO
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_SOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_SOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b, NORMAL);
                cout<< "DEFINE_AB_SOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AC_SOLUBLE_CELLS_TO_TYPE")==0)  //DUPLICADO
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_SOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_SOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c, NORMAL);
                cout<< "DEFINE_AC_SOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_BC_SOLUBLE_CELLS_TO_TYPE")==0)  //DUPLICADO
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_SOLUBLE_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_SOLUBLE_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(atom_type,dim_a, dim_b, dim_c, side_b, side_c, NORMAL);
                cout<< "DEFINE_BC_SOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }



            /////
            /**
            if (keyword.compare("DEFINE_AB_DISLOCATION_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AB_DISLOCATION_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AB_DISLOCATION_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b, INDISLOCATION);
                cout<< "DEFINE_AB_DISLOCATION_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("DEFINE_AC_DISLOCATION_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_AC_DISLOCATION_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_AC_DISLOCATION_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c, INDISLOCATION);
                cout<< "DEFINE_AC_SOLUBLE_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("DEFINE_BC_DISLOCATION_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_BC_DISLOCATION_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in DEFINE_BC_DISLOCATION_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(atom_type,dim_a, dim_b, dim_c, side_b, side_c, INDISLOCATION);
                cout<< "DEFINE_BC_DISLOCATION_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }
*/

            //////////////////////////////////////////////////////////////////

            if (keyword.compare("REMOVE_AB_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AB_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AB_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b, FREEPOSITION);
                cout<< "REMOVE_AB_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("REMOVE_AC_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AC_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AC_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c, FREEPOSITION);
                cout<< "REMOVE_AC_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_BC_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_BC_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_BC_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(atom_type,dim_a, dim_b, dim_c, side_b, side_c, FREEPOSITION);
                cout<< "REMOVE_BC_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_AB_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);

                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AB_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in ADD_AB_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xyplane(atom_type, dim_a, dim_b, dim_c, side_a, side_b, NORMAL);
                cout<< "ADD_AB_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("ADD_AC_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AC_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_AC_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c, NORMAL);
                cout<< "ADD_AC_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_BC_PLANE_BY_CELLS_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_BC_PLANE_BY_CELLS_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_BC_PLANE_BY_CELLS_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_cell_yzplane(atom_type,dim_a, dim_b, dim_c, side_b, side_c, NORMAL);
                cout<< "ADD_BC_PLANE_BY_CELLS_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_CELL_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CELL_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_CELL_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sym_cell(atom_type,dim_a, dim_b, dim_c, NORMAL);
                cout<< "ADD_CELL_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<< endl;

                continue;
            }

            if (keyword.compare("REMOVE_CELL_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CELL_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CELL_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sym_cell(atom_type,dim_a, dim_b, dim_c, FREEPOSITION);
                cout<< "REMOVE_CELL_TO_TYPE to "<<atom_type << " at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<< endl;

                continue;
            }

            if (keyword.compare("ADD_CUBE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CUBE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_CUBE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(atom_type,x_value,y_value,z_value,radius,NORMAL);
                cout<< "ADD_CUBE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_CUBE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);


                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CUBE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CUBE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(atom_type,x_value,y_value,z_value,radius,FREEPOSITION);
                cout<< "REMOVE_CUBE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                    continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_CUBE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_CUBE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }
                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_CUBE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_gap(atom_type,x_value,y_value,z_value,radius,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_CUBE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_SPHERE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_SPHERE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_SPHERE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(atom_type,x_value,y_value,z_value,radius,FREEPOSITION);
                cout<< "REMOVE_SPHERE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_SPHERE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_SPHERE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_SPHERE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(atom_type,x_value,y_value,z_value,radius,NORMAL);
                cout<< "ADD_SPHERE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_SPHERE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_SPHERE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_SPHERE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_sphere(atom_type,x_value,y_value,z_value,radius,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_SPHERE_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_ELLIPSOID_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz,NORMAL);
                cout<< "ADD_ELLIPSOID_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);


                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_ELLIPSOID_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz,FREEPOSITION);
                cout<< "REMOVE_ELLIPSOID_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_ELLIPSOID_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_ELLIPSOID_TO_TYPE to "<<atom_type << " at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("ADD_GENERAL_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_GENERAL_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_type_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,NORMAL);
                cout<< "ADD_GENERAL_ELLIPSOID_TO_TYPE "<< atom_type <<" A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_GENERAL_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);


                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_GENERAL_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                box->change_type_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,FREEPOSITION);
                cout<< "REMOVE_GENERAL_ELLIPSOID_TO_TYPE "<< atom_type <<" A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_GENERAL_ELLIPSOID_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_GENERAL_ELLIPSOID_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_type_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_GENERAL_ELLIPSOID_TO_TYPE "<< atom_type <<" A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("DEFINE_SOLUBLE_PLANE_TO_TYPE")==0)  //DUPLICADO
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_SOLUBLE_PLANE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_SOLUBLE_PLANE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance,NORMAL);
                cout<< "DEFINE_SOLUBLE_PLANE_TO_TYPE "<< atom_type <<" A: " <<A_value<<" B: "<<B_value<<" C: "<<C_value<<" D: "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("DEFINE_INSOLUBLE_PLANE_TO_TYPE")==0)
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in DEFINE_INSOLUBLE_PLANE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in DEFINE_INSOLUBLE_PLANE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance,INSOLUBLE);
                cout<< "DEFINE_INSOLUBLE_PLANE_TO_TYPE "<< atom_type <<" A: " <<A_value<<" B: "<<B_value<<" C: "<<C_value<<" D: "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_PLANE_TO_TYPE")==0)
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_PLANE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_PLANE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance,FREEPOSITION);
                cout<< "REMOVE_PLAN_TO_TYPE "<< atom_type <<" A: " <<A_value<<" B: "<<B_value<<" C: "<<C_value<<" D: "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("ADD_PLANE_TO_TYPE")==0)
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_PLANE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_PLANE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_type_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance,NORMAL);
                cout<< "ADD_PLANE_TO_TYPE "<< atom_type <<" A: " <<A_value<<" B: "<<B_value<<" C: "<<C_value<<" D: "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}


 long long int Data_reader::read_modify_element(string * filename_, Box * box)
 {
    bool readed=false;

    ifstream inputfile(filename_->c_str());

    string keyword;
    string value;
    string value2;
    string value3;
    string value4;
    string value5;
    string value6;
    string value7;
    string value8;
    string value9;
    string value10;
    string value11;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if (!readed) {inputfile >> keyword;}
            readed=false;

            if (keyword.compare("CHANGE_AB_PLANE_TYPE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_AB_PLANE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_AB_PLANE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_xyplane(dim_a, dim_b, dim_c, side_a, side_b,target_element);
                cout<< "CHANGE_AB_PLANE_TYPE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << " to "<< target_element <<endl;

                continue;
            }

            if (keyword.compare("CHANGE_AC_PLANE_TYPE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_AC_PLANE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_AC_PLANE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_xzplane(dim_a, dim_b, dim_c, side_a, side_c,target_element);
                cout<< "CHANGE_AC_PLANE_TYPE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << " to "<< target_element << endl;

                continue;
            }
            if (keyword.compare("CHANGE_BC_PLANE_TYPE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_BC_PLANE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_BC_PLANE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_yzplane(dim_a, dim_b, dim_c, side_b, side_c,target_element);
                cout<< "CHANGE_BC_PLANE_TYPE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << " to "<< target_element << endl;

                continue;
            }
            if (keyword.compare("CHANGE_CUBE_TYPE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_CUBE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_CUBE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_gap(x_value,y_value,z_value,radius,target_element);
                cout<< "CHANGE_CUBE_TYPE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<< " to "<< target_element <<endl;

                continue;
            }

            if (keyword.compare("CHANGE_XY_DISLOCATION_TYPE")==0)
            {
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_XY_DISLOCATION_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_XY_DISLOCATION_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in CHANGE_XY_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in CHANGE_XY_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_xy_heli_dislocation(x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius,target_element);
                    cout<< "CHANGE_XY_DISLOCATION_TYPE at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius << " to "<< target_element <<endl;
                }
                else
                {
                    box->change_element_xy_heli_dislocation(x_value,y_value,angle_xz,angle_yz,radius,target_element);
                    cout<< "CHANGE_XY_DISLOCATION_TYPE at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius << " to "<< target_element <<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_XZ_DISLOCATION_TYPE")==0)
            {
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_XZ_DISLOCATION_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_XZ_DISLOCATION_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in CHANGE_XZ_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in CHANGE_XZ_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_xz_heli_dislocation(x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius,target_element);
                    cout<< "CHANGE_XZ_DISLOCATION_TYPE at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<< " to "<< target_element<<endl;
                }
                else
                {
                    box->change_element_xz_heli_dislocation(x_value,z_value,angle_xy,angle_zy,radius,target_element);
                    cout<< "CHANGE_XZ_DISLOCATION_TYPE at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<< " to "<< target_element<<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_YZ_DISLOCATION_TYPE")==0)
            {
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_YZ_DISLOCATION_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_YZ_DISLOCATION_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in CHANGE_YZ_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in CHANGE_YZ_DISLOCATION_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_yz_heli_dislocation(y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius,target_element);
                    cout<< "CHANGE_YZ_DISLOCATION_TYPE at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<< " to "<< target_element<<endl;
                }
                else
                {
                    box->change_element_yz_heli_dislocation(y_value,z_value,angle_yx,angle_zx,radius,target_element);
                    cout<< "CHANGE_YZ_DISLOCATION_TYPE at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<< " to "<< target_element<<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_PLANE_TYPE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_PLANE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_PLANE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_plane_distance(A_value,B_value,C_value,D_value,distance, target_element);
                cout<< "CHANGE_PLANE_TYPE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<< " to "<< target_element<<endl;

                continue;
            }

           if (keyword.compare("CHANGE_SPHERE_TYPE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_SPHERE_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_SPHERE_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_sphere(x_value,y_value,z_value,radius,target_element);
                cout<< "CHANGE_SPHERE_TYPE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<" to "<< target_element<<endl;

                continue;
            }

            if (keyword.compare("CHANGE_ELLIPSOID_TYPE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;
                string target_element;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_ELLIPSOID_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_ELLIPSOID_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz,target_element);
                cout<< "CHANGE_ELLIPSOID_TYPE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<" to "<< target_element<<endl;

                continue;
            }

            if (keyword.compare("CHANGE_GENERAL_ELLIPSOID_TYPE")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;
                string target_element;


                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;
                inputfile >> target_element;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_GENERAL_ELLIPSOID_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_element_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,target_element);
                cout<< "CHANGE_GENERAL_ELLIPSOID_TYPE A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<" to "<< target_element<<endl;

                continue;
            }

            /////////////////////////////////////////////////
            ///////////////////// TO TYPE ///////////////////
            /////////////////////////////////////////////////

            if (keyword.compare("CHANGE_AB_PLANE_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_AB_PLANE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_AB_PLANE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b,target_element);
                cout<< "CHANGE_AB_PLANE_TYPE_TO_TYPE to "<< atom_type<<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << " to "<< target_element <<endl;

                continue;
            }

            if (keyword.compare("CHANGE_AC_PLANE_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_AC_PLANE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_AC_PLANE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c,target_element);
                cout<< "CHANGE_AC_PLANE_TYPE_TO_TYPE to "<< atom_type<<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << " to "<< target_element << endl;

                continue;
            }

            if (keyword.compare("CHANGE_BC_PLANE_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_BC_PLANE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in CHANGE_BC_PLANE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_yzplane(atom_type, dim_a, dim_b, dim_c, side_b, side_c,target_element);
                cout<< "CHANGE_BC_PLANE_TYPE_TO_TYPE to "<< atom_type<<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << " to "<< target_element << endl;

                continue;
            }
            if (keyword.compare("CHANGE_CUBE_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_CUBE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_CUBE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_gap(atom_type,x_value,y_value,z_value,radius,target_element);
                cout<< "CHANGE_CUBE_TYPE_TO_TYPE  to "<< atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<< " to "<< target_element <<endl;

                continue;
            }

            if (keyword.compare("CHANGE_XY_DISLOCATION_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_XY_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_XY_DISLOCATION_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in CHANGE_XY_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in CHANGE_XY_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_xy_heli_dislocation(atom_type,x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius,target_element);
                    cout<< "CHANGE_XY_DISLOCATION_TYPE_TO_TYPE to "<< atom_type<<"  at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius << " to "<< target_element <<endl;
                }
                else
                {
                    box->change_element_xy_heli_dislocation(atom_type,x_value,y_value,angle_xz,angle_yz,radius,target_element);
                    cout<< "CHANGE_XY_DISLOCATION_TYPE_TO_TYPE  to "<< atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius << " to "<< target_element <<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE")==0)
            {
                string target_type;
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;
                string target_element;

                inputfile >> target_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_xz_heli_dislocation(target_type, x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius,target_element);
                    cout<< "CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE  to "<< target_type<<"  at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<< " to "<< target_element<<endl;
                }
                else
                {
                    box->change_element_xz_heli_dislocation(target_type, x_value,z_value,angle_xy,angle_zy,radius,target_element);
                    cout<< "CHANGE_XZ_DISLOCATION_TYPE_TO_TYPE   to "<< target_type<<" at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<< " to "<< target_element<<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE")==0)
            {
                string target_type;
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;
                string target_element;

                inputfile >> target_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> target_element;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->change_element_yz_heli_dislocation(target_type,y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius,target_element);
                    cout<< "CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE  to "<< target_type<<" at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<< " to "<< target_element<<endl;
                }
                else
                {
                    box->change_element_yz_heli_dislocation(target_type,y_value,z_value,angle_yx,angle_zx,radius,target_element);
                    cout<< "CHANGE_YZ_DISLOCATION_TYPE_TO_TYPE  to "<< target_type<<" at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<< " to "<< target_element<<endl;
                }

                continue;
            }

            if (keyword.compare("CHANGE_PLANE_TYPE_TO_TYPE")==0)
            {
                string target_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;
                string target_element;


                inputfile >> target_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> target_element;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_PLANE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_PLANE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_plane_distance(target_type,A_value,B_value,C_value,D_value,distance, target_element);
                cout<< "CHANGE_PLANE_TYPE_TO_TYPE to "<< target_type<<" with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<< " to "<< target_element<<endl;

                continue;
            }

           if (keyword.compare("CHANGE_SPHERE_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_SPHERE_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_SPHERE_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_sphere(atom_type,x_value,y_value,z_value,radius,target_element);
                cout<< "CHANGE_SPHERE_TYPE_TO_TYPE  to "<< atom_type <<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<" to "<< target_element<<endl;

                continue;
            }
            if (keyword.compare("CHANGE_ELLIPSOID_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;
                string target_element;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> target_element;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_ELLIPSOID_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in CHANGE_ELLIPSOID_TYPE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->change_element_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz,target_element);
                cout<< "CHANGE_ELLIPSOID_TYPE_TO_TYPE to "<< atom_type <<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<" to "<< target_element<<endl;

                continue;
            }


            if (keyword.compare("CHANGE_GENERAL_ELLIPSOID_TYPE_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;
                string target_element;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;
                inputfile >> target_element;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CHANGE_GENERAL_ELLIPSOID_TYPE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->change_element_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_,target_element);
                cout<< "CHANGE_GENERAL_ELLIPSOID_TYPE_TO_TYPE to "<< atom_type <<" A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<" to "<< target_element<<endl;

                continue;
            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}



 long long int Data_reader::read_remove_add_from_surface(string * filename_, Box * box)
 {
    bool readed=false;

    ifstream inputfile(filename_->c_str());

    string keyword;
    string value;
    string value2;
    string value3;
    string value4;
    string value5;
    string value6;
    string value7;
    string value8;
    string value9;
    string value10;
    string value11;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if (!readed) {inputfile >> keyword;}
            readed=false;

            if (keyword.compare("REMOVE_AB_PLANE_FROM_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AB_PLANE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AB_PLANE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_xyplane(dim_a, dim_b, dim_c, side_a, side_b);
                cout<< "REMOVE_AB_PLANE_FROM_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("REMOVE_AC_PLANE_FROM_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AC_PLANE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AC_PLANE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_xzplane(dim_a, dim_b, dim_c, side_a, side_c);
                cout<< "REMOVE_AC_PLANE_FROM_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_BC_PLANE_FROM_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_BC_PLANE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_BC_PLANE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_yzplane(dim_a, dim_b, dim_c, side_b, side_c);
                cout<< "REMOVE_BC_PLANE_FROM_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }
            if (keyword.compare("ADD_AB_PLANE_TO_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AB_PLANE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in ADD_AB_PLANE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_xyplane(dim_a, dim_b, dim_c, side_a, side_b);
                cout<< "ADD_AB_PLANE_TO_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("ADD_AC_PLANE_TO_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AC_PLANE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_AC_PLANE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_xzplane(dim_a, dim_b, dim_c, side_a, side_c);
                cout<< "ADD_AC_PLANE_TO_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_BC_PLANE_TO_SURFACE")==0)
            {
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_BC_PLANE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_BC_PLANE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_yzplane(dim_a, dim_b, dim_c, side_b, side_c);
                cout<< "ADD_BC_PLANE_TO_SURFACE at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_CUBE_FROM_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CUBE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CUBE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_gap(x_value,y_value,z_value,radius);
                cout<< "REMOVE_CUBE_FROM_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_CUBE_TO_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CUBE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_CUBE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_gap(x_value,y_value,z_value,radius);
                cout<< "ADD_CUBE_TO_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_XY_DISLOCATION_FROM_SURFACE")==0)
            {
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_XY_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_XY_DISLOCATION_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in REMOVE_XY_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in REMOVE_XY_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_xy_heli_dislocation(x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "REMOVE_XY_DISLOCATION_FROM_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_xy_heli_dislocation(x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "REMOVE_XY_DISLOCATION_FROM_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_XY_DISLOCATION_TO_SURFACE")==0)
            {
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XY_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XY_DISLOCATION_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in ADD_XY_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in ADD_XY_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_xy_heli_dislocation(x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_xy_heli_dislocation(x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_XZ_DISLOCATION_FROM_SURFACE")==0)
            {
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_XZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_XZ_DISLOCATION_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in REMOVE_XZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in REMOVE_XZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_xz_heli_dislocation(x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "REMOVE_XZ_DISLOCATION_FROM_SURFACE at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_xz_heli_dislocation(x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "REMOVE_XZ_DISLOCATION_FROM_SURFACE at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_XZ_DISLOCATION_TO_SURFACE")==0)
            {
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XZ_DISLOCATION_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in ADD_XZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in ADD_XZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_xz_heli_dislocation(x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_SURFACE at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_xz_heli_dislocation(x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_SURFACE at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_YZ_DISLOCATION_FROM_SURFACE")==0)
            {
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_YZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_YZ_DISLOCATION_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in REMOVE_YZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in REMOVE_YZ_DISLOCATION_FROM_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_yz_heli_dislocation(y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "REMOVE_YZ_DISLOCATION_FROM_SURFACE at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_yz_heli_dislocation(y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "REMOVE_YZ_DISLOCATION_FROM_SURFACE at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_YZ_DISLOCATION_TO_SURFACE")==0)
            {
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_YZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_YZ_DISLOCATION_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in ADD_YZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in ADD_YZ_DISLOCATION_TO_SURFACE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_yz_heli_dislocation(y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_SURFACE at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_yz_heli_dislocation(y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_SURFACE at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_PLANE_FROM_SURFACE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_PLANE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_PLANE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_plane_distance(A_value,B_value,C_value,D_value,distance);
                cout<< "REMOVE_PLANE_FROM_SURFACE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("ADD_PLANE_TO_SURFACE")==0)
            {
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_PLANE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_PLANE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_plane_distance(A_value,B_value,C_value,D_value,distance);
                cout<< "ADD_PLANE_TO_SURFACE with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;
                continue;
            }

            if (keyword.compare("REMOVE_SPHERE_FROM_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_SPHERE_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_SPHERE_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_sphere(x_value,y_value,z_value,radius);
                cout<< "REMOVE_SPHERE_FROM_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;
                continue;
            }

            if (keyword.compare("ADD_SPHERE_TO_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_SPHERE_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_SPHERE_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_sphere(x_value,y_value,z_value,radius);
                cout<< "ADD_SPHERE_TO_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;
                continue;
            }

            if (keyword.compare("REMOVE_ELLIPSOID_FROM_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_ELLIPSOID_FROM_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_ELLIPSOID_FROM_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz);
                cout<< "REMOVE_ELLIPSOID_FROM_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;
                continue;
            }

            if (keyword.compare("ADD_ELLIPSOID_TO_SURFACE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_ELLIPSOID_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_ELLIPSOID_TO_SURFACE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_ellipsoid(x_value,y_value,z_value,radiusx,radiusy,radiusz);
                cout<< "ADD_ELLIPSOID_TO_SURFACE at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_GENERAL_ELLIPSOID_TO_SURFACE")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;


                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_GENERAL_ELLIPSOID_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->remove_from_surface_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_);
                cout<< "REMOVE_GENERAL_ELLIPSOID_TO_SURFACE A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("ADD_GENERAL_ELLIPSOID_TO_SURFACE")==0)
            {
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;


                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_GENERAL_ELLIPSOID_TO_SURFACE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->add_to_surface_general_ellipsoid(A_,B_,C_,D_,E_,F_,G_,H_,J_,K_);
                cout<< "ADD_GENERAL_ELLIPSOID_TO_SURFACE A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }



//////////////////////////////////TO_TYPE
//////////////////////////////////TO_TYPE
//////////////////////////////////TO_TYPE
//////////////////////////////////TO_TYPE


            if (keyword.compare("REMOVE_AB_PLANE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AB_PLANE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AB_PLANE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b);
                cout<< "REMOVE_AB_PLANE_FROM_SURFACE_TO_TYPE to "<< atom_type <<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("REMOVE_AC_PLANE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_AC_PLANE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_AC_PLANE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_xzplane(atom_type, dim_a, dim_b, dim_c, side_a, side_c);
                cout<< "REMOVE_AC_PLANE_FROM_SURFACE_TO_TYPE "<< atom_type <<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_BC_PLANE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_BC_PLANE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in REMOVE_BC_PLANE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_yzplane(atom_type, dim_a, dim_b, dim_c, side_b, side_c);
                cout<< "REMOVE_BC_PLANE_FROM_SURFACE_TO_TYPE to "<< atom_type <<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }
            if (keyword.compare("ADD_AB_PLANE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_b;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_b=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AB_PLANE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_b<0)
                {
                    cout<<endl<<"Not valid number in ADD_AB_PLANE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_xyplane(atom_type,dim_a, dim_b, dim_c, side_a, side_b);
                cout<< "ADD_AB_PLANE_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_b "
                << side_b << endl;

                continue;
            }

            if (keyword.compare("ADD_AC_PLANE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_a;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_a=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_AC_PLANE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_a<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_AC_PLANE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_xzplane(atom_type,dim_a, dim_b, dim_c, side_a, side_c);
                cout<< "ADD_AC_PLANE_TO_SURFACE_TO_TYPE to "<< atom_type <<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_a "<< side_a << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("ADD_BC_PLANE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                long long int dim_a;
                long long int dim_b;
                long long int dim_c;
                long long int side_b;
                long long int side_c;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    dim_a=stoll(value);
                    dim_b=stoll(value2);
                    dim_c=stoll(value3);
                    side_b=stoll(value4);
                    side_c=stoll(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_BC_PLANE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (dim_a<0 || dim_b<0 || dim_c<0 || side_b<0 || side_c<0)
                {
                    cout<<endl<<"Not valid number in ADD_BC_PLANE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_yzplane(atom_type,dim_a, dim_b, dim_c, side_b, side_c);
                cout<< "ADD_BC_PLANE_TO_SURFACE_TO_TYPE to "<<atom_type<<" at " <<dim_a<<"  "<<dim_b<<"  "<<dim_c<<" and side_b "<< side_b << " side_c "
                << side_c << endl;

                continue;
            }

            if (keyword.compare("REMOVE_CUBE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_CUBE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_CUBE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_gap(atom_type,x_value,y_value,z_value,radius);
                cout<< "REMOVE_CUBE_FROM_SURFACE_TO_TYPE to "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_CUBE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_CUBE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_CUBE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_gap(atom_type,x_value,y_value,z_value,radius);
                cout<< "ADD_CUBE_TO_SURFACE_TO_TYPE to "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with side "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_xy_heli_dislocation(atom_type, x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_xy_heli_dislocation(atom_type, x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "REMOVE_XY_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double y_value;
                double z_bot;
                double z_top;
                double radius;
                double angle_xz=90.0;
                double angle_yz=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Z_TO_Z")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        z_bot=stof(value);
                        z_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Z_TO_Z in ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }


                if (keyword.compare("ANGLE_XZ_ANGLE_YZ")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xz=stof(value);
                        angle_yz=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XZ_ANGLE_YZ in ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_xy_heli_dislocation(atom_type,x_value,y_value,z_bot,z_top,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<" with z_bot "<<z_bot<<" z_top "<<z_top
                    <<" angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_xy_heli_dislocation(atom_type,x_value,y_value,angle_xz,angle_yz,radius);
                    cout<< "ADD_XY_DISLOCATION_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<
                    " with angle_xz "<<angle_xz<<" angle_yz "<<angle_yz<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_xz_heli_dislocation(atom_type,x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_xz_heli_dislocation(atom_type,x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "REMOVE_XZ_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double x_value;
                double z_value;
                double y_bot;
                double y_top;
                double radius;
                double angle_xy=90.0;
                double angle_zy=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;


                try {
                    x_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_Y_TO_Y")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        y_bot=stof(value);
                        y_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_Y_TO_Y in ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;

                }

                if (keyword.compare("ANGLE_XY_ANGLE_ZY")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_xy=stof(value);
                        angle_zy=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_XY_ANGLE_ZY in ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_xz_heli_dislocation(atom_type,x_value,z_value,y_bot,y_top,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<z_value<<" with y_bot "<<y_bot<<" y_top "<<y_top
                    <<" angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_xz_heli_dislocation(atom_type,x_value,z_value,angle_xy,angle_zy,radius);
                    cout<< "ADD_XZ_DISLOCATION_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<z_value<<
                    " with angle_xy "<<angle_xy<<" angle_zy "<<angle_zy<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->remove_from_surface_yz_heli_dislocation(atom_type,y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<"  at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->remove_from_surface_yz_heli_dislocation(atom_type,y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "REMOVE_YZ_DISLOCATION_FROM_SURFACE_TO_TYPE "<<atom_type<<" at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                bool fromto=false;
                double y_value;
                double z_value;
                double x_bot;
                double x_top;
                double radius;
                double angle_yx=90.0;
                double angle_zx=90.0;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;

                try {
                    y_value=stof(value);
                    z_value=stof(value2);
                    radius=stof(value3);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }


                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }

                inputfile >> keyword;
                readed=true;

                if (keyword.compare("FROM_X_TO_X")==0)
                {
                    fromto=true;
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        x_bot=stof(value);
                        x_top=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in FROM_X_TO_X in ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    inputfile >> keyword;
                    readed=true;
                }

                if (keyword.compare("ANGLE_YX_ANGLE_ZX")==0)
                {
                    inputfile >> value;
                    inputfile >> value2;

                    try {
                        angle_yx=stof(value);
                        angle_zx=stof(value2);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in ANGLE_YX_ANGLE_ZX in ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                        inputfile.close();
                        exit(1);
                    }
                    readed=false;
                }

                if (fromto)
                {
                    box->add_to_surface_yz_heli_dislocation(atom_type,y_value,z_value,x_bot,x_top,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE "<<atom_type<<"  at " <<y_value<<"  "<<z_value<<" with x_bot "<<x_bot<<" x_top "<<x_top
                    <<" angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }
                else
                {
                    box->add_to_surface_yz_heli_dislocation(atom_type,y_value,z_value,angle_yx,angle_zx,radius);
                    cout<< "ADD_YZ_DISLOCATION_TO_SURFACE_TO_TYPE  "<<atom_type<<" at " <<y_value<<"  "<<z_value<<"  "<<
                    " with angle_yx "<<angle_yx<<" angle_zx "<<angle_zx<<" and radius "<< radius<<endl;
                }

                continue;
            }

            if (keyword.compare("REMOVE_PLANE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_PLANE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_PLANE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance);
                cout<< "REMOVE_PLANE_FROM_SURFACE_TO_TYPE "<<atom_type<<"  with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("ADD_PLANE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double A_value;
                double B_value;
                double C_value;
                double D_value;
                double distance;

                inputfile >> atom_type;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;

                try {
                    A_value=stof(value);
                    B_value=stof(value2);
                    C_value=stof(value3);
                    D_value=stof(value4);
                    distance=stof(value5);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_PLANE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (distance<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_PLANE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_plane_distance(atom_type,A_value,B_value,C_value,D_value,distance);
                cout<< "ADD_PLANE_TO_SURFACE_TO_TYPE  "<<atom_type<<" with A " <<A_value<<" B "<<B_value<<" C "<<C_value<<" D "<<D_value<<" and distance "<< distance<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_SPHERE_FROM_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_SPHERE_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_SPHERE_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_sphere(atom_type,x_value,y_value,z_value,radius);
                cout<< "REMOVE_SPHERE_FROM_SURFACE_TO_TYPE  "<<atom_type<<"  at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("ADD_SPHERE_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radius;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radius=stof(value4);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_SPHERE_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radius<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_SPHERE_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_sphere(atom_type,x_value,y_value,z_value,radius);
                cout<< "ADD_SPHERE_TO_SURFACE_TO_TYPE  "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<" with radius "<<radius<<endl;

                continue;
            }

            if (keyword.compare("REMOVE_ELLIPSOID_FROM_SURFACE_TO_TYPE")==0)
            {
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;
                string atom_type;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_ELLIPSOID_FROM_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in REMOVE_ELLIPSOID_FROM_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->remove_from_surface_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz);
                cout<< "REMOVE_ELLIPSOID_FROM_SURFACE_TO_TYPE "<<atom_type<<"  at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("ADD_ELLIPSOID_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double x_value;
                double y_value;
                double z_value;
                double radiusx;
                double radiusy;
                double radiusz;

                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;

                try {
                    x_value=stof(value);
                    y_value=stof(value2);
                    z_value=stof(value3);
                    radiusx=stof(value4);
                    radiusy=stof(value5);
                    radiusz=stof(value6);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_ELLIPSOID_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (radiusx<0.0 || radiusy<0.0 || radiusz<0.0)
                {
                    cout<<endl<<"Not valid number in ADD_ELLIPSOID_TO_SURFACE_TO_TYPE"<<endl;
                    inputfile.close();
                    exit(1);
                }
                box->add_to_surface_ellipsoid(atom_type,x_value,y_value,z_value,radiusx,radiusy,radiusz);
                cout<< "ADD_ELLIPSOID_TO_SURFACE_TO_TYPE "<<atom_type<<" at " <<x_value<<"  "<<y_value<<"  "<<z_value<<
                " with radius x "<<radiusx <<" radius y "<<radiusy <<" and radius z "<<radiusz <<endl;

                continue;
            }

            if (keyword.compare("REMOVE_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in REMOVE_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->remove_from_surface_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_);
                cout<< "REMOVE_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE  "<<atom_type<<"  A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }

            if (keyword.compare("ADD_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE")==0)
            {
                string atom_type;
                double A_;
                double B_;
                double C_;
                double D_;
                double E_;
                double F_;
                double G_;
                double H_;
                double J_;
                double K_;


                inputfile >> atom_type;
                inputfile >> value;
                inputfile >> value2;
                inputfile >> value3;
                inputfile >> value4;
                inputfile >> value5;
                inputfile >> value6;
                inputfile >> value7;
                inputfile >> value8;
                inputfile >> value9;
                inputfile >> value10;

                try {
                    A_=stof(value);
                    B_=stof(value2);
                    C_=stof(value3);
                    D_=stof(value4);
                    E_=stof(value5);
                    F_=stof(value6);
                    G_=stof(value7);
                    H_=stof(value8);
                    J_=stof(value9);
                    K_=stof(value10);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in ADD_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                box->add_to_surface_general_ellipsoid(atom_type,A_,B_,C_,D_,E_,F_,G_,H_,J_,K_);
                cout<< "ADD_GENERAL_ELLIPSOID_TO_SURFACE_TO_TYPE  "<<atom_type<<" A: " << A_ << "B: "<< B_ << "C: "<< C_ << "D: "<< D_ <<
                "E: "<< E_ << "F: "<< F_ << "G: "<< G_ << "H: "<< H_ << "J: "<< J_ <<   "K: "<< K_ <<endl;

                continue;
            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}


long long int Data_reader::read_input_file_mass(string * filename_, Box * box)
{
    ifstream inputfile(filename_->c_str());

    string keyword;
    string value;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            inputfile >> keyword;

            if (keyword.compare("SET_MASS")==0)
            {
            string type_s;
            string mass_s;

            double massvalue;

            inputfile >> type_s;
            inputfile >> mass_s;

            try {
                massvalue=stof(mass_s);
            }
            catch (const invalid_argument& ia) {
                cout << "Invalid argument in SET_MASS: " << ia.what() << endl;
                inputfile.close();
                exit(1);
            }

            if (massvalue<0.0)
            {
                cout<<endl<<"Not valid number in SET_MASS"<<endl;
                inputfile.close();
                exit(1);
            }
                box->set_mass(type_s, massvalue);
                cout<< "SET_MASS to "<< type_s <<" of "<< massvalue << " u"<<endl;
            }
        }
    }
    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}

long long int Data_reader::read_kimera_box(string * filename_, Box * box)
{
    bool readed=false;

    ifstream inputfile(filename_->c_str());


    string keyword;
    string value;

    string value2;
    string value3;
    string value4;
    string value5;
    string value6;

    if (inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            if (!readed) {inputfile >> keyword;}
            readed=false;

            if (keyword.compare("CELL")==0)
            {
                long long int cell_id;

                inputfile>>value;

                try {
                    cell_id=stoll(value);
                }
                catch (const invalid_argument& ia) {
                    cout << "Invalid argument in CELL_ID: " << ia.what() << endl;
                    inputfile.close();
                    exit(1);
                }

                if (cell_id<0)
                {
                    cout<<endl<<"Not valid number in CELL_ID"<<endl;
                    inputfile.close();
                    exit(1);
                }

                Cell* cell1=box->add_cell_to_box(cell_id);

                inputfile>>keyword;
                readed=true;

                 long long int id;
                string atom_type;
                long long int type;
                double pos_x;
                double pos_y;
                double pos_z;
                long long int insurface;

                while(keyword.compare("CELL")!=0)
                {
                     if (keyword.compare("NEIGHBORS")==0) break;

                    inputfile>>value;
                    inputfile>>value2;
                    inputfile>>value3;
                    inputfile>>value4;
                    inputfile>>value5;
                    inputfile>>value6;
                    try {
                        id=stoll(keyword);
                        atom_type=value;
                        type=stoll(value2);
                        pos_x=stod(value3);
                        pos_y=stod(value4);
                        pos_z=stod(value5);
                        insurface=stoll(value6);
                    }
                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in adding ATOM: " << ia.what() << endl;
                        cout << "Atom id  "<< id <<endl;
                        inputfile.close();
                        exit(1);
                    }
                    Atom * atomtoenter=box->add_atom(id,atom_type,type,pos_x,pos_y,pos_z);
                    if (insurface==1) {box->add_surface_atom(atomtoenter); }
                    if (insurface==0) {atomtoenter->set_insurface(false);}
                    cell1->add_atom(atomtoenter);

                    inputfile>>keyword;
                    readed=true;
                }
            }

            if (keyword.compare("NEIGHBORS")==0)
            {
                string line;
                string value;
                string value_neigh;
                string value_distance;
                 long long int atomid;
                 long long int neigh_id;
                double distance;

                getline(inputfile,line);

                while(!inputfile.eof())
                {
                    getline(inputfile,line);

                    stringstream iline(line);

                    iline>>value;

                    if (value.compare("AFFECTED")==0){keyword="AFFECTED";break;}
                    if (value.compare("LINKED")==0){keyword="LINKED";break;}

                    iline>>value_neigh;
                    iline>>value_distance;

                    try {
                        atomid=stoll(value);
                    }

                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in NEIGHBOURS: " << ia.what() << endl;
                        cout << "wrong id "<< atomid <<endl;
                        inputfile.close();
                        exit(1);
                    }

                    Atom * targetatom=box->select_atom_by_id(atomid);

                    if (!iline.eof()){
                        try {

                            neigh_id=stoll(value_neigh);
                            distance=stof(value_distance);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in NEIGHBOURS: " << ia.what() << endl;

                            inputfile.close();
                            exit(1);
                        }

                        targetatom->add_neighbour(box->select_atom_by_id(neigh_id));
                        targetatom->add_distance_to_neighbours(distance);

                        while(true)
                        {
                            iline>>value_neigh;
                            iline>>value_distance;

                            if(!iline.eof())
                            {
                                try {
                                    neigh_id=stoll(value_neigh);
                                    distance=stof(value_distance);
                                }

                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in NEIGHBOURS: " << ia.what() << endl;
                                    cout << "id del que falla "<< atomid <<endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                targetatom->add_neighbour(box->select_atom_by_id(neigh_id));
                                targetatom->add_distance_to_neighbours(distance);

                            }
                            else break;
                        }
                    }
                    targetatom->set_neigh_record();
                }
            }

            if (keyword.compare("AFFECTED")==0)
            {
                string line;
                string value;
                string value_affected;
                 long long int atomid;
                 long long int affected_id;

                while(!inputfile.eof())
                {
                    getline(inputfile,line);

                    stringstream iline(line);

                    iline>>value;

                    if (value.compare("LINKED")==0){keyword="LINKED";break;}

                    iline>>value_affected;

                    try {
                        atomid=stoll(value);
                    }
                    catch (const invalid_argument& ia) {
                        inputfile.close();
                        return 0;
                    }

                    Atom * targetatom=box->select_atom_by_id(atomid);

                    if (!iline.eof()){
                        try {
                            affected_id=stoll(value_affected);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in AFFECTED: " << ia.what() << endl;
                            inputfile.close();
                            return 0;
                        }

                        targetatom->add_affected(box->select_atom_by_id(affected_id));

                        while(true)
                        {
                            iline>>value_affected;

                            if(!iline.eof())
                            {
                                try {
                                    affected_id=stoll(value_affected);
                                }

                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in AFFECTED: " << ia.what() << endl;
                                    inputfile.close();
                                    exit(1);
                                }

                                targetatom->add_affected(box->select_atom_by_id(affected_id));
                            }
                            else break;
                        }
                    }
                }
            }

            if (keyword.compare("LINKED")==0)
            {
                bool readed2=false;
                string line;
                string value;
                string valuetype;
                string value_id;

                string Linkingword;

                long long int atom_id;
                long long int link_type;
                long long int link_id;

                //getline(inputfile,line);

                while(!inputfile.eof())
                {
                    long long int towhom=0;

                    getline(inputfile,line);

                    stringstream iline(line);

                    iline>>value;

                    try {
                        atom_id=stoll(value);
                    }

                    catch (const invalid_argument& ia) {
                        cout << "Invalid argument in LINKED: " << ia.what() << endl;
                        cout << "wrong id: "<< value <<endl;
                        inputfile.close();
                        exit(1);
                    }

                    Atom * targetatom=box->select_atom_by_id(atom_id);

                    while(!iline.eof())
                    {
                        if (!readed2){iline>>Linkingword;}
                        readed2=false;
                        iline>>valuetype;


                        try {

                            link_type=stoll(valuetype);
                        }

                        catch (const invalid_argument& ia) {
                            cout << "Invalid argument in LINKED: " << ia.what() << endl;

                            inputfile.close();
                            exit(1);
                        }

                        targetatom->add_linked_type(link_type);
                        targetatom->add_void_vector_in_linked();

                        //despues, si estamos annadiendo linkeds, tiene que haber una L_




                        while(!iline.eof())
                        {
                            iline>>value_id;


                            if (value_id.compare("L")==0)
                            {
                                readed2=true;
                                towhom++;
                                break;
                            }

                            if(!iline.eof())
                            {
                                try {
                                    link_id=stoll(value_id);
                                }

                                catch (const invalid_argument& ia) {
                                    cout << "Invalid argument in LINKED: " << ia.what() << endl;
                                    cout << "wrong id: "<< value_id <<endl;
                                    inputfile.close();
                                    exit(1);
                                }


                                targetatom->add_linked_to(box->select_atom_by_id(link_id),towhom);

                            }
                            else break;
                        }
                    }
                }
            }
        }
    }


    else
    {
        cout<< "Not available input file" <<endl;
        inputfile.close();
        exit(1);
    }
    return 0;
}


