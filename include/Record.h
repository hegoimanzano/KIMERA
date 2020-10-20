#ifndef RECORD_H
#define RECORD_H


#include <string>

using namespace std;

class Record
{
    public:
        Record(string type_, double distance_);
        virtual ~Record();

        long long int get_number();
        void set_number(long long int number_);
        long long int increment_number_by_1();
        long long int reduce_number_by_1();


        string get_type();
        void set_type(string type_);

        double get_distance();
        void set_distance(double distance_);


    protected:

    private:

        long long int number;
        string type;
        double distance;
};

#endif // RECORD_H
