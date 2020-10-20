#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H


#include <iostream>
#include <random>
#include <chrono>

class RandomGenerator
{
public:
    static RandomGenerator& Instance() {
        static RandomGenerator s;
        return s;
    }
    std::mt19937 & get();

private:
    RandomGenerator();
    ~RandomGenerator() {}

    RandomGenerator(RandomGenerator const&) = delete;
    RandomGenerator& operator= (RandomGenerator const&) = delete;

    std::mt19937 mt;
};

#endif // RANDOMGENERATOR_H
