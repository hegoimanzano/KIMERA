#include "RandomGenerator.h"

RandomGenerator::RandomGenerator() {
    std::random_device rd;

    if (rd.entropy() != 0) {
        mt.seed(rd());
    }
    else {
        auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        mt.seed(seed);
    }
}

std::mt19937 & RandomGenerator::get() {
    return mt;
}

