#ifndef MIDPOINT_H
#define MIDPOINT_H

#include <vector>

struct MidpointResult
{
    int j;
    int type;
    double cost;

    // For displaying the complete midpoint table
    std::vector<double> type1Values;
    std::vector<double> type2Values;
};

MidpointResult findMidpoint(
    const std::vector<double>& CC,
    const std::vector<double>& DD,
    const std::vector<double>& RR,
    const std::vector<double>& SS,
    double g
);

#endif
