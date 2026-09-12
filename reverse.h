#ifndef REVERSE_H
#define REVERSE_H

#include <string>
#include <vector>

struct ReverseResult
{
    std::vector<double> RR;
    std::vector<double> SS;
};

ReverseResult reverseCosts(
    const std::string& A,
    const std::string& B,
    double g,
    double h,
    double te
);

#endif
