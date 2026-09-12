#ifndef FORWARD_H
#define FORWARD_H

#include <string>
#include <vector>

struct ForwardResult
{
    std::vector<double> CC;
    std::vector<double> DD;
};

ForwardResult forwardCosts(
    const std::string& A,
    const std::string& B,
    double g,
    double h,
    double tb
);

#endif
