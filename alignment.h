#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>

struct AlignmentResult
{
    std::string alignedA;
    std::string alignedB;
    double cost;
};

AlignmentResult alignMyersMiller(
    const std::string& A,
    const std::string& B,
    double g,
    double h,
    bool trace
);

#endif
