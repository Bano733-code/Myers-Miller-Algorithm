#include "scoring.h"

double substitutionCost(char a, char b)
{
    return (a == b) ? 0.0 : 1.0;
}

double gapCost(int length, double g, double h)
{
    if (length <= 0)
        return 0.0;

    return g + h * length;
}
