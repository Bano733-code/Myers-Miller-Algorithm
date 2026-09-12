#include "midpoint.h"

#include <algorithm>
#include <limits>

using namespace std;

MidpointResult findMidpoint(
    const vector<double>& CC,
    const vector<double>& DD,
    const vector<double>& RR,
    const vector<double>& SS,
    double g)
{
    int N = static_cast<int>(CC.size()) - 1;

    MidpointResult result;

    result.j = 0;
    result.type = 1;
    result.cost = numeric_limits<double>::infinity();

    result.type1Values.resize(N + 1);
    result.type2Values.resize(N + 1);

    for (int j = 0; j <= N; ++j)
    {
        // Type 1 midpoint
        double type1 =
            CC[j] + RR[N - j];

        // Type 2 midpoint
        double type2 =
            DD[j] + SS[N - j] - g;

        result.type1Values[j] = type1;
        result.type2Values[j] = type2;

        if (type1 < result.cost)
        {
            result.cost = type1;
            result.j = j;
            result.type = 1;
        }

        if (type2 < result.cost)
        {
            result.cost = type2;
            result.j = j;
            result.type = 2;
        }
    }

    return result;
}
