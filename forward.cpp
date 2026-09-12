#include "forward.h"
#include "scoring.h"

#include <algorithm>
#include <limits>

using namespace std;

ForwardResult forwardCosts(
    const string& A,
    const string& B,
    double g,
    double h,
    double tb)
{
    int M = static_cast<int>(A.size());
    int N = static_cast<int>(B.size());

    const double INF = numeric_limits<double>::infinity();

    vector<double> CC(N + 1, INF);
    vector<double> DD(N + 1, INF);

    // Initial row
    CC[0] = 0.0;

    for (int j = 1; j <= N; ++j)
    {
        CC[j] = g + h * j;
        DD[j] = CC[j] + g;
    }

    // Remaining rows
    for (int i = 1; i <= M; ++i)
    {
        double s = CC[0];

        double c = tb + h * i;
        CC[0] = c;

        double e = INF;

        for (int j = 1; j <= N; ++j)
        {
            // Insert state
            e = min(e, c + g) + h;

            // Delete state
            DD[j] = min(DD[j], CC[j] + g) + h;

            // Overall C state
            c = min({
                DD[j],
                e,
                s + substitutionCost(A[i - 1], B[j - 1])
            });

            s = CC[j];
            CC[j] = c;
        }
    }

    ForwardResult result;
    result.CC = CC;
    result.DD = DD;

    return result;
}
