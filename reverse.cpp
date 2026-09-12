#include "reverse.h"
#include "forward.h"

#include <algorithm>

using namespace std;

ReverseResult reverseCosts(
    const string& A,
    const string& B,
    double g,
    double h,
    double te)
{
    string revA = A;
    string revB = B;

    reverse(revA.begin(), revA.end());
    reverse(revB.begin(), revB.end());

    ForwardResult temp =
        forwardCosts(revA, revB, g, h, te);

    ReverseResult result;

    result.RR = temp.CC;
    result.SS = temp.DD;

    return result;
}
