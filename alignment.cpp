#include "alignment.h"

#include "forward.h"
#include "reverse.h"
#include "midpoint.h"
#include "scoring.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;


/*
    Small full-DP alignment.

    This is used only for small terminal subproblems.
    The main divide-and-conquer part uses linear-space
    forward/reverse calculations.
*/
AlignmentResult smallAlignment(
    const string& A,
    const string& B,
    double g,
    double h)
{
    int M = static_cast<int>(A.size());
    int N = static_cast<int>(B.size());

    const double INF = numeric_limits<double>::infinity();

    vector<vector<double>> C(
        M + 1,
        vector<double>(N + 1, INF)
    );

    vector<vector<double>> D(
        M + 1,
        vector<double>(N + 1, INF)
    );

    vector<vector<double>> I(
        M + 1,
        vector<double>(N + 1, INF)
    );

    C[0][0] = 0.0;

    // First row
    for (int j = 1; j <= N; ++j)
    {
        C[0][j] = g + h * j;
        I[0][j] = C[0][j];
    }

    // First column
    for (int i = 1; i <= M; ++i)
    {
        C[i][0] = g + h * i;
        D[i][0] = C[i][0];
    }

    // DP
    for (int i = 1; i <= M; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            D[i][j] =
                min(
                    D[i - 1][j],
                    C[i - 1][j] + g
                ) + h;

            I[i][j] =
                min(
                    I[i][j - 1],
                    C[i][j - 1] + g
                ) + h;

            C[i][j] =
                min({
                    D[i][j],
                    I[i][j],
                    C[i - 1][j - 1] +
                        substitutionCost(
                            A[i - 1],
                            B[j - 1]
                        )
                });
        }
    }

    // Backtracking
    string alignedA;
    string alignedB;

    int i = M;
    int j = N;

    while (i > 0 || j > 0)
    {
        if (i > 0 && j > 0)
        {
            double substitution =
                C[i - 1][j - 1] +
                substitutionCost(
                    A[i - 1],
                    B[j - 1]
                );

            if (abs(C[i][j] - substitution) < 1e-9)
            {
                alignedA += A[i - 1];
                alignedB += B[j - 1];

                --i;
                --j;
                continue;
            }
        }

        if (i > 0 &&
            abs(C[i][j] - D[i][j]) < 1e-9)
        {
            int oldI = i;

            while (i > 0)
            {
                if (i == oldI &&
                    abs(
                        D[i][j] -
                        (C[i - 1][j] + g + h)
                    ) < 1e-9)
                {
                    alignedA += A[i - 1];
                    alignedB += '-';

                    --i;
                    break;
                }

                if (i > 0 &&
                    abs(
                        D[i][j] -
                        (D[i - 1][j] + h)
                    ) < 1e-9)
                {
                    alignedA += A[i - 1];
                    alignedB += '-';

                    --i;
                }
                else
                {
                    break;
                }
            }

            continue;
        }

        if (j > 0)
        {
            alignedA += '-';
            alignedB += B[j - 1];
            --j;
        }
    }

    reverse(alignedA.begin(), alignedA.end());
    reverse(alignedB.begin(), alignedB.end());

    AlignmentResult result;

    result.alignedA = alignedA;
    result.alignedB = alignedB;
    result.cost = C[M][N];

    return result;
}


/*
    Print the forward vectors.
*/
void printForward(
    const ForwardResult& F)
{
    cout << "\n========== FORWARD ==========\n";

    cout << "CC: ";

    for (double x : F.CC)
        cout << fixed << setprecision(1) << x << " ";

    cout << "\n";

    cout << "DD: ";

    for (double x : F.DD)
        cout << fixed << setprecision(1) << x << " ";

    cout << "\n";
}


/*
    Print the reverse vectors.
*/
void printReverse(
    const ReverseResult& R)
{
    cout << "\n========== REVERSE ==========\n";

    cout << "RR: ";

    for (double x : R.RR)
        cout << fixed << setprecision(1) << x << " ";

    cout << "\n";

    cout << "SS: ";

    for (double x : R.SS)
        cout << fixed << setprecision(1) << x << " ";

    cout << "\n";
}


/*
    Print Type 1 / Type 2 midpoint calculations.
*/
void printMidpointTable(
    const MidpointResult& M)
{
    cout << "\n========== MIDPOINT ==========\n";

    cout << left
         << setw(8) << "j"
         << setw(14) << "Type 1"
         << setw(14) << "Type 2"
         << "\n";

    cout << "--------------------------------\n";

    for (size_t j = 0;
         j < M.type1Values.size();
         ++j)
    {
        cout << left
             << setw(8) << j
             << setw(14)
             << fixed << setprecision(1)
             << M.type1Values[j]
             << setw(14)
             << M.type2Values[j]
             << "\n";
    }

    cout << "\nChosen midpoint:\n";
}


/*
    Recursive Myers-Miller-style alignment.
*/
AlignmentResult recursiveAlign(
    const string& A,
    const string& B,
    double g,
    double h,
    bool trace,
    int depth)
{
    int M = static_cast<int>(A.size());
    int N = static_cast<int>(B.size());

    // Empty A
    if (M == 0)
    {
        AlignmentResult result;

        result.alignedA = string(N, '-');
        result.alignedB = B;
        result.cost = gapCost(N, g, h);

        return result;
    }

    // Empty B
    if (N == 0)
    {
        AlignmentResult result;

        result.alignedA = A;
        result.alignedB = string(M, '-');
        result.cost = gapCost(M, g, h);

        return result;
    }

    /*
        Small terminal problem.
    */
    if (M <= 2 || N <= 2)
    {
        return smallAlignment(A, B, g, h);
    }

    /*
        Divide A at its midpoint.
    */
    int mid = M / 2;

    string Aleft =
        A.substr(0, mid);

    string Aright =
        A.substr(mid);

    /*
        Forward calculation.
    */
    ForwardResult F =
        forwardCosts(
            Aleft,
            B,
            g,
            h,
            0.0
        );

    /*
        Reverse calculation.
    */
    ReverseResult R =
        reverseCosts(
            Aright,
            B,
            g,
            h,
            0.0
        );

    /*
        Only print the top-level trace.

        Otherwise recursive calls would produce a huge
        amount of output.
    */
    bool printThisLevel = trace && depth == 0;

    if (printThisLevel)
    {
        printForward(F);
        printReverse(R);
    }

    /*
        Find midpoint.
    */
    MidpointResult Mpoint =
        findMidpoint(
            F.CC,
            F.DD,
            R.RR,
            R.SS,
            g
        );

    if (printThisLevel)
    {
        printMidpointTable(Mpoint);

        cout << "i* = "
             << mid
             << "\n";

        cout << "j* = "
             << Mpoint.j
             << "\n";

        cout << "Type = "
             << Mpoint.type
             << "\n";

        cout << "Midpoint cost = "
             << fixed << setprecision(1)
             << Mpoint.cost
             << "\n";
    }

    /*
        Split B at j*.
    */
    int splitJ = Mpoint.j;

    string Bleft =
        B.substr(0, splitJ);

    string Bright =
        B.substr(splitJ);

    /*
        Recursively solve left half.
    */
    AlignmentResult leftResult =
        recursiveAlign(
            Aleft,
            Bleft,
            g,
            h,
            trace,
            depth + 1
        );

    /*
        Recursively solve right half.
    */
    AlignmentResult rightResult =
        recursiveAlign(
            Aright,
            Bright,
            g,
            h,
            trace,
            depth + 1
        );

    /*
        Combine both alignments.
    */
    AlignmentResult result;

    result.alignedA =
        leftResult.alignedA +
        rightResult.alignedA;

    result.alignedB =
        leftResult.alignedB +
        rightResult.alignedB;

    /*
        Recalculate the cost from the final alignment.
    */
    double totalCost = 0.0;

    int k = 0;

    while (k < static_cast<int>(result.alignedA.size()))
    {
        char a = result.alignedA[k];
        char b = result.alignedB[k];

        /*
            Gap in B = deletion from A.
        */
        if (b == '-')
        {
            int length = 0;

            while (
                k + length <
                    static_cast<int>(
                        result.alignedA.size()
                    )
                &&
                result.alignedB[k + length] == '-'
            )
            {
                ++length;
            }

            totalCost +=
                gapCost(length, g, h);

            k += length;
        }

        /*
            Gap in A = insertion into A.
        */
        else if (a == '-')
        {
            int length = 0;

            while (
                k + length <
                    static_cast<int>(
                        result.alignedA.size()
                    )
                &&
                result.alignedA[k + length] == '-'
            )
            {
                ++length;
            }

            totalCost +=
                gapCost(length, g, h);

            k += length;
        }

        /*
            Match / mismatch.
        */
        else
        {
            totalCost +=
                substitutionCost(a, b);

            ++k;
        }
    }

    result.cost = totalCost;

    return result;
}


/*
    Public function.
*/
AlignmentResult alignMyersMiller(
    const string& A,
    const string& B,
    double g,
    double h,
    bool trace)
{
    AlignmentResult result =
        recursiveAlign(
            A,
            B,
            g,
            h,
            trace,
            0
        );

    if (trace)
    {
        cout << "\n========== FINAL ==========\n";

        cout << "Optimal alignment:\n";

        cout << result.alignedA << "\n";
        cout << result.alignedB << "\n";

        cout << "\nMinimum cost = "
             << fixed << setprecision(1)
             << result.cost
             << "\n";
    }

    return result;
}
