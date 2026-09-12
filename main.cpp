#include <iostream>
#include <iomanip>
#include <string>

#include "alignment.h"

using namespace std;

int main()
{
    string A, B;

    cout << "========================================\n";
    cout << " Myers-Miller Linear-Space Alignment\n";
    cout << "========================================\n\n";

    cout << "Enter sequence A: ";
    cin >> A;

    cout << "Enter sequence B: ";
    cin >> B;

    // Scoring parameters from the paper-style example
    double g = 2.0;
    double h = 0.5;

    cout << "\nScoring parameters:\n";
    cout << "Match cost    = 0\n";
    cout << "Mismatch cost = 1\n";
    cout << "Gap opening   = " << g << '\n';
    cout << "Gap extension = " << h << '\n';

    /*
        true = show algorithm trace
        false = only return final result
    */
    AlignmentResult result =
        alignMyersMiller(
            A,
            B,
            g,
            h,
            true
        );

    return 0;
}
