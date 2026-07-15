#include <iostream>
#include <string>
#include <fstream>

#include "B_Spline.h"

using namespace std;

void main()
{
    int div = 100;

    B_Spline B;

    B.saveData();
    for (int h = 1; h <= div; h++)
    {
        B.define_d(h);
        B.define_u();
        B.define_zeta();
        B.de_boor_algorithm();
        B.savePOC(h);
    }
    B.showPOC();
    B.writePS();

}
