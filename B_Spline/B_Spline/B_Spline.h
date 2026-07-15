#pragma once
#include <fstream>
#include "Point2D.h"

using namespace std;

class B_Spline
{
private:
   ifstream readFile;
   ofstream psFile;
   int div;

public:
    /*--- constructor ---*/
    B_Spline();
    /*--- destructor ---*/
    ~B_Spline();

    /*--- Variable ---*/
    int Degree = 0;
    int Number_of_CP = 0;
	int Number_of_Knots = 0;
    int R = 0;
    double parameter = 0;

    double *CP = nullptr;
	double *Knots = nullptr;

    double **d = nullptr;
    double **u = nullptr;
    double** zeta = nullptr;

    Point2D* POC;

    /*--- Function ---*/
    void saveData();
    void define_d(int h);
    void define_u();
    void define_zeta();
    void de_boor_algorithm();
    void savePOC(int h);
    void showPOC();
    void writePS();
};
