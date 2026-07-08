#include "B_Spline.h"
#include "Point2D.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

B_Spline::B_Spline()
{
    /*--- Input Data ---*/
    Degree = 3;
    Number_of_CP = Degree + 2;
    div = 100;

    POC = new Point2D[div + 1];
}

B_Spline::~B_Spline()
{

}

void B_Spline::saveData()
{
	/*--- CP.txt File Open ---*/
    readFile.open("CP.txt");
    if (!readFile)
    {
        cout << "This CP.txt file can't Open" << endl;
    }
    /*--- Dynamic allocation CP ---*/
    CP = new double[Number_of_CP];
    /*--- Save CP Data ---*/
    cout << "===== CP =====" << endl;
    for (int i = 0; i < Number_of_CP; i++)
    {
        readFile >> CP[i];
        cout << CP[i] << endl;
    }
    readFile.close();

	/*--- Knots.txt File Open ---*/
    readFile.open("Knots.txt");
    if (!readFile)
    {
        cout << "This Knots.txt file can't Open" << endl;
    }
	/*--- Save Number of Knots ---*/
	readFile >> Number_of_Knots;
	cout << "===== Number of Knots =====" << endl;
	cout << Number_of_Knots << endl;
	/*--- Dynamic allocation Knots ---*/
	Knots = new double[Number_of_Knots];
	cout << "===== Knots =====" << endl;
	for (int i = 0; i < Number_of_Knots; i++)
	{
		readFile >> Knots[i];
		cout << Knots[i] << endl;
	}
    readFile.close();
}

void B_Spline::define_d(int h)
{
    /*--- Dynamic allocation Matrix d ---*/
    d = new double*[Degree + 1];
    for (int k = 0; k <= Degree; k++)
    {
        d[k] = new double[Degree + Number_of_CP];
    }
    /*--- Save CP Data to Matrix d ---*/
    cout << "==== Matrix d at k = 0 ====" << endl;
    for (int i = 0; i < Number_of_CP; i++)
    {
        d[0][i] = CP[i];
        cout << d[0][i] << endl;
    }

    parameter = Knots[0] + (Knots[Number_of_Knots - 1] - Knots[0]) / div * h;
}

void B_Spline::define_u()
{
    /*--- Dynamic allocation Matrix u ---*/
    u = new double* [Degree + 1];
    for (int k = 0; k <= Degree; k++)
    {
        u[k] = new double[Degree + Number_of_Knots];
    }
    /*--- Save CP Data to Matrix u ---*/
    cout << "==== Matrix u at k = 0 ====" << endl;
    for (int i = 0; i < Number_of_Knots; i++)
    {
        u[0][i] = Knots[i];
        cout << u[0][i] << endl;
    }

    /*--- Add knots ---*/
    for (int k = 1; k <= Degree; k++)
    {
        for (int i = 0; i < Number_of_Knots; i++)
        {
            u[k][i] = Knots[i];
        }
        for (int i = Number_of_Knots; i < Number_of_Knots + k; i++)
        {
            u[k][i] = parameter;
        }

        /*--- 재정렬 ---*/
        sort(u[k] + 0, u[k] + (Number_of_Knots + k));

        cout << "==== k" << " = " << k << " ====" << endl;
        for (int i = 0; i < Number_of_Knots + k; i++)
        {
            cout << "u[" << k << "][" << i << "]" << " = " << u[k][i] << endl;
        }
    }
}

void B_Spline::define_zeta()
{
    zeta = new double* [Degree + 1];
    for (int k = 0; k <= Number_of_CP + Degree; k++)
    {
        zeta[k] = new double[Number_of_CP + Degree];
    }

    /*--- Moving Average ---*/
    cout << "==== zeta[1][i] ====" << endl;
    for (int i = 0; i < Number_of_CP; i++)
    {
        zeta[1][i] = (Knots[i] + Knots[i + 1] + Knots[i + 2]) / Degree;
        cout << "zeta[1][" << i << "]" << " " << "=" << " " << zeta[1][i] << endl;
    }

    for (int k = 2; k <= Degree; k++)
    {
        cout << "==== zeta[" << k <<"][i] ====" << endl;
        for (int i = 0; i <= Number_of_CP + k - 2; i++)
        {
            zeta[k][i] = (u[k - 1][i] + u[k - 1][i + 1] + u[k - 1][i + 2]) / Degree;
            cout << "zeta[" << k << "][" << i << "]" << " " << " = " << " " << zeta[k][i] << endl;
        }
    }
}

void B_Spline::de_boor_algorithm()
{
    for (int k = 1; k <= Degree; k++)
    {
        int m_old = Number_of_Knots + (k - 1);   // u[k-1]의 길이
        int n_new = Number_of_CP + k;            // 이번 단계 제어점 개수

        /*--- parameter가 속한 knot 구간 j 찾기 ---*/
        int j = 0;
        for (int s = 0; s < m_old - 1; s++)
        {
            if (u[k - 1][s] <= parameter && u[k - 1][s] < u[k - 1][m_old - 1])
                j = s;
        }

        for (int i = 0; i < n_new; i++)
        {
            if (i <= j - 2)                      // 앞쪽: 영향 없음
                d[k][i] = d[k - 1][i];
            else if (i >= j + 2)                 // 뒤쪽: 한 칸 밀림
                d[k][i] = d[k - 1][i - 1];
            else                                 // i = j-1, j, j+1 만 보간
            {
                double denom = u[k - 1][i + 2] - u[k - 1][i - 1];
                if (denom == 0) { d[k][i] = d[k - 1][i - 1]; continue; }
                double alpha = (parameter - u[k - 1][i - 1]) / denom;
                d[k][i] = (1 - alpha) * d[k - 1][i - 1] + alpha * d[k - 1][i];
            }
        }
    }
}

void B_Spline::savePOC(int h)
{
    int istar = 0;
    for (int i = 0; i < Number_of_Knots + Degree; i++)
        if (u[Degree][i] < parameter) istar++;

    POC[h].x = parameter;
    POC[h].y = d[Degree][istar];

    cout << POC[h].x << " " << POC[h].y << endl;
}

void B_Spline::showPOC()
{
    cout << "==== All POC ====" << endl;
    for (int i = 1; i <= div; i++)
    {
        cout << POC[i].x << " " << POC[i].y << endl;
    }
}
