#pragma once
double Tetta(double sigma);
struct Coord_4
{
	double x, y, z;
	Coord_4(double x, double y, double z);
};
double Get_Index_4(double x, double y, double z, int ROW, int COL, int LIST, double L_side);
Coord_4 Get_Coord_4(int R, int C, int L, int Row, int Col, int List, double L_side);
double F_Tettta_AB(double x, double y, double z, double t);
double F_Tettta_BC(double x, double y, double z, double t);
double F_4(double x, double y, double z, double t);
double Fi_4(double x, double y, double z, double t);
double Analit_Task4(double x, double y, double z, double t, double type);
void Build_Tesk4(int Row, int Col, int List, double& Lenght, double* T, double* Num_obl, int* check_oblast, double* ro, int* type_bondary,
    double (**arrF)(double, double, double, double), double(**arrFi)(double, double, double, double));

