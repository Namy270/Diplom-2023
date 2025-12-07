#include "Tesk4.h"
#include <fstream>
#include <cmath>
#include <iostream>
using namespace std;

Coord_4 Get_Coord_4(int R, int C, int L, int Row, int Col, int List, double L_side)
{
	double x, y, z;
	double x_length = (List - 2) * L_side;
	z = L_side * R - L_side / 2;
	y = L_side * C - L_side / 2;
	x = x_length -  (L_side * L - L_side / 2);
	return Coord_4(x, y, z);

}
double Get_Index_4(double x, double y, double z, int ROW, int COL, int LIST, double L_side)
{
	int inner_cell = ROW * COL + ROW + 1;
	double X_index = ((LIST - 2) * L_side - x) / L_side;
	double Y_index = y / L_side;
	double Z_index = z / L_side;
	int Index = ROW * COL * int(ceil(X_index)) + ROW * int(ceil(Y_index)) + int(ceil(Z_index));
	return Index;
}
double Tetta(double sigma)
{
	double pog = 1, x0, x2, x3, Del;
	if (sigma <= 1e-12)
		return 0;
	else
	{
		x0 = sqrt(6 * sigma);
		for (int i = 1; i < 1000000; i++)
		{
			Del = (3 * sigma - x0 + log(1 + x0)) / x0;
			pog = max(abs(Del), abs(Del * x0));
			if (pog <= 1e-10)
				return pow(x0, 0.3333333333333);
			x2 = x0 + (1 + x0) * Del;
			x3 = max(x2, 0.5 * x0);
			x0 = min(x3, 2 * x0);
		}
		return pow(x0, 0.3333333333333);
	}
}
double Analit_Task4(double x, double y, double z, double t, double type)
{
	double ai, bi, ni, si, ui, s0i, ti;
	/*x -= 0.2;
	y += 98;*/
	if (type == 1)
	{
		ti = 0;
		s0i = 141.367176507;
		ui = 10.0;
		ai = -0.995037190210;
		bi = 0.0995037190210;
		ni = -9.95037190210;
		si = 0.995037190210;
	}
	else if (type == 2)
	{
		ti = -3.98014876084;
		s0i = 127.073287136;
		ui = 8.98881202170;
		ai = -0.894427191000;
		bi = 0.447213595500;
		ni = -1.99007438042;
		si = 0.995037190210;
	}
	else if (type == 3)
	{
		//ti = -501.448992006;
		//s0i = 2.84087696303;
		//ui = 0.200957324976;
		//ai = -0.0199960011996;
		//bi = 0.999100059980;
		//ni = -0.0199007438042;
		//si = 0.995037190210;
		ti = -501.449;
		s0i = 2.84087;
		ui = 0.200957;
		ai = -0.019996;
		bi = 0.9998;
	}
	//else
	//	std::cout << "Error\n";
	return 2.15443769003 * Tetta(-(ai * (x - 0.2) + bi * (y + 98) - ui * (t - ti)) / s0i);
}

double F_Tettta_AB(double x, double y, double z, double t)
{
	double sigma = t / 14.1367176507;

	return Tetta(sigma) * pow(10., 1. / 3);
}
double F_Tettta_BC(double x, double y, double z, double t)
{
	double sigma = 0.0707377783663 * t - 0.00351933601136 * (y + 98) + 0.351933601136;

	return Tetta(sigma) * pow(10., 1. / 3);
}
double F_4(double x, double y, double z, double t)
{
	return 0.0;
}
double Fi_4(double x, double y, double z, double t)
{
	return 0.0;
}
void Build_Tesk4(int Row, int Col, int List, double& Length, double* T, double* Num_obl, int* check_oblast, double* ro, int* type_bondary,
    double (**arrF)(double, double, double, double), double(**arrFi)(double, double, double, double))
{
	for (int i = 0; i < Row; i++) // по z
	{
		for (int j = 0; j < Col; j++)  // по у от 0 до lenght * Col
		{
			for (int k = 0; k < List; k++) // по x от lenght* List до 0
			{
				int index = i + j * Row + k * Row * Col;
				Coord_4 xyz = Get_Coord_4(i, j, k, Row, Col, List, Length);
				T[index] = 0;
				type_bondary[index] = 0;
				check_oblast[index] = 1;
				arrF[index] = nullptr;
				arrFi[index] = nullptr;
				if (i == 0 || i == Row - 1)
				{
					type_bondary[index] = 2;
					check_oblast[index] = -1;
					Num_obl[index] = 0;
					arrF[index] = F_4;
					arrFi[index] = Fi_4;

				}
				else
				{
					if (j * Length < 2) // Заполняем 1 область
					{
						Num_obl[index] = 1;
						ro[index] = 1.;

						if (j * Length <= 1 && k * Length > 0.2) // Рассматриваем нижнюю часть
						{
							double aA_y = 1 - xyz.x / 10.;
							if (xyz.y > aA_y) // внутренняя часть
							{
								//T[index] = 3;
							}
							else
							{
								if (aA_y - xyz.y > Length) // внешняя часть
								{
									check_oblast[index] = 0;
									//T[index] = 6;
								}
								else // граница
								{
									type_bondary[index] = 2;
									check_oblast[index] = -1;
									arrF[index] = nullptr;
									arrFi[index] = nullptr;

									//T[index] = 0;
								}
							}
						}
						else
						{
							if (k == List - 1) // Граница А'D'
							{
								type_bondary[index] = 2;
								check_oblast[index] = -1;
								Num_obl[index] = -1;
								arrF[index] = nullptr;
								arrFi[index] = nullptr;

								//T[index] = 0;
							}
							else
							{
								// y = 10x - 100 | x = 10 + y/10
								double AB_x = xyz.y / 10. + 10;
								if (xyz.x < AB_x)			// Внутренняя часть
								{
									//T[index] = 3;
								}
								else
								{
									if (xyz.x - AB_x > Length) // Внешняя часть
									{
										check_oblast[index] = 0;
										//T[index] = 6;
									}
									else						// Граница
									{
										type_bondary[index] = 1;
										check_oblast[index] = -1;
										arrF[index] = F_Tettta_AB;
										arrFi[index] = F_Tettta_AB;
										//T[index] = 1;
									}
								}
							}
						}
					}
					else if (j * Length <= 3)  // Заполняем вторую область
					{
						if (k == List - 1) // Граница А'D'
						{
							type_bondary[index] = 2;
							check_oblast[index] = -1;
							Num_obl[index] = -1;
							arrF[index] = nullptr;
							arrFi[index] = nullptr;
							//T[index] = 0;
						}
						else
						{
							Num_obl[index] = 2;
							ro[index] = 2.;
							if (k == 0) // если граница BC
							{
								type_bondary[index] = 1;
								check_oblast[index] = -1;
								arrF[index] = F_Tettta_BC;
								arrFi[index] = F_Tettta_BC;

								//T[index] = 2;
							}
							else
							{
								//T[index] = 4;
							}
						}
					}
					else  // Заполняем третью ообласть
					{
						if (k == List - 1) // Граница А'D'
						{
							type_bondary[index] = 2;
							check_oblast[index] = -1;
							arrF[index] = F_4;
							arrFi[index] = nullptr;

							//T[index] = 0;
						}
						else
						{
							Num_obl[index] = 3;
							ro[index] = 10.;
							if (Length * k < 0.006)
							{
								type_bondary[index] = 2;
								check_oblast[index] = -1;
								arrF[index] = F_4;
								arrFi[index] = nullptr;

								//T[index] = 0;

								//// y = 513 - 50x | x = 10.26 - y/50
								//double DC_x = 10.26 - xyz.y / 50.;
								//if (xyz.x < DC_x)
								//{
								//	type_bondary[index] = 0;
								//	check_oblast[index] = 1;
								//	arrF[index] = nullptr;
								//	arrFi[index] = nullptr;
								//	T[index] = 5;
								//}
								//else
								//{
								//	if (xyz.x - DC_x > Length)
								//	{
								//		type_bondary[index] = 0;
								//		check_oblast[index] = 0;
								//		arrF[index] = nullptr;
								//		arrFi[index] = nullptr;
								//		T[index] = 6;
								//	}
								//	else
								//	{
								//		type_bondary[index] = 1;
								//		check_oblast[index] = -1;
								//		arrF[index] = F_4;
								//		arrFi[index] = nullptr;
								//		T[index] = 0;
								//	}
								//}
							}
							else
							{
								if (j == Col - 1)
								{
									type_bondary[index] = 2;
									check_oblast[index] = -1;
									arrF[index] = F_4;
									arrFi[index] = nullptr;

									//T[index] = 0;
								}
								else
								{
									//T[index] = 5;
								}
							}
						}
					}
				}
				//ro[index] = 1.;
				//type_bondary[index] = 2;
				//Num_obl[index] = 0;
				//T[index] = 1;

			}
		}
	}
	//std::ofstream file2("Grafik 2d.txt");
	//if (!file2.is_open())
	//{
	//	std::cout << "Error open file test for Grafik 2d!\n";
	//	exit(0);
	//}
	//else
	//{
	//	for (int i = 0; i < List; i++)
	//	{
	//		for (int j = 0; j < Col; j++)
	//		{
	//			int index = 1 + j * Row + i * Row * Col;
	//			Coord_4 xyz = Get_Coord_4(1, j, i, Row, Col, List, Length);
	//			file2 << xyz.x << "\t" << xyz.y << "\t" << T[index] << "\n";
	//		}
	//		file2 << "\n";
	//	}

	//}
	//file2.close();
}

Coord_4::Coord_4(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}
