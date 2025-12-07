#include <iostream>
#include "Grid.h"
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <string>
#include <functional>
#include <chrono>
#include "Tesk4.h"
using namespace std;

void SetColor(int text, int background)
{
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, (WORD)((background << 4) | text));
}

double F_(double x, double y, double z, double t)
{
	return 0.0;
}
double Fi_(double x, double y, double z, double t)
{
	return 0.00;
}
double Analit_1(double x, double y, double z, double t, double type)
{
	double pi = 3.1415926535e0;
	double sumX{}, sumY{}, sumZ{};
	for (int k = 0; k < 3000; k++)
	{
		sumX += exp(-std::pow((2.0 * k + 1) * pi, 2) * t) / (pi * (2.0 * k + 1)) * sin(pi * (2.0 * k + 1) * x);
		sumY += exp(-std::pow((2.0 * k + 1) * pi, 2) * t) / (pi * (2.0 * k + 1)) * sin(pi * (2.0 * k + 1) * y);
		sumZ += exp(-std::pow((2.0 * k + 1) * pi, 2) * t) / (pi * (2.0 * k + 1)) * sin(pi * (2.0 * k + 1) * z);
	}
	return sumX * sumY * sumZ * 64;
}
int main(int argc, const char** argv)
{
	setlocale(0, "");
	double Length, L_side;
	int Task;
	int sizeGrid, Row, Col, List, iters;
	double time_step;
	double Time;
	double Gradient;
	string solv;
	auto begin_time = std::chrono::steady_clock::now();
	auto end_time = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time);
	// Чтение данных
	std::ifstream file_comand("comand.txt");
	if (!file_comand.is_open())
	{
		std::cout << "Error open file comand.txt" << std::endl;
		exit(0);
	}
	else
	{
		std::string comand_s;
		file_comand >> comand_s >> Task;
		file_comand >> comand_s >> solv >> iters;
		file_comand >> comand_s >> sizeGrid;
		file_comand >> comand_s >> time_step;
		file_comand >> comand_s >> Time;
		file_comand >> comand_s >> Gradient;
	}
	file_comand.close();


	// Заполняем данные
	if (Task == 1)
	{
		Length = 1;
		Col = List = Row = sizeGrid;
		L_side = Length / (sizeGrid - 2);

	}
	else if (Task == 4)
	{

		Row = 3, Col = 33 * sizeGrid + 2, List = 102 * sizeGrid + 2;
		Length = 10.2;
		L_side = 3.3 / (Col - 2);
	}
	int size = Row * Col * List;
	double* T, * dE, * Number_obl, * ro;
	int* type_bondary;
	typedef double (*funcF)(double, double, double, double);
	funcF* arrF, * arrFi;
	double (*Analit)(double, double, double, double, double);
	int* check_oblast;
	T = new double[size];
	ro = new double[size];
	dE = new double[size];
	Number_obl = new double[size];
	check_oblast = new int[size];
	arrF = new funcF[size];
	arrFi = new funcF[size];
	type_bondary = new int[size];

	////////////////////////////////////////////////////////////////////////////////////////
	if (Task == 4)
	{
		Build_Tesk4(Row, Col, List, L_side, T, Number_obl, check_oblast, ro, type_bondary, arrF, arrFi);
		Analit = Analit_Task4;
		/*int Finish_iter = round(Time / time_step);
		for (int k = Finish_iter; k <= Finish_iter; k++)
		{
			//Sleep(5);

			std::ofstream file2("Grafik 2d.txt");
			if (!file2.is_open())
			{
				std::cout << "Error open file test for Grafik 2d!\n";
				exit(0);
			}
			else
			{
				double center_z = ((Row - 2) / 2.0) * L_side;
				for (int i = 0; i < List; i++)
				{
					for (int j = 0; j < Col; j++)
					{
						int index = 1 + j * Row + i * Row * Col;
						Coord_4 xyz = Get_Coord_4(1, j, i,Row, Col, List, L_side);
						if(check_oblast[index] != 0)
							file2 << xyz.x << "\t" << xyz.y << "\t" << Analit_Task4(xyz.x, xyz.y, center_z, k * time_step, Number_obl[index]) << "\n";
						else
							file2 << xyz.x << "\t" << xyz.y << "\t" << 0 << "\n";

					}
					file2 << "\n";
				}

			}
			file2.close();
		}

		std::ofstream file11("T = 0.5, X = 6, X = 8, X = 10.txt");
		if (!file11.is_open())
		{
			printf("Error open file11!");
			exit(1);
		}
		else
		{
			file11 << "\t  Y  \t|\t  X = 5.8  \t|\t  X = 7.8  \t|\tX = 9.8\n";
			file11 << "______________________________________________________________________________________\n";
			for (int iii = Col - 1; iii >= 0; iii--)
			{
				double center_z = ((Row - 2) / 2.0) * L_side;
				double center_x = ((List - 2) / 2.0) * L_side;

				Coord_4 Y = Get_Coord_4(center_z, iii, center_x, Row, Col, List, L_side);
				int Ind = Get_Index_4(center_x, Y.y, center_z, Row, Col, List, L_side);
				file11 << std::fixed << std::setprecision(6) << Y.y + 98 << "\t|\t" << Analit_Task4(6, Y.y, Y.z, Time, Number_obl[Ind])
					<< "\t|\t" << Analit_Task4(8, Y.y, Y.z, Time, Number_obl[Ind]) << "\t|\t" << Analit_Task4(10, Y.y, Y.z, Time, Number_obl[Ind]) << "\n";

			}
		}
		file11.close();
		return 0;*/
	}
	else if (Task == 1)
	{
		Analit = Analit_1;
		for (int i = 0; i < Row; i++)
		{
			for (int j = 0; j < Col; j++)
			{
				for (int k = 0; k < List; k++)
				{
					int index = i + j * Row + k * Row * Col;
					check_oblast[index] = -1;
					T[index] = 0;
					if (j == 0) // Левая сторона
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else if (j == Col - 1) // Правая сторона
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else if (i == 0) // Верхняя сторона 
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else if (i == Row - 1) // Нижняя сторона
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else if (k == 0) // Передняя сторона
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else if (k == List - 1) // Задняя сторона
					{
						type_bondary[index] = 1;
						arrF[index] = F_;
					}
					else
					{
						type_bondary[index] = 0;
						T[index] = 1;
						check_oblast[index] = 1;
						arrF[index] = nullptr;
						arrFi[index] = nullptr;
					}
					/*					if ((k == 1 && j == 1 )||( k == 1 && j == 1))
										{
											type_bondary[index] = 1;
											check_oblast[index] = -1;
											arrF[index] = Fi_;
										}*/
					Number_obl[index] = 0;
					ro[index] = 1;
				}
			}
		}
	}
	Grid Test1(Row, Col, List);
	Test1.init(T, Number_obl, check_oblast, ro, L_side, type_bondary, arrF, arrFi, Gradient);

	double An, Im, sumPogr{}, sumPogrOtn{};
	begin_time = std::chrono::steady_clock::now();
	Test1.Calculate(time_step, Time, solv, iters); // Решаем
	end_time = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time);
	int time_n = elapsed_ms.count() / 1000;
	int hours = time_n / 3600;
	int minuts = (time_n % 3600) / 60;
	int seconds = (time_n % 3600) % 60;

	// Выводим результаты
	ofstream file("test1.txt");
	if (!file.is_open())
	{
		printf("Error open file1!");
		exit(1);
	}
	else
	{
		int Index;
		double center_z = ((Row - 2) / 2.0) * L_side;
		double center_x = ((List - 2) / 2.0) * L_side;
		double center_y = ((Col - 2) / 2.0) * L_side;
		double Linde_rez_y, Linde_rez_x;
		std::cout << "\n";
		std::cout << "Time Build adapt = " << std::fixed << std::setprecision(3) << Test1.Time_program.adapt_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.adapt_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "Time memory SLAY = " << std::fixed << std::setprecision(3) << Test1.Time_program.memory_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.memory_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "Time build SLAY  = " << std::fixed << std::setprecision(3) << Test1.Time_program.build_slay_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.build_slay_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "Time build CSR   = " << std::fixed << std::setprecision(3) << Test1.Time_program.CSR_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.CSR_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "Time solve SLAY  = " << std::fixed << std::setprecision(3) << Test1.Time_program.Slay_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.Slay_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "Time recording T = " << std::fixed << std::setprecision(3) << Test1.Time_program.result_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.result_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time Build adapt = " << std::fixed << std::setprecision(3) << Test1.Time_program.adapt_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.adapt_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time memory SLAY = " << std::fixed << std::setprecision(3) << Test1.Time_program.memory_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.memory_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time build SLAY  = " << std::fixed << std::setprecision(3) << Test1.Time_program.build_slay_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.build_slay_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time build CSR   = " << std::fixed << std::setprecision(3) << Test1.Time_program.CSR_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.CSR_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time solve SLAY  = " << std::fixed << std::setprecision(3) << Test1.Time_program.Slay_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.Slay_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		file << "Time recording T = " << std::fixed << std::setprecision(3) << Test1.Time_program.result_time << "\t|\t" << std::setprecision(2) << Test1.Time_program.result_time / Test1.Time_program.sum_time * 100 << "%" << std::endl;
		std::cout << "\nTimer programm: " << hours << " h : " << minuts << " min : " << seconds << " sec\n\n";
		std::cout << "Coordinate (" << center_x << ", y, " << center_z << ") |\tImplicit\t|\tAnalit  \t|\tAbs. pogr\t|\tRel. pogr (%)\n";
		std::cout << "____________________________________________________________________________________________________________________" << std::endl;
		file << "\nTimer programm: " << hours << " h : " << minuts << " min : " << seconds << " sec\n\n";
		file << "Coordinate   \t|\tImplicit\t|\tAnalit  \t|\tAbs. pogr\t|\tRel. pogr (%)\n";
		file << "____________________________________________________________________________________________________________________" << std::endl;

		for (int i = 0; i < List; i++)
		{
			//Linde_rez_y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1);
			Linde_rez_x = Length / (List - 2) / 2.0 + Length / (List - 2.0) * (i - 1);
			Index = Get_Index_4(Linde_rez_x, center_y, center_z, Row, Col, List, L_side);
			if (i == 0 || i == Col - 1)
			{
				std::cout << fixed << setprecision(8) << "\t" << 0.0 << "\t|\t";
				//SetColor(14, 0);
				std::cout << 0.0 << "\t|\t";
				//SetColor(10, 0);
				std::cout << 0.0 << "\t|\t";
				//SetColor(15, 0);
				std::cout << 0.0 << "\t|\t" << 0.0 << "\n";

				file << fixed << setprecision(8) << "\t" << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t" << 0.0 << "\n";
			}
			else
			{
				std::cout << fixed << setprecision(8) << "\t" << Linde_rez_x << "\t|\t";
				//SetColor(14, 0);
				Im = Test1.Get_adapt_T(Linde_rez_x, center_y, center_z);
				std::cout << Im << "\t|\t";
				//SetColor(10, 0);
				An = Analit(Linde_rez_x, center_y, center_z, Time, Number_obl[Index]);
				std::cout << An << "\t|\t";
				//SetColor(15, 0);
				std::cout << abs(An - Im) << "\t|\t" << abs((An - Im) / Im) * 100 << "\n";

				file << fixed << setprecision(8) << "\t" << Linde_rez_x << "\t|\t";
				file << Im << "\t|\t";
				file << An << "\t|\t";
				file << abs(An - Im) << "\t|\t" << abs((An - Im) / Im) * 100 << "\n";

				sumPogr += abs(An - Im);
				sumPogrOtn += abs((An - Im) / Im);
			}
			//for(int j = 0; j < Row; j++)
			//	for (int k = 0; k < List; k++)
			//	{
			//		sumPogr += abs(Analit(Length / (List - 2) / 2.0 + Length / (List - 2.0) * (k - 1),
			//			Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1),
			//			Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (j - 1), 0.02) -
			//			Test1.GetT(Length / (List - 2) / 2.0 + Length / (List - 2.0) * (k - 1),
			//				Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1),
			//				Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (j - 1)));
			//	}
		}
		std::cout << endl;
		std::cout << "\nСредняя абсолютная погрешность = " << sumPogr / Col;
		std::cout << "\nСредняя оносительная погрешность = " << sumPogrOtn / Col * 100 <<" %" << endl << endl;
		file << endl;
		file << "\nСредняя абсолютная погрешность = " << sumPogr / Col;
		file << "\nСредняя оносительная погрешность = " << sumPogrOtn / Col * 100 << " %" << endl << endl;

		//////////////////////////////////////////////////////////////
		// Сравниваем на адаптивной сетке и на равномерной сетке
		//////////////////////////////////////////////////////////////

		// Вывод координат:
		for (int i = 0; i < Col; i++)
		{
			Linde_rez_y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1);
			//Linde_rez_x = Length / (List - 2) / 2.0 + Length / (List - 2.0) * (i - 1);
			if (i == 0 || i == Col - 1)
			{
				std::cout << fixed << setprecision(6) << 0.0 << "\t";
				file << fixed << setprecision(6) << 0.0 << "\t";
			}
			else
			{
				std::cout << fixed << setprecision(6) << Linde_rez_y << "\t";
				file << fixed << setprecision(6) << Linde_rez_y << "\t";
			}
		}
		std::cout << endl;
		// Вывод результата программы
		SetColor(14, 0);
		std::cout << "Implicit: (" << center_x << ", y, " << center_z << ") : \n";
		file << "Implicit: (" << center_x << ", y, " << center_z << ") : \n";
		for (int i = 0; i < Col; i++)
		{
			Linde_rez_y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1);
			//Linde_rez_x = Length / (List - 2) / 2.0 + Length / (List - 2.0) * (i - 1);
			if (i == 0 || i == Col - 1)
			{
				std::cout << fixed << setprecision(6) << Test1.GetT(center_x, 0, center_z) << "\t";
				file << fixed << setprecision(6) << Test1.GetT(center_x, 0, center_z) << "\t";
			}
			else
			{
				std::cout << fixed << setprecision(6) << Test1.GetT(center_x, Linde_rez_y, center_z) << "\t";
				file << fixed << setprecision(6) << Test1.GetT(center_x, Linde_rez_y, center_z) << "\t";
			}
		}
		std::cout << endl;
		// Вывод точного результата
		SetColor(10, 0);
		std::cout << "Analit: (" << center_x << ", y, " << center_z << ") : \n";
		file << "\nAnalit: (" << center_x << ", y, " << center_z << ") : \n";
		for (int i = 0; i < Col; i++)
		{
			Linde_rez_y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1);
			//Linde_rez_x = Length / (List - 2) / 2.0 + Length / (List - 2.0) * (i - 1);
			Index = Get_Index_4(center_x, Linde_rez_y, center_z, Row, Col, List, L_side);

			if (i == 0 || i == Col - 1)
			{
				std::cout << fixed << setprecision(6) << 0.0 << "\t";
				file << fixed << setprecision(6) << 0.0 << "\t";
			}
			else
			{
				std::cout << fixed << setprecision(6) << Analit(center_x, Linde_rez_y, center_z, Time, Number_obl[Index]) << "\t";
				file << fixed << setprecision(6) << Analit(center_x, Linde_rez_y, center_z, Time, Number_obl[Index]) << "\t";
			}
		}
		std::cout << endl;
		SetColor(15, 0);

		std::cout << "Coordinate (0.5,y,0.5)  |\tImplicit\t|\tImplicit_t\t|\tAnalit  \t|\tAbs. pogr\t|\tAbs_t. pogr\n";
		std::cout << "____________________________________________________________________________________________________________________" << std::endl;
		file << "Coordinate (0.5,y,0.5)  |\tImplicit\t|\tImplicit_t\t|\tAnalit  \t|\tAbs. pogr\t|\tAbs_t. pogr\n";
		file << "____________________________________________________________________________________________________________________" << std::endl;
		sumPogr = 0, sumPogrOtn = 0;
		double Im_t, sumPogr_t = 0;
		for (int i = 0; i <= 20; i++)
		{
			Linde_rez_y = Length / 20 * i;
			//Linde_rez_x = Length / (List - 2) / 2.0 + Length / (List - 2.0) * (i - 1);
			Index = Get_Index_4(center_x, Linde_rez_y, center_z, Row, Col, List, L_side);

			if (i == 0 || i == 20)
			{
				std::cout << fixed << setprecision(8) << "\t" << 0.0 << "\t|\t";
				SetColor(14, 0);
				std::cout << 0.0 << "\t|\t";
				std::cout << 0.0 << "\t|\t";
				SetColor(10, 0);
				std::cout << 0.0 << "\t|\t";
				SetColor(15, 0);
				std::cout << 0.0 << "\t|\t" << 0.0 << "\n";

				file << fixed << setprecision(8) << "\t" << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t";
				file << 0.0 << "\t|\t" << 0.0 << "\n";
			}
			else
			{
				std::cout << fixed << setprecision(8) << "\t" << Length / 20 * i << "\t|\t";
				SetColor(14, 0);
				Im = Test1.GetT(center_x, Length / 20 * i, center_y);
				std::cout << Im << "\t|\t";
				Im_t = Test1.Get_adapt_T(center_x, Linde_rez_y, center_z);
				std::cout << Im_t << "\t|\t";
				SetColor(10, 0);
				An = Analit(center_x, Linde_rez_y, center_z, Time, Number_obl[Index]);
				std::cout << An << "\t|\t";
				SetColor(15, 0);
				std::cout << abs(An - Im) << "\t|\t" << abs(An - Im_t) << "\n";

				file << fixed << setprecision(8) << "\t" << Length / 20 * i << "\t|\t";
				file << Im << "\t|\t";
				file << Im_t << "\t|\t";
				file << An << "\t|\t";
				file << abs(An - Im) << "\t|\t" << abs(An - Im_t) << "\n";

				sumPogr += abs(An - Im);
				sumPogr_t += abs(An - Im_t);
			}
		}
		std::cout << endl;
		std::cout << "\nСредняя абсолютная погрешность = " << sumPogr / 20;
		std::cout << "\nСредняя абсолютная погрешность адаптивной сетки = " << sumPogr_t / 20 << endl << endl;
		file << endl;
		file << "\nСредняя абсолютная погрешность = " << sumPogr / 20;
		file << "\nСредняя оносительная погрешность = " << sumPogr_t / 20 << endl << endl;

		// Вывод среднеквадратичной погрешности
		if (Task == 1)
		{
			double sredPogr{};
			double x_length, y_length, z_length;
			x_length = 1;
			double shag = 1. / 19;
			Coord_4 xyz(0, 0, 0);
			for (int i = 1; i < 20; i++)
			{
				xyz.z = shag * i;
				for (int j = 1; j < 20; j++)
				{
					xyz.y = shag * j;
					for (int k = 1; k < 20; k++)
					{
						xyz.x = shag * k;
						sredPogr += pow(Test1.Get_adapt_T(xyz.x, xyz.y, xyz.z) - Analit(xyz.x, xyz.y, xyz.z, Time, 0), 2.);
					}

				}
			}
			sredPogr = sqrt(sredPogr / pow(19, 6));
			std::cout << endl;
			std::cout << "\nСреднеквадратичная погрешность = " << sredPogr << endl;
			file << "\nСреднеквадратичная погрешность = " << sredPogr << endl;
		}
		else if (Task == 4)
		{
			double center_z = ((Row - 2) / 2.0) * L_side;
			double sredPogr{};
			double x_length, y_length, z_length;
			x_length = 1;
			double shag_x = 10.2 / 20;
			double shag_y = 3.3 / 20;
			Coord_4 xyz(0, 0, 0);
			for (int j = 1; j < 20; j++)
			{
				xyz.y = shag_y * j;
				for (int k = 1; k < 20; k++)
				{
					xyz.x = shag_x * k;
					sredPogr += pow(Test1.Get_adapt_T(xyz.x, xyz.y, center_z) - Analit(xyz.x, xyz.y, center_z, Time, 0), 2.);
				}

			}
			sredPogr = sqrt(sredPogr / pow(19, 3));
			std::cout << endl;
			std::cout << "\nСреднеквадратичная погрешность = " << sredPogr << endl;
			file << "\nСреднеквадратичная погрешность = " << sredPogr << endl;
		}
		//for (int i = 0; i < 20; i++)
		//{
		//	for (int j = 0; j < 20; j++)
		//	{

		//	}
		//}

	}
	file.close();

	/////////	Вывод для графиков:		///////

	// 1d
	ofstream file1("Grafik 1d.txt");
	if (!file1.is_open())
	{
		std::cout << "Error open file test for Grafik 1d!\n";
		exit(0);
	}
	else
	{
		for (int i = 0; i < Col; i++)
		{
			if (i == 0)
				file1 << 0.0 << "\t" << 0.0 << '\n';
			else if (i == Col - 1)
				file1 << 1.0 << "\t" << 0.0 << '\n';
			else
				file1 << Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1) << "\t" << Test1.GetT(0.5, Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (i - 1), 0.5) << '\n';
		}
	}
	file1.close();

	// 2d
	ofstream file2("Grafik 2d.txt");
	if (!file2.is_open())
	{
		std::cout << "Error open file test for Grafik 2d!\n";
		exit(0);
	}
	else
	{
		double center_z = ((Row - 2) / 2.0) * L_side;
		for (int i = 0; i < List; i++)
		{
			for (int j = 0; j < Col; j++)
			{
				int index = 1 + j * Row + i * Row * Col;
				Coord_4 xyz = Get_Coord_4(1, j, i, Row, Col, List, L_side);
				if (Test1.Tipe_Obl(xyz.x, xyz.y, center_z) == 0)
					file2 << xyz.x << "\t" << xyz.y << "\t" << "NaN" << "\n";
				else
					file2 << xyz.x << "\t" << xyz.y << "\t" << Test1.GetT(xyz.x, xyz.y, center_z) << "\n";
			}
			file2 << "\n";
		}
	}
	file2.close();

	// 3d
	if (Task != 1)
		return 0;
	ofstream file3("Grafik 3d.txt");
	if (!file3.is_open())
	{
		std::cout << "Error open file test for Grafik 3d!\n";
		exit(0);
	}
	else
	{
		// Разрез по середине:
		/*for (int i = 0; i < Row; i++)
		{
			double x;
			if (i == 0)
				x = 0;
			else if (i == Row - 1)
				x = 1;
			else
				x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
			for (int j = 0; j < Col; j++)
			{
				double y;
				if (j == 0)
					y = 0;
				else if (j == Col - 1)
					y = 1;
				else
					y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
				for (int k = 0; k < List; k++)
				{
					double z;
					if (k == 0)
						z = 0;
					else if (k == List - 1)
						z = 1;
					else
						z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
					if ((i == 0 || i == Row - 1) && j <= Col / 2)
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) << "\n";
					}

				}
				file3 << "\n";
			}
			file3 << "\n";
		}


		for (int j = 0; j < Col; j++)
		{
			double y;
			if (j == 0)
				y = 0;
			else if (j == Col - 1)
				y = 1;
			else
				y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
			for (int k = 0; k < List; k++)
			{
				double z;
				if (k == 0)
					z = 0;
				else if (k == List - 1)
					z = 1;
				else
					z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
				for (int i = 0; i < Row; i++)
				{
					double x;
					if (i == 0)
						x = 0;
					else if (i == Row - 1)
						x = 1;
					else
						x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
					if (j == 0 || j == Col / 2)
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) << "\n";
					}
				}
				file3 << "\n";
			}
			file3 << "\n";
		}

		for (int k = 0; k < List; k++)
		{
			double z;
			if (k == 0)
				z = 0;
			else if (k == List - 1)
				z = 1;
			else
				z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
			for (int j = 0; j < Col; j++)
			{
				double y;
				if (j == 0)
					y = 0;
				else if (j == Col - 1)
					y = 1;
				else
					y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
				for (int i = 0; i < Row; i++)
				{
					double x;
					if (i == 0)
						x = 0;
					else if (i == Row - 1)
						x = 1;
					else
						x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
					if ((k == List - 1) && j <= Col / 2)
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) << "\n";
					}

				}
				file3 << "\n";
			}
			file3 << "\n";
		}*/

		for (int i = 0; i < Row; i++)
		{
			double x;
			if (i == 0)
				x = 0;
			else if (i == Row - 1)
				x = 1;
			else
				x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
			for (int j = 0; j < Col; j++)
			{
				double y;
				if (j == 0)
					y = 0;
				else if (j == Col - 1)
					y = 1;
				else
					y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
				for (int k = 0; k < List; k++)
				{
					double z;
					if (k == 0)
						z = 0;
					else if (k == List - 1)
						z = 1;
					else
						z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
					if ((i == Row - 1 && j <= Col / 2) || (i == Row / 2 && j >= Col / 2) || i == 0)
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) + 0.03 << "\n";
					}

				}
				file3 << "\n";
			}
			file3 << "\n";
		}


		for (int j = 0; j < Col; j++)
		{
			double y;
			if (j == 0)
				y = 0;
			else if (j == Col - 1)
				y = 1;
			else
				y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
			for (int k = 0; k < List; k++)
			{
				double z;
				if (k == 0)
					z = 0;
				else if (k == List - 1)
					z = 1;
				else
					z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
				for (int i = 0; i < Row; i++)
				{
					double x;
					if (i == 0)
						x = 0;
					else if (i == Row - 1)
						x = 1;
					else
						x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
					if ((j == Col - 1 && i <= Row / 2) || j == 0 || (j == Col / 2 && i >= Row / 2))
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) + 0.03 << "\n";
					}
				}
				file3 << "\n";
			}
			file3 << "\n";
		}

		for (int k = 0; k < List; k++)
		{
			double z;
			if (k == 0)
				z = 0;
			else if (k == List - 1)
				z = 1;
			else
				z = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (k - 1);
			for (int j = 0; j < Col; j++)
			{
				double y;
				if (j == 0)
					y = 0;
				else if (j == Col - 1)
					y = 1;
				else
					y = Length / (Col - 2) / 2.0 + Length / (Col - 2.0) * (j - 1);
				for (int i = 0; i < Row; i++)
				{
					double x;
					if (i == 0)
						x = 0;
					else if (i == Row - 1)
						x = 1;
					else
						x = Length / (Row - 2) / 2.0 + Length / (Row - 2.0) * (i - 1);
					if ((k == List - 1) && !(i > Row / 2 && j > Col / 2))
					{
						file3 << x << "\t" << y << "\t" << z << "\t" << Test1.GetT(x, y, z) + 0.03 << "\n";
					}

				}
				file3 << "\n";
			}
			file3 << "\n";
		}

	}
	file3.close();

	return 0;
}