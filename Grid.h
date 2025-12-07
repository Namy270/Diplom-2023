#pragma once
#include<vector>
#include <iomanip>
#include <chrono>
#include "SparseMatrix.h"
#include "cuSolve.cuh"
#include "SolveHost.h"
#include <string>
#include <omp.h>
// Структура, запоминающее ремя работы программы
struct Timer
{
	double adapt_time = 0;          // таймер дробления ячеек
	double memory_time = 0;         // таймер выделения памяти
	double build_slay_time = 0;     // таймер построения коэффициентов слау
	double Slay_time = 0;           // таймер решения слау
	double CSR_time = 0;            // таймер преобразование в формат СSR
	double result_time = 0;         // таймер записи результатов
	double sum_time = 0;
};
// Структура, описывающаяя ячейку сетки:
struct Cell
{
	int type_of_boundary_conditions = 0;				// Тип граничных условий (1 - первая краевая; 2 - Жесткая стенка; 3 - W(t) = fi; 4 - Чёрное тело)
	double (*F_edge)(double, double, double, double) = nullptr;  // Функция источника (обычно для граничной ячейки с условием 1-го рода)
	double (*Fi_edge)(double, double, double, double) = nullptr; // Функция для граничной ячейки с выделением энергии W(t)
	double Number_obl;						// Номер области
	double T;								// Температура
	double E0;								// Внутренняя энергия
	double KAPPA;							// Коэффициент теплопроводности
	double L_side;		// Радиус ячейки (по умолчанию = 1)
	double S_side;		// Площадь грани (по умолчанию = 1)
	double ro;								// Плотность ячейки
	double M;								// Масса ячейки
	int lvl = 1;							// Уровень адаптивности ячейки (по умолчанию = 1)
	int check_oblast;						// Тип области:  0 - внешняя ячейка; -1 граничная ячейка;
											//				 1 - первая область
											//				....
	int id;
	int parent_id;
	std::vector<Cell> next_LVL;
	bool split_cell = false;
};


// Класс сетки:
class Grid
{
private:
	double Max_Gradient;
	int max_lvl = 2;
	int sum_cell;		//Количество Ячеек (для подсчета с адаптивностью)
	int ROW, COL, LIST;						// Размерность сетки (количество строк, столбцов, слоев)
	std::vector<Cell> Array;			// Массив ячеек
	SparseMatrix Matrix;				// Матрица коэффициентов СЛАУ
	std::vector<double> vector_B;	   // Вектор свободных членов СЛАУ
	double dt, T;		  // Шаг по времени и конечное время соответсвенно

public:
	Timer Time_program;
	Grid(int ROW, int COL, int LIST); // Конструктор с заданием размерности сетки
	void init(double* T, double* Num_obl, int* check_oblast, double* ro, double L_side, int* type_bondary,
		double (**arrF)(double, double, double, double), double(**arrFi)(double, double, double, double), double gradient); // Задание начальных данных и информации об областях
	double Get_Kappa(double T, int type);	  //расчёт коэф. теплопроводности
	double Get_E(double T, int type);		 // расчет внутренней энергии
	double Get_fi(double T, int type);		// рассчет производной фи
	void Calculate(double dt, double T, std::string solv, int iters); // решть с заданным шагом по времени, типом решателя(процессор или видеокарта) и количеством итераций в решателе СЛАУ
	void Build_SLAY_coef(double t_);									 // Строим систему линейных уравнений
	void Build_Adapt_grid();								// Дробим сетку по заданному градиенту
	double KAPPA_edge(int index_1, int index_2, int parant_1, int parant_2);		   // Получить коэффициент теплопроводности на границе двух ячеек
	double S_edge(int index_1, int index_2);			  // Получить площадь общей грани двух ячеек	
	double L_(int index_1, int index_2);				 // Получить расстояние между центрами двух ячейк 
	double GetT(double x, double y, double z);
	double Get_adapt_T(double x, double y, double z);
	int Tipe_Obl(double x, double y, double z);
	void Grafik_adapt();
	void Grafik_2d();
	~Grid(); // Деструктор
};

