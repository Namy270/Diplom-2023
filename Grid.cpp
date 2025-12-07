#include "Grid.h"
#include <fstream>
#include <iomanip>
#include<iostream>
#define sigma 1370e-10 //Постоянная Стефана-Больцмана
#define C 3e10		 //Скорость света

struct Coord
{
    double x, y, z;
    Coord(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};
Coord Get_Coord(int Row, int Col, int List, int Rows, int Cols, int Lists, double L_side);
double Get_Index(double x, double y, double z, int Row, int Col, int List, double L_side);
void show_loading(int i, int Finish_iter, time_t start_time_1);

Grid::Grid(int ROW, int COL, int LIST)
{
    this->ROW = ROW;
    this->COL = COL;
    this->LIST = LIST;
    int size = ROW * COL * LIST;
    Max_Gradient = 10;
    Array.resize(size);
    dt = 0;
    sum_cell = 0;
    T = 0;
}

void Grid::init(double* T, double* Num_obl, int* check_oblast, double* ro, double L_side, int* type_bondary,
    double (**arrF)(double, double, double, double), double(**arrFi)(double, double, double, double), double gradient)
{
    Max_Gradient = gradient;
    for (int i = 0; i < ROW * COL * LIST; i++)
    {
        this->Array[i].Number_obl = Num_obl[i];
        this->Array[i].F_edge = arrF[i];
        this->Array[i].Fi_edge = arrFi[i];
        this->Array[i].type_of_boundary_conditions = type_bondary[i];
        this->Array[i].check_oblast = check_oblast[i];
        this->Array[i].T = T[i];
        this->Array[i].E0 = Get_E(Array[i].T, Array[i].Number_obl);
        this->Array[i].KAPPA = Get_Kappa(Array[i].T, Array[i].Number_obl);
        this->Array[i].L_side = L_side;
        this->Array[i].S_side = L_side * L_side;
        this->Array[i].ro = ro[i];
        this->Array[i].M = ro[i] * this->Array[i].S_side * L_side;
        this->Array[i].lvl = 1;
    }
}

double Grid::Get_Kappa(double T, int type)
{
    if (type == 0)
        return 1;
    else if (type == 1)
        return 35.3417911268 * pow(T, 6) * 4;
    else if (type == 2)
        return 7.06135882537 * pow(T, 6) * 4;
    else
        return 0.0706135882537 * pow(T, 6) * 4;
}
double Grid::Get_E(double T, int type)
{
    if (type == 0)
        return T;
    else if (type == 1)
        return 10 * T + pow(T, 4);
    else if (type == 2)
        return 1.23762376238 * T + 0.123762376238 * pow(T, 4);
    else
        return 4.952475247524 * T + 0.4952475247524 * pow(T, 4);
}
double Grid::Get_fi(double T, int type)
{
    if (type == 0)
        return 1;
    else if (type == 1)
        return 10 + 4 * pow(T, 3);
    else if (type == 2)
        return 1.23762376238 + 4 * 0.123762376238 * pow(T, 3);
    else
        return 4.952475247524 + 4 * 0.4952475247524 * pow(T, 3);
}

void Grid::Calculate(double dt, double T, std::string solv, int iters)
{
    this->dt = dt;
    this->T = T;
    int size = ROW * COL * LIST;
    double* vector_T;

    // Для замера скорости
    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    int show; // для вывода на экран (не заморачивайся)
    time_t start_time_1;
    time(&start_time_1);
    int Finish_iter = round(T / dt);
    for (int i = 0; i <= Finish_iter; i++)
    {
        double time_solve = i * dt;
        /* if (time_solve >= 0.4 && time_solve <= 0.5)
         {
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
                 for (int iii = 0; iii < COL; iii++)
                 {
                     double center_z = ((ROW - 2) / 2.0) * Array[0].L_side;
                     double center_x = ((LIST - 2) / 2.0) * Array[0].L_side;

                     Coord Y = Get_Coord(center_z, iii, center_x, ROW, COL, LIST, Array[0].L_side);

                     int Index_6 =  Get_Index(6, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);
                     int Index_8 =  Get_Index(8, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);
                     int Index_10 = Get_Index(10, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);

                     file11 << std::fixed << std::setprecision(6) << Y.y + 98 << "\t|\t" << Array[Index_6].T << "\t|\t" << Array[Index_8].T << "\t|\t" << Array[Index_10].T << "\n";

                 }
             }
             file11.close();
         }*/
         // Выводим статистику загрузки
        if (ROW >= 50 && Finish_iter >= 100)
            show = Finish_iter / 100;
        else
            show = Finish_iter / 20;
        if (i % show == 0)
            show_loading(i, Finish_iter, start_time_1);

        // Дробим сетку
        begin = std::chrono::steady_clock::now();
        Build_Adapt_grid();
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        Time_program.adapt_time += elapsed_ms.count() / 1000.;

        for (int w = 0; w < 1; w++) // итерации по нелинейности
        {
            // Выделяем память пот СЛАУ
            begin = std::chrono::steady_clock::now();
            Matrix.clear();
            vector_B.clear();
            Matrix.resize(sum_cell);
            vector_B.resize(sum_cell);
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            Time_program.memory_time += elapsed_ms.count() / 1000.;

            // Cобираем коэффициенты 
            begin = std::chrono::steady_clock::now();
            Build_SLAY_coef(time_solve);
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            Time_program.build_slay_time += elapsed_ms.count() / 1000.;

            // Преобразуем в формат CSR
            begin = std::chrono::steady_clock::now();
            Matrix.CSR();
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            Time_program.CSR_time += elapsed_ms.count() / 1000.;

            // Решаем
            begin = std::chrono::steady_clock::now();
            if (solv == "Host")
                vector_T = Solve_SLAY_host(Matrix.getVal(), Matrix.getRowPtr(), Matrix.getColInd(), vector_B, iters);
            else
                vector_T = Solve_SLAY_device(Matrix.getVal(), Matrix.getRowPtr(), Matrix.getColInd(), vector_B, iters);
            for (int j = 0; j < size; j++)
            {
                if (Array[j].next_LVL.empty())
                {
                    Array[j].T = vector_T[j];
                    Array[j].KAPPA = Get_Kappa(Array[j].T, Array[j].Number_obl); // Задать вопрос научнику
                }
                else
                {
                    Array[j].T = 0;
                    for (int ii = 0; ii < 8; ii++)
                    {
                        Array[j].next_LVL[ii].T = vector_T[Array[j].next_LVL[ii].id];
                        Array[j].next_LVL[ii].KAPPA = Get_Kappa(Array[j].next_LVL[ii].T, Array[j].next_LVL[ii].Number_obl);
                        Array[j].T += Array[j].next_LVL[ii].T;
                    }
                    Array[j].T /= 8;
                    Array[j].KAPPA = Get_Kappa(Array[j].T, Array[j].Number_obl);
                }
            }
            end = std::chrono::steady_clock::now();
            elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            Time_program.Slay_time += elapsed_ms.count() / 1000.;
        }

        //std::ofstream file12("T = 1, X = 6, X = 8, X = 10.txt");
        //if (!file12.is_open())
        //{
        //    printf("Error open file11!");
        //    exit(1);
        //}
        //else
        //{
        //    file12 << "\t  Y  \t|\t  X = 5.8  \t|\t  X = 7.8  \t|\tX = 9.8\n";
        //    file12 << "______________________________________________________________________________________\n";
        //    for (int iii = 0; iii < COL; iii++)
        //    {
        //        double center_z = ((ROW - 2) / 2.0) * Array[0].L_side;
        //        double center_x = ((LIST - 2) / 2.0) * Array[0].L_side;

        //        Coord Y = Get_Coord(center_z, iii, center_x, ROW, COL, LIST, Array[0].L_side);

        //        int Index_6 = Get_Index(6, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);
        //        int Index_8 = Get_Index(8, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);
        //        int Index_10 = Get_Index(10, Y.y, center_z, ROW, COL, LIST, Array[0].L_side);

        //        file12 << std::fixed << std::setprecision(6) << Y.y + 98 << "\t|\t" << Array[Index_6].T << "\t|\t" << Array[Index_8].T << "\t|\t" << Array[Index_10].T << "\n";

        //    }
        //}
        //file12.close();

        Grafik_adapt();
        Grafik_2d();


        // Записываем результаты
        begin = std::chrono::steady_clock::now();
        for (int j = 0; j < size; j++)
        {
            if (Array[j].next_LVL.empty())
            {
                Array[j].T = vector_T[j];
                Array[j].KAPPA = Get_Kappa(Array[j].T, Array[j].Number_obl);
                Array[j].E0 = Get_E(Array[j].T, Array[j].Number_obl);
            }
            else
            {
                Array[j].T = 0;
                for (int ii = 0; ii < 8; ii++)
                {
                    Array[j].next_LVL[ii].T = vector_T[Array[j].next_LVL[ii].id];
                    Array[j].next_LVL[ii].KAPPA = Get_Kappa(Array[j].next_LVL[ii].T, Array[j].next_LVL[ii].Number_obl);
                    Array[j].next_LVL[ii].E0 = Get_E(Array[j].next_LVL[ii].T, Array[j].next_LVL[ii].Number_obl);
                    Array[j].T += Array[j].next_LVL[ii].T;
                }
                Array[j].T /= 8;
                Array[j].KAPPA = Get_Kappa(Array[j].T, Array[j].Number_obl);
                Array[j].E0 = Get_E(Array[j].T, Array[j].Number_obl);
            }
        }
        end = std::chrono::steady_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        Time_program.result_time += elapsed_ms.count() / 1000.;
    }
	// Запись матрицы в файл
// Запись матрицы в файл
//        if (i == 1)
//        {
//            std::vector<std::vector<double>> matrix_(vector_B.size());
//            for (int ii = 0; ii < vector_B.size(); ii++)
//                matrix_[ii].resize(vector_B.size());
//            matrix_ = Build_SLAY_coef_2(time_solve);
//            std::ofstream file_matrix("Matrix_test.txt");
//            if (!file_matrix.is_open())
//            {
//                std::cout << "Error file2!\n";
//                exit(1);
//            }
//            else
//            {
//                file_matrix << vector_B.size() << "\n";
//                for (int jj = 0; jj < vector_B.size(); jj++)
//                {
//                    file_matrix << vector_B[jj] << "\t";
//                }
//                file_matrix << "\n";
//                for (int jj = 0; jj < vector_B.size(); jj++)
//                {
//                    for (int kk = 0; kk < vector_B.size(); kk++)
//                        file_matrix << matrix_[jj][kk]<< "\t";
//                    file_matrix << "\n";
//                }
//            }
//        }
//    }

    Time_program.sum_time = Time_program.Slay_time + Time_program.build_slay_time + Time_program.adapt_time + Time_program.result_time + Time_program.memory_time + Time_program.CSR_time;
    delete[] vector_T;
}

void Grid::Build_Adapt_grid()
{
    int UPPER, DOWN, LEFT, RIGHT, FRONT, BACK, counter;
    int id = ROW * COL * LIST;
    for (int k = 0; k < LIST; k++)
    {
        for (int j = 0; j < COL; j++)
        {
            for (int i = 0; i < ROW; i++)
            {
                counter = i + j * ROW + k * ROW * COL; // номер рассматриваемой ячейки
                DOWN = counter - 1;
                UPPER = counter + 1;
                LEFT = counter - ROW;
                RIGHT = counter + ROW;
                FRONT = counter - ROW * COL;
                BACK = counter + ROW * COL;

                if (Array[counter].check_oblast < 1) // Рассматриваем только внутренние ячейки
                    continue;

                // Идем по ячейкам и смотрим градиент рассматриваемой ячейки с верхней, задней и правой 
                if (Array[UPPER].check_oblast > 0)
                {
                    if (abs(Array[counter].T - Array[UPPER].T) > Max_Gradient)  //Если gradient м\у текущей и верхней ячейкой > заданного
                    {
                        Array[UPPER].split_cell = true;
                        Array[counter].split_cell = true;

                    }
                    else if (!(Array[counter].next_LVL.empty() && Array[UPPER].next_LVL.empty()))
                    {
                        if (Array[counter].next_LVL.empty() > Array[UPPER].next_LVL.empty()) // если сверху раздробленная ячейка, а рассматриваемая нет
                        {
                            for (int ii = 0; ii < 8; ii += 2)
                                if (abs(Array[counter].T - Array[UPPER].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[UPPER].split_cell = true;
                                    goto met1;
                                }
                        }
                        else if (Array[counter].next_LVL.empty() < Array[UPPER].next_LVL.empty()) // если рассматриваемая  ячейка раздроблена, а сверху нет
                        {
                            for (int ii = 1; ii < 8; ii += 2)
                                if (abs(Array[counter].next_LVL[ii].T - Array[UPPER].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[UPPER].split_cell = true;
                                    goto met1;
                                }
                        }
                        else  // если сверху раздробленная ячейка и рассматриваемая тоже
                        {
                            for (int ii = 0; ii < 8; ii += 2)
                                if (abs(Array[counter].next_LVL[ii + 1].T - Array[UPPER].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[UPPER].split_cell = true;
                                    goto met1;
                                }
                        }
                    }
                }
            met1:
                if (Array[RIGHT].check_oblast > 0)
                {
                    if (abs(Array[counter].T - Array[RIGHT].T) > Max_Gradient)  //Если gradient м\у текущей и правой ячейкой > заданного
                    {
                        Array[RIGHT].split_cell = true;
                        Array[counter].split_cell = true;
                    }
                    else if (!(Array[counter].next_LVL.empty() && Array[RIGHT].next_LVL.empty()))
                    {
                        if (Array[counter].next_LVL.empty() > Array[RIGHT].next_LVL.empty()) // если справа раздробленная ячейка, а рассматриваемая нет
                        {
                            for (int ii = 0; ii < 6; ii += 1 + 2 * (ii == 1))
                                if (abs(Array[counter].T - Array[RIGHT].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[RIGHT].split_cell = true;
                                    goto met2;
                                }
                        }
                        else if (Array[counter].next_LVL.empty() < Array[RIGHT].next_LVL.empty()) // если рассматриваемая  ячейка раздроблена, а сверху нет
                        {
                            for (int ii = 2; ii < 8; ii += 1 + 2 * (ii == 3))
                                if (abs(Array[counter].next_LVL[ii].T - Array[RIGHT].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[RIGHT].split_cell = true;
                                    goto met2;
                                }
                        }
                        else  // если сверху раздробленная ячейка и рассматриваемая тоже
                        {
                            for (int ii = 0; ii < 6; ii += 1 + 2 * (ii == 1))
                                if (abs(Array[counter].next_LVL[ii + 2].T - Array[RIGHT].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[RIGHT].split_cell = true;
                                    goto met2;
                                }
                        }
                    }
                }
            met2:
                if (Array[BACK].check_oblast > 0)
                {
                    if (abs(Array[counter].T - Array[BACK].T) > Max_Gradient)  //Если gradient м\у текущей и задней ячейкой > заданного
                    {
                        Array[BACK].split_cell = true;
                        Array[counter].split_cell = true;
                    }
                    else if (!(Array[counter].next_LVL.empty() && Array[BACK].next_LVL.empty()))
                    {
                        if (Array[counter].next_LVL.empty() > Array[BACK].next_LVL.empty()) // если сзади раздробленная ячейка, а рассматриваемая нет
                        {
                            for (int ii = 0; ii < 4; ii++)
                                if (abs(Array[counter].T - Array[BACK].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[BACK].split_cell = true;
                                    goto met3;
                                }
                        }
                        else if (Array[counter].next_LVL.empty() < Array[BACK].next_LVL.empty()) // если рассматриваемая  ячейка раздроблена, а сзади нет
                        {
                            for (int ii = 4; ii < 8; ii++)
                                if (abs(Array[counter].next_LVL[ii].T - Array[BACK].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[BACK].split_cell = true;
                                    goto met3;
                                }
                        }
                        else  // если сзади раздробленная ячейка и рассматриваемая тоже
                        {
                            for (int ii = 0; ii < 4; ii++)
                                if (abs(Array[counter].next_LVL[ii + 4].T - Array[BACK].next_LVL[ii].T) > Max_Gradient)
                                {
                                    Array[counter].split_cell = true;
                                    Array[BACK].split_cell = true;
                                    goto met3;
                                }
                        }
                    }
                }
            met3:

                // Разбиваем ячейку, если можно
                if (Array[counter].split_cell)   // ДОРАБОТАТЬ ПРИ УСЛОВИИ, ЧТО они уже разбиты!!!!!
                {
                    if (Array[counter].lvl < max_lvl)
                    {
                        if (Array[counter].next_LVL.empty())
                        {
                            Array[counter].next_LVL.resize(8);
                            for (int ii = 0; ii < 8; ii++)
                            {
                                Array[counter].next_LVL[ii].T = Array[counter].T;
                                Array[counter].next_LVL[ii].Number_obl = Array[counter].Number_obl;
                                Array[counter].next_LVL[ii].KAPPA = Array[counter].KAPPA;
                                Array[counter].next_LVL[ii].E0 = Array[counter].E0;
                                Array[counter].next_LVL[ii].L_side = Array[counter].L_side / 2;
                                Array[counter].next_LVL[ii].S_side = Array[counter].S_side / 4;
                                Array[counter].next_LVL[ii].lvl = Array[counter].lvl + 1;
                                Array[counter].next_LVL[ii].ro = Array[counter].ro;
                                Array[counter].next_LVL[ii].M = Array[counter].M / 8;
                                Array[counter].next_LVL[ii].check_oblast = Array[counter].check_oblast;
                                Array[counter].next_LVL[ii].id = id;
                                Array[counter].next_LVL[ii].parent_id = counter;
                                id++;
                            }
                        }
                        else
                        {
                            for (int ii = 0; ii < 8; ii++)
                            {
                                Array[counter].next_LVL[ii].id = id;
                                id++;
                            }
                        }
                    }
                    Array[counter].split_cell = false;
                }
                else
                {
                    if (!Array[counter].next_LVL.empty())
                    {
                        Array[counter].next_LVL.clear();
                    }
                }
            }
        }
    }
    sum_cell = id;
}

void Grid::Build_SLAY_coef(double t_)
{
    int UPPER, DOWN, LEFT, RIGHT, FRONT, BACK, counter;
    for (int k = 0; k < LIST; k++)
    {
        for (int j = 0; j < COL; j++)
        {
            for (int i = 0; i < ROW; i++)
            {
                counter = i + j * ROW + k * ROW * COL; // номер рассматриваемой ячейки
                DOWN = counter - 1;
                UPPER = counter + 1;
                LEFT = counter - ROW;
                RIGHT = counter + ROW;
                FRONT = counter - ROW * COL;
                BACK = counter + ROW * COL;
                if (Array[counter].check_oblast == 0)         // Если ячейка внешняя
                {
                    Matrix.add(counter, counter, 1);
                    vector_B[counter] = 0;
                }
                else if (Array[counter].check_oblast < 0)     // Если ячейка граничная
                {
                    switch (Array[counter].type_of_boundary_conditions) // Выбираем построение коэф исходя из заданного гран условия
                    {
                    case(1):    // ПЕРВАЯ КРАЕВАЯ
                        Matrix.add(counter, counter, 1);
                        vector_B[counter] = Array[counter].F_edge((LIST - k - 1) * Array[0].L_side, j * Array[0].L_side, i * Array[0].L_side, t_); // Подумать!!!
                        break;
                    case(2):    // ЖЁСКАЯ СТЕНКА
                        Matrix.add(counter, counter, 1);
                        vector_B[counter] = 0;
                        break;
                    case(3):    // W(t) = fi
                        Matrix.add(counter, counter, 1);
                        vector_B[counter] = 0;
                        break;
                    case(4):    // ИЗЛУЧЕНИЕ ЧЁРНОГО ТЕЛА //    ТУТ НУЖНО ПОДУМАТЬ ЕЩЁ!!!!!!!!!!!!!!
                        int index_inside;
                        if (UPPER >= 0)
                            if (Array[UPPER].check_oblast > 0)
                                index_inside = UPPER;
                        //Chek_side = "UP";
                        if (DOWN >= 0)
                            if (Array[DOWN].check_oblast > 0)
                                index_inside = DOWN;
                        // Chek_side = "DOWN";
                        if (LEFT >= 0)
                            if (Array[LEFT].check_oblast > 0)
                                index_inside = LEFT;
                        //Chek_side = "LEFT";
                        if (RIGHT >= 0)
                            if (Array[RIGHT].check_oblast > 0)
                                index_inside = RIGHT;
                        //Chek_side = "RIGHT";
                        if (FRONT >= 0)
                            if (Array[FRONT].check_oblast > 0)
                                index_inside = FRONT;
                        //Chek_side = "FRONT";
                        if (BACK >= 0)
                            if (Array[BACK].check_oblast > 0)
                                index_inside = BACK;
                        //Chek_side = "BACK";
                        //Matrix(counter, counter) = sigma * C * pow(Array[counter].T, 3) + KAPPA_edge(counter, index_inside, -1, -1);
                        //Matrix(counter, index_inside) = -KAPPA_edge(counter, index_inside, -1, -1) / L_(counter, index_inside);
                        Matrix.add(counter, counter, sigma * C * pow(Array[counter].T, 3) + KAPPA_edge(counter, index_inside, -1, -1));
                        Matrix.add(counter, index_inside, -KAPPA_edge(counter, index_inside, -1, -1) / L_(counter, index_inside));
                        vector_B[counter] = 75e-2 * sigma * C * pow(Array[counter].T, 4);
                        break;
                    default:
                        printf("\nIndefinite boundary condition!\n");
                        exit(1);
                    }
                }
                else                                    //Если ячейка внутренняя
                {
                    Coord xyz = Get_Coord(i, j, k, ROW, COL, LIST, Array[counter].L_side);
                    int aU = 1, aD = 1, aR = 1, aL = 1, aB = 1, aF = 1;
                    double add_B = 0;
                    {   // Проверка на пригроничность 
                            //!!! Тут можно добавить флаг, чтобы не дробить приграничные ячейки (просто если дробить, то это сложно реализовать)
                        if (Array[UPPER].check_oblast < 0)
                        {
                            if (Array[UPPER].type_of_boundary_conditions == 2)
                                aU = 0;
                            else if (Array[UPPER].type_of_boundary_conditions == 3)
                            {
                                aU = 0;
                                add_B += Array[UPPER].Fi_edge(xyz.x, xyz.y, xyz.z + Array[counter].L_side / 2., t_); // Над этим надо подумать, это неправильно
                            }
                        }

                        if (Array[DOWN].check_oblast < 0)
                        {
                            if (Array[DOWN].type_of_boundary_conditions == 2)
                                aD = 0;
                            else if (Array[DOWN].type_of_boundary_conditions == 3)
                            {
                                aD = 0;
                                add_B += Array[DOWN].Fi_edge(xyz.x, xyz.y, xyz.z - Array[counter].L_side / 2., t_);
                            }
                        }

                        if (Array[RIGHT].check_oblast < 0)
                        {
                            if (Array[RIGHT].type_of_boundary_conditions == 2)
                                aR = 0;
                            else if (Array[RIGHT].type_of_boundary_conditions == 3)
                            {
                                aR = 0;
                                add_B += Array[RIGHT].Fi_edge(xyz.x, xyz.y + Array[counter].L_side / 2., xyz.z, t_);
                            }
                        }

                        if (Array[LEFT].check_oblast < 0)
                        {
                            if (Array[LEFT].type_of_boundary_conditions == 2)
                                aL = 0;
                            else if (Array[LEFT].type_of_boundary_conditions == 3)
                            {
                                aL = 0;
                                add_B += Array[LEFT].Fi_edge(xyz.x, xyz.y - Array[counter].L_side / 2., xyz.z, t_);
                            }
                        }

                        if (Array[BACK].check_oblast < 0)
                        {
                            if (Array[BACK].type_of_boundary_conditions == 2)
                                aB = 0;
                            else if (Array[BACK].type_of_boundary_conditions == 3)
                            {
                                aB = 0;
                                add_B += Array[BACK].Fi_edge(xyz.x - Array[counter].L_side / 2., xyz.y, xyz.z, t_);
                            }
                        }

                        if (Array[FRONT].check_oblast < 0)
                        {
                            if (Array[FRONT].type_of_boundary_conditions == 2)
                                aF = 0;
                            else if (Array[FRONT].type_of_boundary_conditions == 3)
                            {
                                aF = 0;
                                add_B += Array[FRONT].Fi_edge(xyz.x + Array[counter].L_side / 2., xyz.y, xyz.z, t_);
                            }
                        }
                    }
                    ////////////////////////////////////////////////////////////////////////////////
                    /// Собираем коэффициенты!
                    ////////////////////////////////////////////////////////////////////////////////
                    if (Array[counter].next_LVL.empty()) // Если текущая ячейка целая
                    {
                        double help_add;
                        double MatrixAdd_counter_counter = 0;
                        //  1)
                        if (!Array[UPPER].next_LVL.empty()) // Если верхняя ячейка раздроблена
                        {
                            //Matrix(counter, UPPER) = 0;    // ????? уточнить у научника
                            for (int ii = 0; ii < 8; ii += 2)
                            {
                                help_add = -(dt * S_edge(counter, UPPER) * KAPPA_edge(counter, ii, -1, UPPER)) / L_(counter, UPPER);
                                Matrix.add(counter, Array[UPPER].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }
                        }
                        else    // Если верхняя ячейка целая
                        {
                            help_add = aU * (-(dt * S_edge(counter, UPPER) * KAPPA_edge(counter, UPPER, -1, -1)) / L_(counter, UPPER));
                            Matrix.add(counter, UPPER, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }
                        //  2)
                        if (!Array[DOWN].next_LVL.empty()) // Если нижняя ячейка раздроблена
                        {
                            for (int ii = 1; ii < 8; ii += 2)
                            {
                                help_add = -(dt * S_edge(counter, DOWN) * KAPPA_edge(counter, ii, -1, DOWN)) / L_(counter, DOWN);
                                Matrix.add(counter, Array[DOWN].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }

                        }
                        else    // Если нижняя ячейка целая
                        {
                            help_add = aD * (-(dt * S_edge(counter, DOWN) * KAPPA_edge(counter, DOWN, -1, -1)) / L_(counter, DOWN));
                            Matrix.add(counter, DOWN, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }
                        //  3)
                        if (!Array[RIGHT].next_LVL.empty()) // Если правая ячейка раздроблена
                        {
                            for (int ii = 0; ii < 6; ii += 1 + 2 * (ii == 1))
                            {
                                help_add = -(dt * S_edge(counter, RIGHT) * KAPPA_edge(counter, ii, -1, RIGHT)) / L_(counter, RIGHT);
                                Matrix.add(counter, Array[RIGHT].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }
                        }
                        else    // Если правая ячейка целая
                        {
                            help_add = aR * (-(dt * S_edge(counter, RIGHT) * KAPPA_edge(counter, RIGHT, -1, -1)) / L_(counter, RIGHT));
                            Matrix.add(counter, RIGHT, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }

                        //  4)
                        if (!Array[LEFT].next_LVL.empty()) // Если левая ячейка раздроблена
                        {
                            for (int ii = 2; ii < 8; ii += 1 + 2 * (ii == 3))
                            {
                                help_add = -(dt * S_edge(counter, LEFT) * KAPPA_edge(counter, ii, -1, LEFT)) / L_(counter, LEFT);
                                Matrix.add(counter, Array[LEFT].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }
                        }
                        else    // Если левая ячейка целая
                        {
                            help_add = aL * (-(dt * S_edge(counter, LEFT) * KAPPA_edge(counter, LEFT, -1, -1)) / L_(counter, LEFT));
                            Matrix.add(counter, LEFT, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }

                        //  5)
                        if (!Array[BACK].next_LVL.empty()) // Если задняя ячейка раздроблена
                        {
                            for (int ii = 0; ii < 4; ii++)
                            {
                                help_add = -(dt * S_edge(counter, BACK) * KAPPA_edge(counter, ii, -1, BACK)) / L_(counter, BACK);
                                Matrix.add(counter, Array[BACK].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }
                        }
                        else    // Если зядняя ячейка целая
                        {
                            help_add = aB * (-(dt * S_edge(counter, BACK) * KAPPA_edge(counter, BACK, -1, -1)) / L_(counter, BACK));
                            Matrix.add(counter, BACK, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }

                        //  6)
                        if (!Array[FRONT].next_LVL.empty()) // Если передняя ячейка раздроблена
                        {
                            for (int ii = 4; ii < 8; ii++)
                            {
                                help_add = -(dt * S_edge(counter, FRONT) * KAPPA_edge(counter, ii, -1, FRONT)) / L_(counter, FRONT);
                                Matrix.add(counter, Array[FRONT].next_LVL[ii].id, help_add);
                                MatrixAdd_counter_counter -= help_add;
                            }
                        }
                        else    // Если левая ячейка целая
                        {
                            help_add = aF * (-(dt * S_edge(counter, FRONT) * KAPPA_edge(counter, FRONT, -1, -1)) / L_(counter, FRONT));
                            Matrix.add(counter, FRONT, help_add);
                            MatrixAdd_counter_counter -= help_add;
                        }
                        MatrixAdd_counter_counter += Array[counter].M * Get_fi(Array[counter].T, Array[counter].Number_obl);
                        Matrix.add(counter, counter, MatrixAdd_counter_counter);
                        vector_B[counter] = Array[counter].M * (Array[counter].T * Get_fi(Array[counter].T, Array[counter].Number_obl) + Array[counter].E0 - Get_E(Array[counter].T, Array[counter].Number_obl)) + add_B * dt;
                    }
                    else
                    {
                        int down, upper, left, right, front, back, upper_2, down_2, right_2, left_2, front_2, back_2;
                        int start_id = Array[counter].next_LVL[0].id;
                        int now_id;
                        //Matrix(counter, UPPER) = aU * (-(dt * S_edge(counter, UPPER) * KAPPA_edge(counter, UPPER)) / L_(counter, UPPER));
                        //Matrix(counter, DOWN) = aD * (-(dt * S_edge(counter, DOWN) * KAPPA_edge(counter, DOWN)) / L_(counter, DOWN));
                        //Matrix(counter, RIGHT) = aR * (-(dt * S_edge(counter, RIGHT) * KAPPA_edge(counter, RIGHT)) / L_(counter, RIGHT));
                        //Matrix(counter, LEFT) = aL * (-(dt * S_edge(counter, LEFT) * KAPPA_edge(counter, LEFT)) / L_(counter, LEFT));
                        //Matrix(counter, BACK) = aB * (-(dt * S_edge(counter, BACK) * KAPPA_edge(counter, BACK)) / L_(counter, BACK));
                        //Matrix(counter, FRONT) = aF * (-(dt * S_edge(counter, FRONT) * KAPPA_edge(counter, FRONT)) / L_(counter, FRONT));

                        for (int ii = 0; ii < 8; ii++) // Проходим по каждой раздробленной ячейки
                        {
                            double help_sum = 0;
                            double MatrixAdd_nowid_nowid = 0;
                            now_id = start_id + ii;
                            upper = now_id + 1;
                            down = now_id - 1;
                            right = now_id + 2;
                            left = now_id - 2;
                            front = now_id - 4;
                            back = now_id + 4;

                            upper_2 = ii + 1;
                            down_2 = ii - 1;
                            right_2 = ii + 2;
                            left_2 = ii - 2;
                            front_2 = ii - 4;
                            back_2 = ii + 4;

                            // Коэффициент при верхней ячейке
                            if (ii == 0 || ii == 2 || ii == 4 || ii == 6)  // Если верхняя ячейка раздроблена внутри текущей 
                            {
                                help_sum = aU * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, upper_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, upper, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[UPPER].next_LVL.empty()) // Если верхняя ячейка не раздроблена
                                {
                                    help_sum = aU * (-(dt * S_edge(counter, UPPER) * KAPPA_edge(ii, UPPER, counter, -1)) / L_(counter, UPPER));
                                    Matrix.add(now_id, UPPER, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else // Если верхняя ячейка раздроблена и находится в верхней ячейке
                                {
                                    upper = Array[UPPER].next_LVL[down_2].id;
                                    help_sum = aU * (-(dt * S_edge(counter, UPPER) * KAPPA_edge(ii, down_2, counter, UPPER)) / L_(counter, UPPER));
                                    Matrix.add(now_id, upper, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            // Коэффициент при нижней ячейке
                            if (ii == 1 || ii == 3 || ii == 5 || ii == 7)  // Если нижняя ячейка раздроблена внутри текущей 
                            {
                                help_sum = aD * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, down_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, down, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[DOWN].next_LVL.empty())   // Если нижняя ячейка не раздроблена
                                {
                                    help_sum = aD * (-(dt * S_edge(counter, DOWN) * KAPPA_edge(ii, DOWN, counter, -1)) / L_(counter, DOWN));
                                    Matrix.add(now_id, DOWN, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else    // Если нижняя ячейка раздроблена и находится в нижней ячейке
                                {
                                    down = Array[DOWN].next_LVL[upper_2].id;
                                    help_sum = aD * (-(dt * S_edge(counter, DOWN) * KAPPA_edge(ii, upper_2, counter, DOWN)) / L_(counter, DOWN));
                                    Matrix.add(now_id, down, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            // Коэффициент при правой ячейке
                            if (ii == 0 || ii == 1 || ii == 4 || ii == 5)    // Если правая ячейка раздроблена внутри текущей 
                            {
                                help_sum = aR * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, right_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, right, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[RIGHT].next_LVL.empty()) // Если правая ячейка не раздроблена
                                {
                                    help_sum = aR * (-(dt * S_edge(counter, RIGHT) * KAPPA_edge(ii, RIGHT, counter, -1)) / L_(counter, RIGHT));
                                    Matrix.add(now_id, RIGHT, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else    // Если правая ячейка раздроблена и находится в правой ячейке
                                {
                                    right = Array[RIGHT].next_LVL[left_2].id;
                                    help_sum = aR * (-(dt * S_edge(counter, RIGHT) * KAPPA_edge(ii, left_2, counter, RIGHT)) / L_(counter, RIGHT));
                                    Matrix.add(now_id, right, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            // Коэффициент при левой ячейке
                            if (ii == 2 || ii == 3 || ii == 6 || ii == 7)  // Если левая ячейка раздроблена внутри текущей 
                            {
                                help_sum = aL * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, left_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, left, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[LEFT].next_LVL.empty())   // Если левая ячейка не раздроблена
                                {
                                    help_sum = aL * (-(dt * S_edge(counter, LEFT) * KAPPA_edge(ii, LEFT, counter, -1)) / L_(counter, LEFT));
                                    Matrix.add(now_id, LEFT, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else    // Если левая ячейка раздроблена и находится в левой ячейке
                                {
                                    left = Array[LEFT].next_LVL[right_2].id;
                                    help_sum = aL * (-(dt * S_edge(counter, LEFT) * KAPPA_edge(ii, right_2, counter, LEFT)) / L_(counter, LEFT));
                                    Matrix.add(now_id, left, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            // Коэффициент при задней ячейке
                            if (ii == 0 || ii == 1 || ii == 2 || ii == 3) // Если задняя ячейка раздроблена внутри текущей 
                            {
                                help_sum = aB * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, back_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, back, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[BACK].next_LVL.empty())   // Если задняя ячейка не раздроблена
                                {
                                    help_sum = aB * (-(dt * S_edge(counter, BACK) * KAPPA_edge(ii, BACK, counter, -1)) / L_(counter, BACK));
                                    Matrix.add(now_id, BACK, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else     // Если задняя ячейка раздроблена и находится в задней ячейке
                                {
                                    back = Array[BACK].next_LVL[front_2].id;
                                    help_sum = aB * (-(dt * S_edge(counter, BACK) * KAPPA_edge(ii, front_2, counter, BACK)) / L_(counter, BACK));
                                    Matrix.add(now_id, back, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            // Коэффициент при передней ячейке
                            if (ii == 4 || ii == 5 || ii == 6 || ii == 7)  // Если передняя ячейка раздроблена внутри текущей 
                            {
                                help_sum = aF * (-(dt * S_edge(counter, counter) * KAPPA_edge(ii, front_2, counter, counter)) / L_(counter, counter));
                                Matrix.add(now_id, front, help_sum);
                                MatrixAdd_nowid_nowid -= help_sum;
                            }
                            else
                            {
                                if (Array[FRONT].next_LVL.empty())  // Если передняя ячейка не раздроблена
                                {
                                    help_sum = aF * (-(dt * S_edge(counter, FRONT) * KAPPA_edge(ii, FRONT, counter, -1)) / L_(counter, FRONT));
                                    Matrix.add(now_id, FRONT, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }
                                else     // Если передняя ячейка раздроблена и находится в передней ячейке
                                {
                                    front = Array[FRONT].next_LVL[back_2].id;
                                    help_sum = aF * (-(dt * S_edge(counter, FRONT) * KAPPA_edge(ii, back_2, counter, FRONT)) / L_(counter, FRONT));
                                    Matrix.add(now_id, front, help_sum);
                                    MatrixAdd_nowid_nowid -= help_sum;
                                }

                            }
                            MatrixAdd_nowid_nowid += Array[counter].next_LVL[ii].M * Get_fi(Array[counter].next_LVL[ii].T, Array[counter].next_LVL[ii].Number_obl);
                            Matrix.add(now_id, now_id, MatrixAdd_nowid_nowid);
                            vector_B[now_id] = add_B * dt + Array[counter].next_LVL[ii].M *
                                (Array[counter].next_LVL[ii].T * Get_fi(Array[counter].next_LVL[ii].T, Array[counter].next_LVL[ii].Number_obl)
                                    + Array[counter].next_LVL[ii].E0 - Get_E(Array[counter].next_LVL[ii].T, Array[counter].next_LVL[ii].Number_obl));
                        }
                        Matrix.add(counter, counter, 1);
                        vector_B[counter] = 0;
                    }
                }
            }
        }
    }
}

double Grid::KAPPA_edge(int index_1, int index_2, int parant_1, int parant_2)
{
    if (parant_1 == -1 && parant_2 == -1) // Если у обоих родителей нет
        return (Array[index_1].KAPPA + Array[index_2].KAPPA) / 2;
    else if (parant_1 >= 0 && parant_2 >= 0) // Если у обоих есть родители
        return (Array[parant_1].next_LVL[index_1].KAPPA + Array[parant_2].next_LVL[index_2].KAPPA) / 2;
    else if (parant_1 == -1 && parant_2 >= 0) // Если у второго есть родитель
        return (Array[index_1].KAPPA + Array[parant_2].next_LVL[index_2].KAPPA) / 2;
    else  // Если у первого есть родитель
        return (Array[parant_1].next_LVL[index_1].KAPPA + Array[index_2].KAPPA) / 2;
}

double Grid::S_edge(int index_1, int index_2)
{
    if (Array[index_1].next_LVL.empty() && Array[index_2].next_LVL.empty())
        return Array[index_1].S_side;
    else
        if (!Array[index_1].next_LVL.empty())
            return Array[index_1].next_LVL[0].S_side;
        else
            return Array[index_2].next_LVL[0].S_side;

    //if (abs(index_1 - index_2) == 1)
    //    return Array[index_1].S_sideZ;
    //if (abs(index_1 - index_2) == ROW)
    //    return Array[index_1].S_sideY;
    //if (abs(index_1 - index_2) == ROW * COL)
    //    return Array[index_1].S_sideX;
    //else
    //{
    //    printf("Error of double Grid::S_edge(int index_1, int index_2)");
    //    exit(-1);
    //}
}

double Grid::L_(int index_1, int index_2)
{
    double L = Array[index_1].L_side;
    double gran_dop;
    if (Array[index_1].check_oblast < 0 || Array[index_2].check_oblast < 0)
        gran_dop = L / 2.0;
    else
        gran_dop = 0;

    if (Array[index_1].next_LVL.empty() && Array[index_2].next_LVL.empty())
        return L - gran_dop;
    else if (!Array[index_1].next_LVL.empty() && !Array[index_2].next_LVL.empty())
        return L / 2.0;
    else
        if (Array[index_1].check_oblast < 0 || Array[index_2].check_oblast < 0)
            return L / 4.0;
        else
            return  sqrt((L / 2.0) * (L / 2.0) + (L / 4.0) * (L / 4.0) + (3 * L / 4.0) * (3 * L / 4.0));


    //if (abs(index_1 - index_2) == 1)
    //    return (koef1 * Array[index_1].L_sideZ + koef2 * Array[index_2].L_sideZ) / 2.0;
    //else if (abs(index_1 - index_2) == ROW)
    //    return (koef1 * Array[index_1].L_sideY + koef2 * Array[index_2].L_sideY) / 2.0;
    //else if(abs(index_1 - index_2) == ROW * COL)
    //    return (koef1 * Array[index_1].L_sideX + koef2 * Array[index_2].L_sideX) / 2.0;
    //else
    //{
    //    printf("Error of double Grid::L_edge(int index_1, int index_2");
    //    exit(-2);
    //}
}

double Grid::GetT(double x, double y, double z)
{
    int inner_cell = ROW * COL + ROW + 1;
    double X_index = ((LIST - 2) * Array[inner_cell].L_side - x) / Array[inner_cell].L_side;
    if (x == ((double)LIST - 2) * Array[inner_cell].L_side)
        X_index = LIST - 1.0;
    double Y_index = y / Array[inner_cell].L_side;
    if (y == ((double)COL - 2) * Array[inner_cell].L_side)
        Y_index = COL - 1.0;
    double Z_index = z / Array[inner_cell].L_side;
    if (z == ((double)ROW - 2) * Array[inner_cell].L_side)
        Z_index = ROW - 1.0;
    int Index = ROW * COL * int(ceil(X_index)) + ROW * int(ceil(Y_index)) + int(ceil(Z_index));
    return Array[Index].T;
}

double Grid::Get_adapt_T(double x, double y, double z)
{
    double Lenght = Array[0].L_side * (LIST - 2);
    x = Lenght - x;
    double L = Array[0].L_side;
    int X_index = ceil(x / L);
    /*  if (x == ((double)LIST - 2) * L)
          X_index = LIST - 1.0;*/
    int Y_index = ceil(y / L);
    //if (y == ((double)COL - 2) * L)
    //    Y_index = COL - 1.0;
    int Z_index = ceil(z / L);
    //if (z == ((double)ROW - 2) * L)
    //    Z_index = ROW - 1.0;
    int Index = ROW * COL * X_index + ROW * Y_index + Z_index;
    //return Array[Index].T;

    if (Array[Index].next_LVL.empty())
        return Array[Index].T;
    else
    {
        double x_adapt, y_adapt, z_adapt;
        if (x > L * X_index - L / 2)
            x_adapt = 1;
        else if (x < L * X_index - L / 2)
            x_adapt = 0;
        else
            x_adapt = (x > Lenght / 2);

        if (y > L * Y_index - L / 2)
            y_adapt = 1;
        else if (y < L * Y_index - L / 2)
            y_adapt = 0;
        else
            y_adapt = (y > Lenght / 2);

        if (z > L * Z_index - L / 2)
            z_adapt = 1;
        else if (z < L * Z_index - L / 2)
            z_adapt = 0;
        else
            z_adapt = (z > Lenght / 2);
        double index_adapt = 4 * x_adapt + 2 * y_adapt + z_adapt;
        return  Array[Index].next_LVL[index_adapt].T;
    }
}

int Grid::Tipe_Obl(double x, double y, double z)
{
    int inner_cell = ROW * COL + ROW + 1;
    double X_index = ((LIST - 2) * Array[inner_cell].L_side - x) / Array[inner_cell].L_side;
    if (x == ((double)LIST - 2) * Array[inner_cell].L_side)
        X_index = LIST - 1.0;
    double Y_index = y / Array[inner_cell].L_side;
    if (y == ((double)COL - 2) * Array[inner_cell].L_side)
        Y_index = COL - 1.0;
    double Z_index = z / Array[inner_cell].L_side;
    if (z == ((double)ROW - 2) * Array[inner_cell].L_side)
        Z_index = ROW - 1.0;
    int Index = ROW * COL * int(ceil(X_index)) + ROW * int(ceil(Y_index)) + int(ceil(Z_index));
    return Array[Index].check_oblast;
}

Grid::~Grid()
{
}

void Grid::Grafik_adapt() //Переделать
{
    std::ofstream Grafik_adapt("Grafik_adapt.txt");
    if (!Grafik_adapt.is_open())
    {
        std::cout << "Error open file Grafik_adapt!\n";
        exit(0);
    }
    else
    {
        // y == z; x == y
        double Lenght_x = Array[0].L_side * (LIST - 2);
        double Lenght_y = Array[0].L_side * (COL - 2);
        double center_z = (ROW) / 2.0;
        double L_x = Lenght_x / LIST;
        double L_y = Lenght_y / COL;

        for (int k = 0; k < LIST; k++)
        {
            Grafik_adapt << k * L_x << "\t" << 0 << std::endl;
            Grafik_adapt << k * L_x << "\t" << Lenght_y << "\n\n" << std::endl;

            for (int j = 0; j < COL; j++)
            {
                if (k == 0)
                {
                    Grafik_adapt << 0 << "\t" << j * L_y << std::endl;
                    Grafik_adapt << Lenght_x << "\t" << j * L_y << "\n\n" << std::endl;
                }
                int counter = center_z + j * ROW + COL * ROW * (LIST - k - 1); // номер рассматриваемой ячейки
                if (Array[counter].next_LVL.empty())
                    continue;
                else
                {
                    Grafik_adapt << k * L_x << "\t" << j * L_y + L_y / 2.0 << std::endl;
                    Grafik_adapt << (k + 1) * L_x << "\t" << j * L_y + L_y / 2.0 << "\n\n" << std::endl;

                    Grafik_adapt << k * L_x + L_x / 2.0 << "\t" << j * L_y << std::endl;
                    Grafik_adapt << k * L_x + L_x / 2.0 << "\t" << (j + 1) * L_y << "\n\n" << std::endl;
                }
            }
        }
        Grafik_adapt << Lenght_x << "\t" << 0 << std::endl;
        Grafik_adapt << Lenght_x << "\t" << Lenght_y << "\n\n" << std::endl;


        Grafik_adapt << 0 << "\t" << Lenght_y << std::endl;
        Grafik_adapt << Lenght_x << "\t" << Lenght_y << "\n\n" << std::endl;
    }
    Grafik_adapt.close();
}

void Grid::Grafik_2d()
{
    std::ofstream file2("Grafik 2d.txt");
    if (!file2.is_open())
    {
        std::cout << "Error open file test for Grafik 2d!\n";
        exit(0);
    }
    else
    {
        //double Length = Array[0].L_side * (LIST - 2);
        //for (int i = 0; i < ROW; i++)
        //{
        //    double x;
        //    if (i == 0)
        //        x = 0;
        //    else if (i == ROW - 1)
        //        x = 1;
        //    else
        //        x = Length / (ROW - 2) / 2.0 + Length / (ROW - 2.0) * (i - 1);
        //    for (int j = 0; j < COL; j++)
        //    {
        //        double y;
        //        if (j == 0)
        //            y = 0;
        //        else if (j == COL - 1)
        //            y = 1;
        //        else
        //            y = Length / (COL - 2) / 2.0 + Length / (COL - 2.0) * (j - 1);
        //        file2 << x << "\t" << y << "\t" << GetT(x, y, 0.5) << "\n";
        //    }
        //    file2 << "\n";
        //}
        double center_z = ((ROW - 2) / 2.0) * Array[0].L_side;
        for (int i = 0; i < LIST; i++)
        {
            for (int j = 0; j < COL; j++)
            {
                int index = 1 + j * ROW + i * ROW * COL;
                Coord xyz = Get_Coord(1, j, i, ROW, COL, LIST, Array[index].L_side);
                if (Tipe_Obl(xyz.x, xyz.y, center_z) == 0)
                    file2 << xyz.x << "\t" << xyz.y << "\t" << "NaN" << "\n";
                else
                    file2 << xyz.x << "\t" << xyz.y << "\t" << GetT(xyz.x, xyz.y, center_z) << "\n";
            }
            file2 << "\n";
        }

    }
    file2.close();
}

Coord Get_Coord(int Row, int Col, int List, int Rows, int Cols, int Lists, double L_side)
{
    double x, y, z;
    double addX = 0, addY = 0, addZ = 0;
    if (Row == 0)
        addZ = L_side / 2;
    else if (Row == Rows - 1)
        addZ = -L_side / 2;
    if (Col == 0)
        addY = L_side / 2;
    else if (Col == Cols - 1)
        addY = -L_side / 2;
    if (List == 0)
        addX = L_side / 2;
    else if (List == Lists - 1)
        addX = -L_side / 2;
    double x_length = (Lists - 2) * L_side;
    z = L_side * Row - L_side / 2;
    y = L_side * Col - L_side / 2;
    x = x_length - (L_side * List - L_side / 2);
    return Coord(x, y + addY, z + addZ);
}
double Get_Index(double x, double y, double z, int Row, int Col, int List, double L_side)
{
    int inner_cell = Row * Col + Row + 1;
    double X_index = ((List - 2) * L_side - x) / L_side;
    double Y_index = y / L_side;
    double Z_index = z / L_side;
    int Index = Row * Col * int(ceil(X_index)) + Row * int(ceil(Y_index)) + int(ceil(Z_index));
    return Index;
}

void show_loading(int i, int Finish_iter, time_t start_time_1)
{
    time_t end_time_1;
    time(&end_time_1);
    double seconds = difftime(end_time_1, start_time_1);
    double show2 = double(i) / Finish_iter * 100;

    system("cls");
    std::cout << "____________________  " << round(show2) << "%" << std::endl;
    for (int l = 0; l < round(show2) / 5; l++)
        std::cout << "*";
    int finishtime = seconds / double(i) * (Finish_iter - i);
    int hours = finishtime / 3600;
    int min = (finishtime % 3600) / 60;
    seconds = ((finishtime % 3600) % 60);
    if (i)
        std::cout << "\nFinish: " << hours << " : " << min << " : " << seconds << std::endl;

}