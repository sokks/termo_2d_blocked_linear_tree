#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>


// описание области
namespace Area {
    extern double x_start;
    extern double y_start;
    extern double x_end;
    extern double y_end;

    extern double ro;
    extern double lambda;
    extern double c;
    extern double a;

    double T0(double x, double y);
    double BorderCond1(double x, double y, double t);
    double BorderCond2(double x, double y, double t);
    double BorderCond3(double x, double y, double t);
    double Q(double x, double y, double t);

    enum Border { LEFT, RIGHT, UP, DOWN };
    void get_border_cond(Border b, char *cond_type, double (**cond_func)(double, double, double));

    extern int max_grid_levels;
    extern double refine_thresholds[];

    int Refine1(double x, double y);
    int Refine2(double x, double y);
    int Refine3(double x, double y);
}
