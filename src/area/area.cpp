#include "area.h"

// описание области
double Area::x_start = 0.0;
double Area::y_start = 0.0;
double Area::x_end = 1.0;
double Area::y_end = 1.0;

// double ro = 1190.0, lambda = 0.16, c = 900.0;
double Area::ro = 1, Area::lambda = 0.1, Area::c = 1;
double Area::a = lambda / (ro * c);
// double a = 1;



// начальное распределение температуры и граничные условия
double Area::T0(double x, double y) {
    // double r = 0.05;
    // double c_x = 0.25;
    // double c_y = 0.25;
    // if ( (x-c_x) * (x-c_x) + (y-c_y) * (y-c_y) < r*r ) { //,0078125 < 0,09
    //     return 10;
    // }
    double sigma_x = 0.05;
    double sigma_y = 0.05;
    double mean_x = 0.25;
    double mean_y = 0.25;
    return 0.05 * (1.0 / (sigma_x * sigma_y)) * exp( - ((x - mean_x)*(x - mean_x))/(2*sigma_x*sigma_x) - ((y - mean_y)*(y - mean_y))/(2*sigma_y*sigma_y) );
}

int Area::Refine1(double x, double y) {
    // double r = 0.24;
    // double c_x = 0.25;
    // double c_y = 0.25;
    // // if ( (x-c_x) * (x-c_x) + (y-c_y) * (y-c_y) < r*r ) { //,0078125 < 0,09
    // if ( (fabs(x-c_x) < r) && (fabs(y - c_y) < r) ) {
    //     return 1;
    // }

    return (x < (x_start + (x_end - x_start) / 2)) && (y < (y_start + (y_end - y_start) / 2) );
}

int Area::Refine2(double x, double y) {
    // double r = 0.11;
    // double c_x = 0.25;
    // double c_y = 0.25;
    // // if ( (x-c_x) * (x-c_x) + (y-c_y) * (y-c_y) < r*r ) { //,0078125 < 0,09
    // if ( (fabs(x-c_x) < r) && (fabs(y - c_y) < r) ) {
    //     return 1;
    // }
    // return 0.0;

    return (x < (x_start + (x_end - x_start) / 4)) && (y < (y_start + (y_end - y_start) / 4) );
}

int Area::Refine3(double x, double y) {
    double r = 0.05;
    double c_x = 0.25;
    double c_y = 0.25;
    // if ( (x-c_x) * (x-c_x) + (y-c_y) * (y-c_y) < r*r ) { //,0078125 < 0,09
    if ( (fabs(x-c_x) < r) && (fabs(y - c_y) < r) ) {
        return 1;
    }
    return 0.0;
}

double Area::BorderCond1(double x, double y, double t) { // down & left
    return 0.0;
}

double Area::BorderCond2(double x, double y, double t) { // right
    return -2 * t;
}

double Area::BorderCond3(double x, double y, double t) { // up
    return 100 * t;
}

double Area::BorderCondFixed1(double x, double y, double t) { // up
    return 1;
}

void Area::get_border_cond(Border b, char *cond_type, double (**cond_func)(double, double, double)) {

    switch (b) {
        case LEFT: {
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
        case RIGHT: {
            // *cond_type = 2;
            // *cond_func = &BorderCond2;
            *cond_type = 1;
            // *cond_func = &BorderCond1;
            *cond_func = &BorderCondFixed1;
            break;
        }
        case UP: {
            // *cond_type = 2;
            // *cond_func = &BorderCond3;
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
        case DOWN: {
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
    }
}



double Area::Q(double x, double y, double t) {
    if (t > 0.2) {
        return 0;
    }
    double eps = 0.1;

    double r = 0.05;
    double c_x = 0.25;
    double c_y = 0.25;

    // if ((fabs(x - 0.45) < eps) && (fabs(y - 0.5) < eps)) {
    //     return 0.01;
    // }
    // if ((fabs(x - 0.55) < eps) && (fabs(y - 0.5) < eps)) {
    //     return 0.01;
    // }

    if ((fabs(x - c_x) < eps ) && (fabs(y - c_y) < eps)) {
        return 0.01;
    }

    return 0;
}
