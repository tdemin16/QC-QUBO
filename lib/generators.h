#ifndef GENERATORS_H
#define GENERATORS_H

#include "lib.h"

int number_partitioning_problem(MatrixXf &Q, string file = "../test/npp16.txt");
void quadratic_assignment(MatrixXf &Q, string file);

#endif