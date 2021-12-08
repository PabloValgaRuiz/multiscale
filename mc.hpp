#include "MobMatrix.hpp"
#include "MCDistMult.hpp"
#include "MarkovDistMult.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <future>

static void iteracionesMC(const MobMatrix& T, double &infectados, double beta, double pI, double pC, int tiempo, int i);
int mc(const MobMatrix& T, const std::string state);