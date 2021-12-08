#pragma once

#include "MobMatrix.hpp"
#include "MarkovDistMult.hpp"
#include <Eigen/Core>
std::pair<Eigen::ArrayXd, double> iteracion_p0(const MobMatrix& T, const MarkovDistMult& markov);
std::pair<Eigen::ArrayXd, double> iteracion(const MobMatrix& T, const MarkovDistMult& markov);