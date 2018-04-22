#pragma once

void shiftMSE(const Image &I, const Image &I1, const Image &I2, int shiftV, int shiftH, long long &v1, long long &v2);
void MSE(const Image &I1, const Image &I2, const Image &I3, const int maxShift, int &minStepY1, int &minStepX1, int &minStepY2, int &minStepX2);

long long shiftMSE(const Image &I1, const Image &I2, int shiftV, int shiftH);
void MSE_big(const Image &I1, const Image &I2, const Image &I3, const int maxShift, int &minStepY1, int &minStepX1, int &minStepY2, int &minStepX2);

#include "metrics.hpp"