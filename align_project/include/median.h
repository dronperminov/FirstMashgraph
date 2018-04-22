#pragma once

#include <vector>
#include <deque>
#include "pixel.h"

/*
	Функции для вычисления медианы двумя способами:
	- calcMedianSquare - квадратичная сложность O(r^2)
	- calcMedianLinear - линейная сложность O(r)

	Для вычисления медианы ха константное время требуются дополнительные классы:
	- RGBhistogram - гистограмма из трёх каналов
	- ConstantMedian - класс, на основе гистограммы подсчитывающий медиану за константное время
*/

uint calcMedianSquare(std::vector<uint> &vector);
uint calcMedianLinear(std::vector<uint> &vector);

Pixel calcMedian(const Image &I, uint y, uint x, int radius);

struct RGBhistogram {
    double rgb[3][256];
    std::deque<Pixel> env;

    uint width;
    uint height;

    RGBhistogram(const Image &I);
};

struct ConstantMedian {
    long long rgb[3][256], old[3][256];
    std::vector<RGBhistogram> columns;
    uint pos, radius;

    ConstantMedian(const Image &I, uint rad);

    void shiftRight();
    void shiftRight(const Pixel &p);
    void shiftDown(const Image &row);
    Pixel getMedian();
};

#include "median.hpp"