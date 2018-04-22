#pragma once

/*
	Класс пикселя для более удобной работы с пикселями изображений
	- Создаётся из трёх вещественных значений или из кортежа
	- Умножение пикселя на число - все значения умножатся на него
	- Прибавление пикселя
	- Вычитание числа
	- Изменение каждого канала по-отельности
	- Получение яркости
	- Получение кортежа для изображения
	- Получение каждого канала по-отдельности
*/

class Pixel {
	std::tuple<uint, uint, uint> rgb;

	uint limit(double v);

public:
	Pixel();
	Pixel(const std::tuple<uint, uint, uint> &tuple);
	Pixel(double r, double g, double b);

	Pixel operator*(double v);
	Pixel operator+(const Pixel &p);
	Pixel operator-(double v);

	Pixel scaleRGB(double r, double g, double b);
	uint getBrightness();

	std::tuple<uint, uint, uint> getRGB() const;
	uint R() const;
	uint G() const;
	uint B() const;
};

#include "pixel.hpp"