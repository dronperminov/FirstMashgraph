uint Pixel::limit(double v) {
    if (v > 255)
        return 255;

    if (v < 0)
        return 0;

    return v;
}

Pixel::Pixel() : rgb(std::make_tuple(0, 0, 0)) {}

Pixel::Pixel(const std::tuple<uint, uint, uint> &tuple) : rgb(std::make_tuple(limit(std::get<0>(tuple)), limit(std::get<1>(tuple)), limit(std::get<2>(tuple)))) {}

Pixel::Pixel(double r, double g, double b) : rgb(std::make_tuple(limit(r), limit(g), limit(b))) {}

Pixel Pixel::operator*(double v) {
    double r = std::get<0>(rgb);
    double g = std::get<1>(rgb);
    double b = std::get<2>(rgb);

    return Pixel(r * v, g * v, b * v);
}

Pixel Pixel::operator+(const Pixel &p) {
    double r = std::get<0>(rgb) + std::get<0>(p.rgb);
    double g = std::get<1>(rgb) + std::get<1>(p.rgb);
    double b = std::get<2>(rgb) + std::get<2>(p.rgb);

    return Pixel(r, g, b);
}

Pixel Pixel::operator-(double v) {
    double r = std::get<0>(rgb) - v;
    double g = std::get<1>(rgb) - v;
    double b = std::get<2>(rgb) - v;

    return Pixel(r, g, b);
}

Pixel Pixel::scaleRGB(double r, double g, double b) {
    return Pixel(std::get<0>(rgb) * r, std::get<1>(rgb) * g, std::get<2>(rgb) * b);
}

uint Pixel::getBrightness() {
    return limit(0.2125 * std::get<0>(rgb) + 0.7154 * std::get<1>(rgb) + 0.0721 * std::get<2>(rgb));
}

std::tuple<uint, uint, uint> Pixel::getRGB() const {
    return rgb;
}

uint Pixel::R() const { return std::get<0>(rgb); }
uint Pixel::G() const { return std::get<1>(rgb); }
uint Pixel::B() const { return std::get<2>(rgb); }