#include "align.h"
#include "pixel.h"
#include "median.h"
#include "metrics.h"

#include <string>

using std::string;
using std::cout;
using std::endl;

bool needMirror = true;

/*********************************************************************************************************************************************************************/
void setChannel(Image &dst, const Image &color, uint ch, int stepV, int stepH) {
    uint y0 = stepV > 0 ? stepV : 0;
    uint y1 = dst.n_rows + (stepV > 0 ? 0 : stepV);

    uint x0 = stepH > 0 ? stepH : 0;
    uint x1 = dst.n_cols + (stepH > 0 ? 0 : stepH);

    for (uint i = y0; i < y1; i++) {
        for (uint j = x0; j < x1; j++) {
            if (ch == 0) {
                std::get<0>(dst(i, j)) = std::get<0>(color(i - stepV, j - stepH));
            }
            else if (ch == 1) {
                std::get<1>(dst(i, j)) = std::get<0>(color(i - stepV, j - stepH));
            }
            else if (ch == 2) {
                std::get<2>(dst(i, j)) = std::get<0>(color(i - stepV, j - stepH));
            }
        }
    }
}

double gauss(double x, double y, double sigma) {
    return std::exp(-(x * x + y * y) / (2.0 * sigma * sigma)) / (M_PI * (2.0 * sigma * sigma));
}

/*********************************************************************************************************************************************************************/

Image align(Image src_image, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, bool isInterp, bool isSubpixel, double subScale) {
    uint h = src_image.n_rows / 3;
    uint w = src_image.n_cols;

    Image R = src_image.submatrix(2 * h, 0, h, w);
    Image G = src_image.submatrix(h, 0, h, w);
    Image B = src_image.submatrix(0, 0, h, w);

    int maxShift = 15;
    int stepV1, stepH1;
    int stepV2, stepH2;

    Image RGB(G);

    if (isSubpixel) {
        if (isInterp) {
            RGB = resizeBicubic(RGB, subScale);
            R = resizeBicubic(R, subScale);
            B = resizeBicubic(B, subScale);
        } else {
            RGB = resize(RGB, subScale);
            R = resize(R, subScale);
            B = resize(B, subScale);
        }

        MSE(RGB, R, B, maxShift * subScale, stepV1, stepH1, stepV2, stepH2);

        setChannel(RGB, R, 0, stepV1, stepH1);
        setChannel(RGB, B, 2, stepV2, stepH2);

        if (isInterp)
            RGB = resizeBicubic(RGB, 1.0 / subScale);
        else
            RGB = resize(RGB, 1.0 / subScale);

    } else {
        pyramid(RGB, R, B, maxShift, stepV1, stepH1, stepV2, stepH2, isInterp);

        setChannel(RGB, R, 0, stepV1, stepH1);
        setChannel(RGB, B, 2, stepV2, stepH2);
    }

    // Crop colored borders
    int x0, y0, w0, h0;

    if (stepV1 > 0 && stepV2 > 0) {
        y0 = std::max(stepV1, stepV2);
        h0 = h - y0;
    } else if (stepV1 < 0 && stepV2 < 0) {
        y0 = 0;
        h0 = h + std::min(stepV1, stepV2);
    } else {
        y0 = std::max(stepV1, stepV2);
        h0 = h + std::min(stepV1, stepV2) - y0;
    }

    if (stepH1 > 0 && stepH2 > 0) {
        x0 = std::max(stepH1, stepH2);
        w0 = w - x0;
    } else if (stepH1 < 0 && stepH2 < 0) {
        x0 = 0;
        w0 = w + std::min(stepH1, stepH2);
    } else {
        x0 = std::max(stepH1, stepH2);
        w0 = w + std::min(stepH1, stepH2) - x0;
    }

    RGB = RGB.submatrix(y0, x0, h0, w0); // crop image

    // Post processing image
    if (isPostprocessing) {
        needMirror = isMirror;

        if (postprocessingType == "--gray-world")
            return gray_world(RGB);

        if (postprocessingType == "--unsharp")
            return unsharp(RGB);

        if (postprocessingType == "--autocontrast")
            return autocontrast(RGB, fraction);
    }

    return RGB;
}

void pyramid(const Image &I, const Image &I1, const Image &I2, const int maxShift, int &minStepY1, int &minStepx1, int &minStepY2, int &minStepx2, bool isBicubic) {
    if (std::min(I.n_cols, I.n_rows) > 400) {
        if (isBicubic)
            pyramid(resizeBicubic(I, 0.5), resizeBicubic(I1, 0.5), resizeBicubic(I2, 0.5), maxShift, minStepY1, minStepx1, minStepY2, minStepx2, isBicubic);
        else
            pyramid(resize(I, 0.5), resize(I1, 0.5), resize(I2, 0.5), maxShift, minStepY1, minStepx1, minStepY2, minStepx2, isBicubic);

        MSE_big(I, I1, I2, 2, minStepY1 *= 2, minStepx1 *= 2, minStepY2 *= 2, minStepx2 *= 2);
    } else
        MSE(I, I1, I2, maxShift, minStepY1, minStepx1, minStepY2, minStepx2);
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{ -1, 0, 1},
        { -2, 0, 2},
        { -1, 0, 1}
    };

    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
        { 0,  0,  0},
        { -1, -2, -1}
    };

    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {
        { -1.0 / 6, -2.0 / 3, -1.0 / 6 },
        { -2.0 / 3, 13.0 / 3, -2.0 / 3 },
        { -1.0 / 6, -2.0 / 3, -1.0 / 6 }
    };

    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
    double Sr = 0;
    double Sg = 0;
    double Sb = 0;

    for (uint i = 0; i < src_image.n_rows; i++) {
        for (uint j = 0; j < src_image.n_cols; j++) {
            Sr += std::get<0>(src_image(i, j));
            Sg += std::get<1>(src_image(i, j));
            Sb += std::get<2>(src_image(i, j));
        }
    }

    double S = (Sr + Sg + Sb) / 3.0;
    double minV = 1.0 / 255;

    Sr = Sr < minV ? 0 : S / Sr;
    Sg = Sg < minV ? 0 : S / Sg;
    Sb = Sb < minV ? 0 : S / Sb;

    for (uint i = 0; i < src_image.n_rows; i++)
        for (uint j = 0; j < src_image.n_cols; j++)
            src_image(i, j) = Pixel(src_image(i, j)).scaleRGB(Sr, Sg, Sb).getRGB();

    return src_image;
}

Image resize(Image src_image, double scale) {
    // new size
    uint w0 = src_image.n_cols * scale;
    uint h0 = src_image.n_rows * scale;

    Image dst_image(h0, w0);

    for (uint i = 0; i < h0; i++) {
        for (uint j = 0; j < w0; j++) {
            uint y = std::min(uint(i / scale), src_image.n_rows - 2);
            uint x = std::min(uint(j / scale), src_image.n_cols - 2);

            double dy = (i / scale) - y;
            double dx = (j / scale) - x;

            double b[4];
            b[0] = (1 - dx) * (1 - dy);
            b[1] = dx * (1 - dy);
            b[2] = (1 - dx) * dy;
            b[3] = dx * dy;

            Pixel p[4];
            p[0] = Pixel(src_image(y, x));
            p[1] = Pixel(src_image(y, x + 1));
            p[2] = Pixel(src_image(y + 1, x));
            p[3] = Pixel(src_image(y + 1, x + 1));

            double rgb[3] = { 0, 0, 0 };

            for (uint k = 0; k < 4; k++) {
                rgb[0] += p[k].R() * b[k];
                rgb[1] += p[k].G() * b[k];
                rgb[2] += p[k].B() * b[k];
            }

            dst_image(i, j) = Pixel(rgb[0], rgb[1], rgb[2]).getRGB();
        }
    }

    return dst_image;
}

Image resizeBicubic(Image src_image, double scale) {
    uint w = src_image.n_cols * scale;
    uint h = src_image.n_rows * scale;

    Image dst_image(h, w);

    for (uint i = 0; i < h; i++) {
        for (uint j = 0; j < w; j++) {
            uint x = j / scale;
            uint y = i / scale;

            if (x > src_image.n_cols - 3)
                x = src_image.n_cols - 3;

            if (y > src_image.n_rows - 3)
                y = src_image.n_rows - 3;

            if (x < 1)
                x = 1;

            if (y < 1)
                y = 1;

            double dx = (j / scale) - x;
            double dy = (i / scale) - y;

            double b[16];
            b[0] = 1.0 / 4 * (dx - 1) * (dx - 2) * (dx + 1) * (dy - 1) * (dy - 2) * (dy + 1);
            b[1] = -1.0 / 4 * dx * (dx + 1) * (dx - 2) * (dy - 1) * (dy - 2) * (dy + 1);
            b[2] = -1.0 / 4 * dy * (dx - 1) * (dx - 2) * (dx + 1) * (dy + 1) * (dy - 2);
            b[3] = 1.0 / 4 * dx * dy * (dx + 1) * (dx - 2) * (dy + 1) * (dy - 2);
            b[4] = -1.0 / 12 * dx * (dx - 1) * (dx - 2) * (dy - 1) * (dy - 2) * (dy + 1);
            b[5] = -1.0 / 12 * dy * (dx - 1) * (dx - 2) * (dx + 1) * (dy - 1) * (dy - 2);
            b[6] = 1.0 / 12 * dx * dy * (dx - 1) * (dx - 2) * (dy + 1) * (dy - 2);
            b[7] = 1.0 / 12 * dx * dy * (dx + 1) * (dx - 2) * (dy - 1) * (dy - 2);
            b[8] = 1.0 / 12 * dx * (dx - 1) * (dx + 1) * (dy - 1) * (dy - 2) * (dy + 1);
            b[9] = 1.0 / 12 * dy * (dx - 1) * (dx - 2) * (dx + 1) * (dy - 1) * (dy + 1);
            b[10] = 1.0 / 36 * dx * dy * (dx - 1) * (dx - 2) * (dy - 1) * (dy - 2);
            b[11] = -1.0 / 12 * dx * dy * (dx - 1) * (dx + 1) * (dy + 1) * (dy - 2);
            b[12] = -1.0 / 12 * dx * dy * (dx + 1) * (dx - 2) * (dy - 1) * (dy + 1);
            b[13] = -1.0 / 36 * dx * dy * (dx - 1) * (dx + 1) * (dy - 1) * (dy - 2);
            b[14] = -1.0 / 36 * dx * dy * (dx - 1) * (dx - 2) * (dy - 1) * (dy + 1);
            b[15] = 1.0 / 36 * dx * dy * (dx - 1) * (dx + 1) * (dy - 1) * (dy + 1);

            Pixel p[16];
            p[0] = Pixel(src_image(y, x));
            p[1] = Pixel(src_image(y, x + 1));
            p[2] = Pixel(src_image(y + 1, x));
            p[3] = Pixel(src_image(y + 1, x + 1));
            p[4] = Pixel(src_image(y, x - 1));
            p[5] = Pixel(src_image(y - 1, x));
            p[6] = Pixel(src_image(y + 1, x - 1));
            p[7] = Pixel(src_image(y - 1, x + 1));
            p[8] = Pixel(src_image(y, x + 2));
            p[9] = Pixel(src_image(y + 2, x));
            p[10] = Pixel(src_image(y - 1, x - 1));
            p[11] = Pixel(src_image(y + 1, x + 2));
            p[12] = Pixel(src_image(y + 2, x + 1));
            p[13] = Pixel(src_image(y - 1, x + 2));
            p[14] = Pixel(src_image(y + 2, x - 1));
            p[15] = Pixel(src_image(y + 2, x + 2));

            double rgb[3] = { 0, 0, 0 };

            for (uint k = 0; k < 16; k++) {
                rgb[0] += p[k].R() * b[k];
                rgb[1] += p[k].G() * b[k];
                rgb[2] += p[k].B() * b[k];
            }

            dst_image(i, j) = Pixel(rgb[0], rgb[1], rgb[2]).getRGB();
        }
    }

    return dst_image;
}

Image mirror(const Image &src_image, uint Rh, uint Rw) {
    uint w = src_image.n_cols + 2 * Rw;
    uint h = src_image.n_rows + 2 * Rh;

    Image dst_image(h, w);

    for (uint i = 0; i < h; i++) {
        for (uint j = 0; j < w; j++) {
            uint x = j, y = i;

            if (i < Rh || i > h - Rh - 1 || j < Rw || j > w - Rw - 1) {
                if (i < Rh)
                    y = 2 * Rh - 1 - i;
                else if (i >= h - Rh)
                    y = 2 * (h - Rh) - 1 - i;

                if (j < Rw)
                    x = 2 * Rw - 1 - j;
                else if (j >= w - Rw)
                    x = 2 * (w - Rw) - 1 - j;
            }

            dst_image(i, j) = src_image(y - Rh, x - Rw);
        }
    }

    return dst_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.

    uint Rw = (kernel.n_cols - 1) / 2; // радиус ядра по горизонтали
    uint Rh = (kernel.n_rows - 1) / 2; // радиус ядра по вертикали

    if (needMirror)
        src_image = mirror(src_image, Rh, Rw);

    Image dst_image(src_image.n_rows - 2 * Rh, src_image.n_cols - 2 * Rw);

    for (uint y = Rh; y < src_image.n_rows - Rh; y++) {
        for (uint x = Rw; x < src_image.n_cols - Rw; x++) {
            double r = 0, g = 0, b = 0;

            for (uint i = 0; i < kernel.n_rows; i++) {
                for (uint j = 0; j < kernel.n_cols; j++) {
                    r += std::get<0>(src_image(y - Rh + i, x - Rw + j)) * kernel(i, j);
                    g += std::get<1>(src_image(y - Rh + i, x - Rw + j)) * kernel(i, j);
                    b += std::get<2>(src_image(y - Rh + i, x - Rw + j)) * kernel(i, j);
                }
            }

            dst_image(y - Rh, x - Rw) = Pixel(r, g, b).getRGB();
        }
    }

    return dst_image;
}

Image autocontrast(Image src_image, double fraction) {
    uint w = src_image.n_cols;
    uint h = src_image.n_rows;

    uint hist[256];

    for (uint i = 0; i < 256; i++)
        hist[i] = 0;

    for (uint i = 0; i < h; i++)
        for (uint j = 0; j < w; j++)
            hist[Pixel(src_image(i, j)).getBrightness()]++;

    uint skipped = w * h * fraction;
    uint min_br = 0;
    uint max_br = 255;

    if (skipped > 0) {
        uint left = 0, right = 0;
        bool find = false;

        while (min_br < 255 && max_br > 0 && !find) {
            find = true;

            if (left + hist[min_br] < skipped) {
                left += hist[min_br++];
                find = false;
            }

            if (right + hist[max_br] < skipped) {
                right += hist[max_br--];
                find = false;
            }
        }
    } else {
        while (min_br < 255 && hist[min_br] == 0)
            min_br++;

        while (max_br > 0 && hist[max_br] == 0)
            max_br--;
    }

    double k = 255.0 / (max_br - min_br);

    Image dst_image(h, w);

    for (uint i = 0; i < h; i++)
        for (uint j = 0; j < w; j++)
            dst_image(i, j) = ((Pixel(src_image(i, j)) - min_br) * k).getRGB();

    return dst_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    int size = 2 * radius + 1;
    double sum = 0;

    Matrix<double>kernel(size, size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double value = gauss(i - radius, j - radius, sigma);

            kernel(i, j) = value;
            sum += value;
        }
    }

    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            kernel(i, j) /= sum;

    return custom(src_image, kernel);
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    int size = 2 * radius + 1;

    Matrix<double>kernel_h(1, size);
    Matrix<double>kernel_v(size, 1);

    double sum = 0;

    for (int i = 0; i < size; i++) {
        double value = gauss(0, i - radius, sigma);

        kernel_h(0, i) = value;
        kernel_v(i, 0) = value;

        sum += value;
    }

    for (int i = 0; i < size; i++) {
        kernel_h(0, i) /= sum;
        kernel_v(i, 0) /= sum;
    }

    return custom(custom(src_image, kernel_h), kernel_v);
}

Image median(Image src_image, int radius) {
    src_image = mirror(src_image, radius, radius);

    Image dst_image(src_image.n_rows - 2 * radius, src_image.n_cols - 2 * radius);

    for (uint y = radius; y < src_image.n_rows - radius; y++)
        for (uint x = radius; x < src_image.n_cols - radius; x++)
            dst_image(y - radius, x - radius) = calcMedian(src_image, y, x, radius, false).getRGB();

    return dst_image;
}

Image median_linear(Image src_image, int radius) {
    src_image = mirror(src_image, radius, radius);

    Image dst_image(src_image.n_rows - 2 * radius, src_image.n_cols - 2 * radius);

    for (uint y = radius; y < src_image.n_rows - radius; y++)
        for (uint x = radius; x < src_image.n_cols - radius; x++)
            dst_image(y - radius, x - radius) = calcMedian(src_image, y, x, radius, true).getRGB();

    return dst_image;
}

Image median_const(Image src_image, int radius) {
    src_image = mirror(src_image, radius, radius);

    uint width = src_image.n_cols;
    uint height = src_image.n_rows;

    Image dst_image(height - 2 * radius, width - 2 * radius);

    ConstantMedian med(src_image, radius);

    for (uint j = radius; j < width - radius - 1; j++) {
        dst_image(0, j - radius) = med.getMedian().getRGB();
        med.shiftRight();
    }

    dst_image(0, width - radius - 1 - radius) = med.getMedian().getRGB();

    for (uint i = radius; i < height - radius; i++) {
        med.shiftDown(src_image.submatrix(i + radius, 0, 1, 2 * radius + 1));

        for (uint j = radius; j < width - radius - 1; j++) {
            dst_image(i - radius, j - radius) = med.getMedian().getRGB();

            med.shiftRight(src_image(i + radius, j + radius + 1));
        }

        dst_image(i - radius, width - 2 * radius - 1) = med.getMedian().getRGB();
    }

    return dst_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
