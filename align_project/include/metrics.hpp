void shiftMSE(const Image &I, const Image &I1, const Image &I2, int stepY, int stepX, long long &v1, long long &v2) {
    v1 = 0;
    v2 = 0;

    const uint border_x = 15;
    const uint border_y = 15;

    uint i0 = border_y + (stepY > 0 ? stepY : 0);
    uint i1 = I.n_rows + (stepY > 0 ? 0 : stepY) - border_y;

    uint j0 = border_x + (stepX > 0 ? stepX : 0);
    uint j1 = I.n_cols + (stepX > 0 ? 0 : stepX) - border_x;

    for (uint i = i0; i < i1; i++) {
        for (uint j = j0; j < j1; j++) {
            long v = long(std::get<0>(I(i, j))) - long(std::get<0>(I1(i - stepY, j - stepX)));
            v1 += v * v;

            v = long(std::get<0>(I(i, j))) - long(std::get<0>(I2(i - stepY, j - stepX)));
            v2 += v * v;
        }
    }
}

void MSE(const Image &I1, const Image &I2, const Image &I3, const int maxShift, int &minStepY1, int &minStepX1, int &minStepY2, int &minStepX2) {
    minStepY1 = 0;
    minStepX1 = 0;    
    minStepY2 = 0;
    minStepX2 = 0;

    long long min1, min2;
    long long m1, m2;
    shiftMSE(I1, I2, I3, 0, 0, min1, min2);

    for (int stepY = -maxShift; stepY <= maxShift; stepY++) {
        for (int stepX = -maxShift; stepX <= maxShift; stepX++) {
            shiftMSE(I1, I2, I3, stepY, stepX, m1, m2);

            if (m1 < min1) {
                min1 = m1;
                minStepY1 = stepY;
                minStepX1 = stepX;
            }

            if (m2 < min2) {
                min2 = m2;
                minStepY2 = stepY;
                minStepX2 = stepX;
            }
        }
    }
}

long long shiftMSE(const Image &I1, const Image &I2, int stepY, int stepX) {
    uint width = I1.n_cols;
    uint height = I1.n_rows;
    long long sum = 0;

    uint border_x = width / 32;
    uint border_y = height / 32;

    uint start_i = stepY > 0 ? stepY : 0;
    uint end_i = height + (stepY > 0 ? 0 : stepY);

    uint start_j = stepX > 0 ? stepX : 0;
    uint end_j = width + (stepX > 0 ? 0 : stepX);

    for (uint i = border_y + start_i; i < end_i - border_y; i++) {
        for (uint j = border_x + start_j; j < end_j - border_x; j++) {
            long pixelValue = std::get<0>(I1(i, j));
            pixelValue -= std::get<0>(I2(i - stepY, j - stepX));
            sum += pixelValue * pixelValue;
        }
    }

    return sum;
}

void MSE_big(const Image &I1, const Image &I2, const Image &I3, const int maxShift, int &minStepY1, int &minStepX1, int &minStepY2, int &minStepX2) {
    long long m1, m2;

    long long min1 = shiftMSE(I1, I2, minStepY1, minStepX1);
    long long min2 = shiftMSE(I1, I3, minStepY2, minStepX2);

    int v1 = minStepY1, h1 = minStepX1;
    int v2 = minStepY2, h2 = minStepX2;

    for (int stepY = -maxShift; stepY <= maxShift; stepY++) {
        for (int stepX = -maxShift; stepX <= maxShift; stepX++) {
            m1 = shiftMSE(I1, I2, v1 + stepY, h1 + stepX);

            if (m1 < min1) {
                min1 = m1;
                minStepY1 = v1 + stepY;
                minStepX1 = h1  +stepX;
            }
            
            m2 = shiftMSE(I1, I3, v2 + stepY, h2 + stepX);

            if (m2 < min2) {
                min2 = m2;
                minStepY2 = v2 + stepY;
                minStepX2 = h2 + stepX;
            }
        }
    }
}