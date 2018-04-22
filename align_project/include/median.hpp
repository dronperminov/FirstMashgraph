uint calcMedianSquare(std::vector<uint> &vector) {
    uint size = vector.size();
    uint index = size / 2;

    for (uint i = 0; i <= index; i++) {
        uint min = vector[i];
        uint imin = i;

        for (uint j = i; j < size; j++)
            if (min > vector[j]) {
                min = vector[j];
                imin = j;
            }

        vector[imin] = vector[i];
        vector[i] = min;
    }

    return vector[index];
}

uint calcMedianLinear(std::vector<uint> &vector) {
    uint left = 0;
    uint right = vector.size() - 1;
    uint k = vector.size() / 2;
    uint i, j;

    while (left < right) {
        uint x = vector[k];
        i = left;
        j = right;

        do {
            while (vector[i] < x)
                i++;
            
            while (x < vector[j])
                j--;

            if (i <= j) {
                uint tmp = vector[i];
                vector[i] = vector[j];
                vector[j] = tmp;
                i++;
                j--;
            }
        } while (i <= j);

        if (j < k)
            left = i;

        if (k < i)
            right = j;
    }

    return vector[k];
}

Pixel calcMedian(const Image &I, uint y, uint x, int radius, bool linear = false) {
    std::vector<uint> red;
    std::vector<uint> green;
    std::vector<uint> blue;

    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            red.push_back(std::get<0>(I(y + i, x + j)));
            green.push_back(std::get<1>(I(y + i, x + j)));
            blue.push_back(std::get<2>(I(y + i, x + j)));
        }
    }

    if (linear)
    	return Pixel(calcMedianLinear(red), calcMedianLinear(green), calcMedianLinear(blue));

    return Pixel(calcMedianSquare(red), calcMedianSquare(green), calcMedianSquare(blue));
}

RGBhistogram::RGBhistogram(const Image &I) : env(std::deque<Pixel>()), width(0), height(0) {
    height = I.n_rows;
    width = I.n_cols;

    for (uint i = 0; i < 256; i++)
        rgb[0][i] = rgb[1][i] = rgb[2][i] = 0;

    for (uint i = 0; i < height; i++) {
        for (uint j = 0; j < width; j++) {
        	Pixel p = Pixel(I(i, j));

            env.push_back(p);

            rgb[0][std::get<0>(I(i, j))]++;
            rgb[1][std::get<1>(I(i, j))]++;
            rgb[2][std::get<2>(I(i, j))]++;
        }
    }
}

ConstantMedian::ConstantMedian(const Image &I, uint radius_) : columns(std::vector<RGBhistogram>()), pos(2 * radius_ + 1), radius(radius_) {
    for (uint i = 0; i < 256; i++) {
        rgb[0][i] = 0;
        rgb[1][i] = 0;
        rgb[2][i] = 0;
    }

    for (uint j = 0; j < I.n_cols; j++) {
        columns.push_back(RGBhistogram(I.submatrix(0, j, 2 * radius + 1, 1)));

        if (j < 2 * radius + 1) {
            for (uint i = 0; i < 256 ; i++)
                for (uint ch = 0; ch < 3; ch++)
                    rgb[ch][i] += columns.back().rgb[ch][i];
        }
    }

    for (uint i = 0; i < 256; i++)
        for (uint ch = 0; ch < 3; ch++)
            old[ch][i] = rgb[ch][i];
}

void ConstantMedian::shiftRight() {
    for (uint i = 0; i < 256; i++)
        for (uint ch = 0; ch < 3; ch++)
            rgb[ch][i] += columns[pos].rgb[ch][i] - columns[pos - 2 * radius + 1].rgb[ch][i];

    pos++;
}

void ConstantMedian::shiftRight(const Pixel &p) {
    Pixel front = columns[pos].env.front();

    columns[pos].rgb[0][front.R()]--;
    columns[pos].rgb[1][front.G()]--;
    columns[pos].rgb[2][front.B()]--;

    columns[pos].env.pop_front();

    columns[pos].rgb[0][p.R()]++;
    columns[pos].rgb[1][p.G()]++;
    columns[pos].rgb[2][p.B()]++;

    columns[pos].env.push_back(p);

    for (uint i = 0; i < 256; i++)
        for (uint ch = 0; ch < 3; ch++)
            rgb[ch][i] += columns[pos].rgb[ch][i] - columns[pos - 2 * radius - 1].rgb[ch][i];

    pos++;
}

void ConstantMedian::shiftDown(const Image &row) {
    for (int i = 0; i < 256; i++)
        for (uint ch = 0; ch < 3; ch++)
            rgb[ch][i] = old[ch][i];

    for (uint j = 0; j < 2 * radius + 1; j++) {
        Pixel p = columns[j].env.front();

        columns[j].rgb[0][p.R()]--;
        columns[j].rgb[1][p.G()]--;
        columns[j].rgb[2][p.B()]--;

        columns[j].env.pop_front();

        columns[j].rgb[0][std::get<0>(row(0, j))]++;
        columns[j].rgb[1][std::get<1>(row(0, j))]++;
        columns[j].rgb[2][std::get<2>(row(0, j))]++;

        columns[j].env.push_back(row(0, j));

        rgb[0][p.R()]--;
        rgb[1][p.G()]--;
        rgb[2][p.B()]--;

        rgb[0][std::get<0>(row(0, j))]++;
        rgb[1][std::get<1>(row(0, j))]++;
        rgb[2][std::get<2>(row(0, j))]++;
    }

    pos = 2 * radius + 1;

    for (int i = 0; i < 256; i++)
        for (uint ch = 0; ch < 3; ch++)
            old[ch][i] = rgb[ch][i];        
}

Pixel ConstantMedian::getMedian() {
    const double median_index = 2 * radius * (1 + radius); // (2r + 1)(2r + 1)) / 2

    uint median[3] = { 0, 0, 0 };

    for (uint ch = 0; ch < 3; ch++) {
        long long sum = 0;
        uint i = 0;

        while (i < 256 && sum < median_index)
            sum += rgb[ch][i++];

        median[ch] = i - 1;
    }

    return Pixel(median[0], median[1], median[2]);
}