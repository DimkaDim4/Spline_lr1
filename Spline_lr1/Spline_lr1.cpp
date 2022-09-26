#include <iostream>
#include <vector>
#include <fstream>
#include <random>

struct Point2D {
    double x;
    double y;
    double w;

    Point2D() {
        x = 0.;
        y = 0.;
        w = 1.;
    }

    Point2D(double a, double b, double c) {
        x = a; y = b; w = c;
    }
};

// ksi = [0, 1]
double phi_1(const double &ksi) {
    return 1. - 3. * ksi * ksi + 2. * ksi * ksi * ksi;
}

// ksi = [0, 1]
double phi_2(const double &ksi) {
    return ksi - 2. * ksi * ksi + ksi * ksi * ksi;
}

// ksi = [0, 1]
double phi_3(const double &ksi) {
    return 3. * ksi * ksi - 2. * ksi * ksi * ksi;
}

// ksi = [0, 1]
double phi_4(const double &ksi) {
    return - ksi * ksi + ksi * ksi * ksi;
}

double Spline(const std::vector<Point2D>& input, const double* x_k, const double* h_k, const int M, const double &alpha, double* solve);

int main()
{
    bool Filter = false;
    bool Filter_coeff = 1.5;
    int M = 20000; // количество точек, задавающие конечно-элементную сетку
    int N = M - 1; // количество конечных элементов
    double alpha = .1;

    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution<> distribution(-0.1, 0.1);

    std::ofstream out;
    out.open("Points.txt");
    int count_point = 40000;
    double a0 = 0.;
    double b0 = 10.;
    double h0 = (b0 - a0) / (count_point - 1);

    std::vector<Point2D> input;
    double inp_x, inp_y;
    if (!out.is_open()) { return -1; }
    for (int i = 0; i < count_point; ++i) {
        double r = distribution(generator);
        //if ((i > 19) && (i < 61)) {continue;}

        double rr = (double)rand() / RAND_MAX;

        inp_x = a0 + i * h0;
        //inp_y = sin(inp_x);
        inp_y = sin(inp_x) + r;
        //inp_y = r;

        //if (i == 100) inp_y += 5;
        if (rr > 0.9) inp_y += 2.;

        input.push_back(Point2D(inp_x, inp_y, 1.));
        out << inp_x << "\t" << inp_y << "\n";
    }
    out.close();

    std::vector<std::vector<Point2D>> points;
    points.resize(N);
    double* x_k = new double[M];
    double* h_k = new double[N];
    double* solve = new double[2 * M];

    double hk = (b0 - a0) / N;
    for (int i = 0; i < M; ++i) {
        x_k[i] = a0 + i * hk;
    }
    for (int i = 0; i < N; ++i) {
        h_k[i] = x_k[i + 1] - x_k[i];
    }

    //x_k[5] = x_k[6];
    //h_k[5] += h_k[6];
    //h_k[6] = 0.;

    double disp;
    bool re_compute;
    do {
        re_compute = false;
        disp = Spline(input, x_k, h_k, M, alpha, solve);
        //std::cout << disp << "\n";
        double disp_point = 0.;
        for (int j = 0; j < input.size(); j++) {
            for (int i = 0; i < N; i++) {
                double x = input[j].x;
                double y = input[j].y;
                double w = input[j].w;
                bool is_belongs = i == N - 1 ? ((x >= x_k[i]) && (x <= x_k[i + 1])) : ((x >= x_k[i]) && (x < x_k[i + 1]));
                if (is_belongs) {
                    double ksi = (x - x_k[i]) / h_k[i];
                    double phi1 = phi_1(ksi);
                    double phi2 = h_k[i] * phi_2(ksi);
                    double phi3 = phi_3(ksi);
                    double phi4 = h_k[i] * phi_4(ksi);
                    double _y = solve[2 * i] * phi1 + solve[2 * i + 1] * phi2 + solve[2 * i + 2] * phi3 + solve[2 * i + 3] * phi4;
                    disp_point = abs(_y - y) * w;
                    if (Filter && (disp_point >= Filter_coeff * disp)) {
                        input[j].w *= 0.5;
                        re_compute = true;
                    }
                    break;
                }
            }
        }
    } while (re_compute);

    out.open("Solution.txt");
    if (!out.is_open()) { return -1; } 
    for (double _x = a0; _x <= b0; _x += 0.01) {
        for (int i = 0; i < N; i++) {
            bool is_belongs = i == N - 1 ? ((_x >= x_k[i]) && (_x <= x_k[i + 1])) : ((_x >= x_k[i]) && (_x < x_k[i + 1]));
            //if ((is_belongs) && (points.at(i).size() != 0.)) {
            if (is_belongs) {
                double ksi = (_x - x_k[i]) / h_k[i];
                double phi1 = phi_1(ksi);
                double phi2 = h_k[i] * phi_2(ksi);
                double phi3 = phi_3(ksi);
                double phi4 = h_k[i] * phi_4(ksi);
                double _y = solve[2 * i] * phi1 + solve[2 * i + 1] * phi2 + solve[2 * i + 2] * phi3 + solve[2 * i + 3] * phi4;

                out << _x << "\t" << _y << "\n";
                break;
            }
        }
    }
    out.close();

    delete[] x_k;
    delete[] h_k;
    delete[] solve;
}

double Spline(const std::vector<Point2D>& input, const double* x_k, const double* h_k,const int M, const double& alpha, double* solve) {
    const int N = M - 1;
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution<> distribution(-0.1, 0.1);

    double a0 = input.at(0).x;
    double b0 = input.at(input.size() - 1.).x;
    double h0 = (b0 - a0) / (input.size() - 1.);

    std::vector<std::vector<Point2D>> points;
    points.resize(N);

    double* f = new double[2 * M];
    double* a = new double[2 * M];

    double* b_1 = new double[2 * M - 1];
    double* b_2 = new double[2 * M - 1];

    double* c_1 = new double[2 * M - 2];
    double* c_2 = new double[2 * M - 2];

    double* d_1 = new double[2 * M - 3];
    double* d_2 = new double[2 * M - 3];


    for (int i = 0; i < 2 * M; ++i) {
        a[i] = f[i] = 0.;
    }

    for (int i = 0; i < 2 * M - 1; ++i) {
        b_2[i] = b_1[i] = 0.;
    }

    for (int i = 0; i < 2 * M - 2; ++i) {
        c_2[i] = c_1[i] = 0.;
    }

    for (int i = 0; i < 2 * M - 3; i += 2) {
        d_2[i] = d_1[i] = 0.;
    }

    auto it = input.begin();
    for (int i = 0; i < input.size(); i++, it++) {
        for (int j = 0; j < N; j++) {
            bool is_belongs = j == N - 1 ? ((it->x >= x_k[j]) && (it->x <= x_k[j + 1])) : ((it->x >= x_k[j]) && (it->x < x_k[j + 1]));
            if (is_belongs) {
                points.at(j).push_back(*it);
                break;
            }
        }
    }

    // диагонали глобальной матрицы
    // вид матрицы для трех конечных элементов
    // a0  b01 c01 d01 
    // b02 a1  b11 c11
    // c02 b12 a2  b21 c21 d21
    // d02 c12 b22 a3  b31 c31
    //         c22 b32 a4  b41 c41 d41
    //         d22 c32 b42 a5  b51 c51
    //                 c42 b52 a6  b61
    //                 d42 c52 b62 a7
    //

    std::vector<Point2D>* points_in_element;
    Point2D point;
    double h_el;
    double x_el;
    for (int i = 0; i < 2 * N; i += 2) {
        points_in_element = &points[i / 2];
        h_el = h_k[i / 2];
        x_el = x_k[i / 2];
        for (auto j = points_in_element->begin(); j != points_in_element->end(); j++) {
            point = *j;
            double ksi = (point.x - x_el) / h_el;
            double phi1 = phi_1(ksi);
            double phi2 = h_el * phi_2(ksi);
            double phi3 = phi_3(ksi);
            double phi4 = h_el * phi_4(ksi);

            f[i] += phi1 * point.y * point.w;
            f[i + 1] += phi2 * point.y * point.w;
            f[i + 2] += phi3 * point.y * point.w;
            f[i + 3] += phi4 * point.y * point.w;

            // симметричная матрица
            a[i] += phi1 * phi1 * point.w;
            a[i + 1] += phi2 * phi2 * point.w;
            a[i + 2] += phi3 * phi3 * point.w;
            a[i + 3] += phi4 * phi4 * point.w;

            b_1[i] += phi1 * phi2 * point.w;
            b_1[i + 1] += phi2 * phi3 * point.w;
            b_1[i + 2] += phi3 * phi4 * point.w;

            c_1[i] += phi1 * phi3 * point.w;
            c_1[i + 1] += phi2 * phi4 * point.w;

            d_1[i] += phi1 * phi4 * point.w;
        }

        a[i] += alpha / (1.2 * h_el);
        a[i + 1] += h_el * alpha / 7.5;
        a[i + 2] += alpha / (1.2 * h_el);
        a[i + 3] += h_el * alpha / 7.5;

        b_1[i] += alpha / 10.;
        b_1[i + 1] -= alpha / 10.;
        b_1[i + 2] -= alpha / 10.;

        c_1[i] -= alpha / (1.2 * h_el);
        c_1[i + 1] -= h_el * alpha / 30.;

        d_1[i] += alpha / 10.;
    }

    for (int i = 0; i < 2 * M - 1; ++i) {
        b_2[i] = b_1[i];
    }

    for (int i = 0; i < 2 * M - 2; ++i) {
        c_2[i] = c_1[i];
    }

    for (int i = 0; i < 2 * M - 3; i += 2) {
        d_2[i] = d_1[i];
    }

    // решаем полученную СЛАУ
    // прямой шаг Гаусса
    double coeff;
    double coeff1;
    double vec1[4];

    bool is_compute1;
    bool is_compute2;
    bool is_compute3;
    bool is_compute4;
    bool is_compute5;
    for (int i = 0; i < 2 * N; i += 2) {
        if (a[i] != 0.) {
            coeff1 = a[i];
            vec1[0] = a[i];
            vec1[1] = b_1[i];
            vec1[2] = c_1[i];
            vec1[3] = d_1[i];

            is_compute1 = b_2[i] != 0. ? true : false;
            is_compute2 = c_2[i] != 0. ? true : false;
            is_compute3 = d_2[i] != 0. ? true : false;
        }
        else {
            is_compute1 = false;
            is_compute2 = false;
            is_compute3 = false;
        }


        if (is_compute1) {
            coeff = b_2[i] / coeff1;
            b_2[i] -= vec1[0] * coeff; // = 0
            a[i + 1] -= vec1[1] * coeff;
            b_1[i + 1] -= vec1[2] * coeff;
            c_1[i + 1] -= vec1[3] * coeff;
            f[i + 1] -= f[i] * coeff;
        }

        if (is_compute2) {
            coeff = c_2[i] / coeff1;
            c_2[i] -= vec1[0] * coeff; // = 0
            b_2[i + 1] -= vec1[1] * coeff;
            a[i + 2] -= vec1[2] * coeff;
            b_1[i + 2] -= vec1[3] * coeff;
            f[i + 2] -= f[i] * coeff;
        }

        if (is_compute3) {
            coeff = d_2[i] / coeff1;
            d_2[i] -= vec1[0] * coeff; // = 0
            c_2[i + 1] -= vec1[1] * coeff;
            b_2[i + 2] -= vec1[2] * coeff;
            a[i + 3] -= vec1[3] * coeff;
            f[i + 3] -= f[i] * coeff;
        }

        if (a[i + 1] != 0.) {
            coeff1 = a[i + 1];
            vec1[0] = a[i + 1];
            vec1[1] = b_1[i + 1];
            vec1[2] = c_1[i + 1];

            is_compute4 = b_2[i + 1] != 0. ? true : false;
            is_compute5 = c_2[i + 1] != 0. ? true : false;
        }
        else {
            is_compute4 = false;
            is_compute5 = false;
        }


        if (is_compute4) {
            coeff = b_2[i + 1] / coeff1;
            b_2[i + 1] -= vec1[0] * coeff; // = 0
            a[i + 2] -= vec1[1] * coeff;
            b_1[i + 2] -= vec1[2] * coeff;
            f[i + 2] -= f[i + 1] * coeff;
        }

        if (is_compute5) {
            coeff = c_2[i + 1] / coeff1;
            c_2[i + 1] -= vec1[0] * coeff; // = 0
            b_2[i + 2] -= vec1[1] * coeff;
            a[i + 3] -= vec1[2] * coeff;
            f[i + 3] -= f[i + 1] * coeff;
        }
    }

    if (a[2 * N] != 0.) {
        coeff = b_2[2 * N] / a[2 * N];
        b_2[2 * N] -= a[2 * N] * coeff; // = 0
        a[2 * N + 1] -= b_1[2 * N] * coeff;
        f[2 * N + 1] -= f[2 * N] * coeff;
    }

    // обратный шаг Гаусса
    solve[2 * M - 1] = a[2 * M - 1] != 0. ? f[2 * M - 1] / a[2 * M - 1] : 0.;
    solve[2 * M - 2] = a[2 * M - 2] != 0. ? (f[2 * M - 2] - b_1[2 * M - 2] * solve[2 * M - 1]) / a[2 * M - 2] : 0.;

    for (int i = 2 * M - 3; i > 0; --i) {
        solve[i] = a[i] != 0. ? (f[i] - b_1[i] * solve[i + 1] - c_1[i] * solve[i + 2]) / a[i] : 0.;
        --i;
        solve[i] = a[i] != 0. ? (f[i] - b_1[i] * solve[i + 1] - c_1[i] * solve[i + 2] - d_1[i] * solve[i + 3]) / a[i] : 0.;
    }

    double disp = 0.;
    for (int i = 0; i < input.size(); i++) {
        for (int i = 0; i < N; i++) {
            double x = input[i].x;
            double y = input[i].y;
            bool is_belongs = i == N - 1 ? ((x >= x_k[i]) && (x <= x_k[i + 1])) : ((x >= x_k[i]) && (x < x_k[i + 1]));
            if (is_belongs) {
                double ksi = (x - x_k[i]) / h_k[i];
                double phi1 = phi_1(ksi);
                double phi2 = h_k[i] * phi_2(ksi);
                double phi3 = phi_3(ksi);
                double phi4 = h_k[i] * phi_4(ksi);
                double _y = solve[2 * i] * phi1 + solve[2 * i + 1] * phi2 + solve[2 * i + 2] * phi3 + solve[2 * i + 3] * phi4;
                disp += abs(_y - y);
                break;
            }
        }
    }

    disp /= (input.size() - 1);

    delete[] f;
    delete[] a;

    delete[] b_1;
    delete[] b_2;

    delete[] c_1;
    delete[] c_2;

    delete[] d_1;
    delete[] d_2;

    return disp;
}
