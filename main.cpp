/* Note that you need the GNUPlot.
 *
 * macOS
 *
 * The easiest way to install it thought the Homebrew.
 * If you are not familiar with homebrew, read more about it here: https://brew.sh/
 * To install GNUPlot:
 * brew install gnuplot
 *
 * Linux
 *
 * You know how this works, don't you?
 */

#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <memory>
#include <stdexcept>
//#include <algorithm>


const unsigned N = 1.0e2; //Number of points.
constexpr double eps = 1.0 / N;
const double R = 1;
const double pi = 3.14159265359;

const double E_min = 1.0e5;
const double E_0 = 2.0e6;
const double E_e = 0.5110034e6;


typedef std::tuple<double, double, double> coord;

typedef std::tuple<double, double, double, double> longDoubleTuple;

const std::vector<longDoubleTuple> planes = {std::make_tuple(0, 0, 1, -1), //There's a plane equation presented like Ax+By+Cz+D=0.
                                             std::make_tuple(0, 1, 0, -1),
                                             std::make_tuple(1, 0, 0, 1),
                                             std::make_tuple(0, 1, 0, 1),
                                             std::make_tuple(1, 0, 0, -1)};

const std::tuple<double, double, double> Sigma_air_2 = std::make_tuple(0.0438, 0.0000003, 0.000387);
const double Sigma_air_sum = std::get<0>(Sigma_air_2) + std::get<1>(Sigma_air_2) + std::get<2>(Sigma_air_2);
const std::tuple<double, double, double> Sigma_Pb_2 = std::make_tuple(0.0349, 0.00523, 0.005);
const double Sigma_Pb_sum = std::get<0>(Sigma_Pb_2) + std::get<1>(Sigma_Pb_2) + std::get<2>(Sigma_Pb_2);

std::array<std::vector<coord>, N> interaction_points;

std::vector <coord> polar ();

std::vector <coord> coordinate_transformation(std::vector<coord>& coords);

void data_file_creation (std::string& DataType, std::vector<coord>& xx);

void default_distribution_plot(std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title);

longDoubleTuple beam_direction(double sigma);

coord coordinates_of_the_interaction(longDoubleTuple& beams);

std::string exec(std::string str_obj);

std::tuple<double, double, double> statistical_weight(std::tuple<double, double, double> sigma);

unsigned interaction_type(std::tuple<double, double, double>& Sigma);

std::array<std::vector<double>, N> interactions (std::vector<coord>& points);

double cost (coord& A, coord& B, coord& C);

void plot(std::array<std::vector<coord>, N>& points);

int main() {
    srand(time(nullptr));
    std::vector<coord> borning = std::move(polar());
    std::vector<coord> born_points = std::move(coordinate_transformation(borning));
    std::string name1 = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name1, born_points);
    default_distribution_plot(name1, name1, "x", "y", name1);
    std::array<std::vector<double>, N> Energies = interactions(born_points);
    plot(interaction_points);
    return 0;
}

std::vector <coord> polar () {
    std::vector <coord> coordinates;
    for (unsigned i = 0; i < N; i++) {
        double rho = R * std::sqrt(eps * (rand() % (N + 1)));
        double phi = 2 * pi * eps * (rand() % (N + 1));
        double z = 0;
        coordinates.emplace_back(std::move(std::make_tuple(phi, rho, z)));
    }
    return coordinates;
}

std::vector<coord> coordinate_transformation(std::vector<coord>& coords) {
    std::vector<coord> xOy;
    for(unsigned i = 0; i < coords.size(); i++) {
        double phi = std::get<0>(coords[i]);
        double rho = std::get<1>(coords[i]);
        double x = rho * cos(phi);
        double y = rho * sin(phi);
        xOy.emplace_back(std::move(std::make_tuple(x, y, std::get<2>(coords[i]))));
    }
    return xOy;
}

void data_file_creation (std::string& DataType, std::vector<coord>& xx) {
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << std::get<0>(xx[i]) << '\t' << std::get<1>(xx[i]) << '\t'<< std::get<2>(xx[i]) << '\t' << std::endl;
    fout.close();
}

void default_distribution_plot(std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title) {
    //if you have problems with ".svg", you have to change ".svg" to ".pdf" in strings bellow.
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term svg",
                                      "set out \'" + name + ".svg\'",
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set grid xtics ytics",
                                      "set title \'" + title + "\'",
                                      "plot \'" + data + "\' using 1:2 lw 1 lt rgb 'orange' ti \'Nodes\'",
                                      "set key box top right",
                                      "set terminal pop",
                                      "set output",
                                      "replot"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

longDoubleTuple beam_direction(double sigma) {
    double mu = 2 * eps * (rand() % (N + 1)) - 1; //cos(\phi);
    double L, a, b, d = 10;
    do {
        a = 2 * eps * (rand() % (N + 1)) - 1;
        b = 2 * eps * (rand() % (N + 1)) - 1;
        d = std::pow(a, 2) + std::pow(b, 2);
        L = - log(eps * (rand() % (N + 1))) / sigma;
    } while (d > 1 || std::isfinite(L) == 0);
    double cosPsi = a / std::sqrt(d);
    double sinPsi = b / std::sqrt(d);
    return std::make_tuple(mu, cosPsi, sinPsi, L);
}

coord coordinates_of_the_interaction (longDoubleTuple& beam) {
    double x = std::get<2>(beam) * std::get<3>(beam);
    double y = std::get<1>(beam) * std::get<3>(beam);
    double z = std::get<0>(beam) * std::get<3>(beam);
    return std::make_tuple(x, y, std::abs(z));
}

void cap() {
    std::vector<coord> cage = {std::make_tuple(-1, 1, 1),
                               std::make_tuple(1, 1, 1),
                               std::make_tuple(1, -1, 1),
                               std::make_tuple(-1, -1, 1),
                               std::make_tuple(-1, 1, 1)};
    std::string name = "cap";
    data_file_creation(name, cage);
}

//This function can return the terminal output. We will use it just for input.
std::string exec(std::string str_obj) {
    const char *cmd = &str_obj[0];
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
        throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length() - 1);
    return result;
}

coord definition_of_intersection_points(coord& initial_point, longDoubleTuple& beam) {
    double x, y, z;
    double x_init = std::get<0>(initial_point);
    double y_init = std::get<1>(initial_point);
    double z_init = std::get<2>(initial_point);
    double cos1 = std::get<2>(beam);
    double cos2 = std::get<1>(beam);
    double cos3 = std::get<0>(beam);
    double k = cos2 / cos1;
    double a = y_init - k * x_init;
    double h = cos3 / cos1;
    double b = z_init - h * x_init;
    int j = 0;
    do {
        double A = std::get<0>(planes[j]);
        double B = std::get<1>(planes[j]);
        double C = std::get<2>(planes[j]);
        double D = std::get<3>(planes[j]);
        x = -(B*a + C*b - D) / (A + B*k + C*h);
        y = k*x + a;
        z = h*x + b;
        j++;
    } while (!(x >= -1 && x <= 1 && y >= -1 && y <= 1 && z >= 0 && z <= 1) && !std::isnan(x));
    if (std::get<3>(beam) < std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2)) || std::isnan(x*y*z))
        return std::move(coordinates_of_the_interaction(beam));
    else
        return std::make_tuple(x, y, z);

}

std::tuple<double, double, double> statistical_weight (std::tuple<double, double, double> sigma) {
    double sum = std::get<0>(sigma) + std::get<1>(sigma) + std::get<2>(sigma);
    double p_Compton = std::get<0>(sigma) / sum;
    double p_ph = std::get<1>(sigma) / sum;
    double p_pp = std::get<2>(sigma) / sum;
    return std::make_tuple(p_Compton, p_ph, p_pp);
}

unsigned interaction_type (std::tuple<double, double, double>& p) {
    double gamma = eps * (rand() % (N+1));
    return (gamma <= std::get<0>(p)) ? 1 : (gamma <= std::get<1>(p)) ? 2 : 3;
}

coord vector_creation (coord& A, coord& B) {
    return std::make_tuple(std::get<0>(B) - std::get<0>(A),
                           std::get<1>(B) - std::get<1>(A),
                           std::get<2>(B) - std::get<2>(A));
}

template<typename T, size_t... Is>
auto abs_components_impl(T const& t, std::index_sequence<Is...>) {
    return std::sqrt((std::pow(std::get<Is>(t), 2) + ...));
}

template <class Tuple>
double abs_components(const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return abs_components_impl(t, std::make_index_sequence<size>{});
}

template<typename T, size_t... Is>
auto scalar_prod_components_impl(T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return ((std::get<Is>(t) * std::get<Is>(t1)) + ...);
}

template <class Tuple>
double scalar_prod_components(const Tuple& t, const Tuple& t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return scalar_prod_components_impl(t, t1,  std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}

//Function returns cos(a, b), where a, b -- vectors;
double cost(coord& A, coord& B, coord& C) {
    coord a = std::move(vector_creation(A, B));
    coord b = std::move(vector_creation(B, C));
    return scalar_prod_components(a, b) / (abs_components(a) * abs_components(b));
}

//The function returns energy steps for every particle.
std::array<std::vector<double>, N> interactions (std::vector<coord>& points) {
    std::tuple<double, double, double> p_air = statistical_weight(Sigma_air_2);
    std::tuple<double, double, double> p_Pb = statistical_weight(Sigma_Pb_2);
    std::vector<double> Energy;
    std::array<std::vector<double>, N> Energies;
    for(unsigned i = 0; i < points.size(); i++) {
        double alpha_min = E_min / E_0;
        double alpha = E_0 / E_e;
        unsigned type;
        bool flag = false;
        double sigma_sum;
        double x = std::get<0>(points[i]);
        double y = std::get<1>(points[i]);
        double z = std::get<2>(points[i]);
        coord A, C, B;
        longDoubleTuple direction = beam_direction(Sigma_air_sum);
        coord point_of_intersection = definition_of_intersection_points(points[i], direction);
        do {
            if (flag == 1) {
                if (std::abs(x) < 1 && std::abs(y) < 1 && z < 1) {
                    sigma_sum = Sigma_air_sum;
                    type = interaction_type(p_air);
                } else {
                    sigma_sum = Sigma_Pb_sum;
                    type = interaction_type(p_Pb);
                }
                direction = std::move(beam_direction(sigma_sum));
                C = std::make_tuple(std::get<2>(direction),
                                    std::get<1>(direction),
                                    std::get<0>(direction));
                double cos_ab = cost(A, B, C);
                A = B;
                B = std::move(definition_of_intersection_points(B, direction));
                x = std::get<0>(B);
                y = std::get<1>(B);
                z = std::get<2>(B);
                alpha /= 1 + (1 - cos_ab)*alpha;
            } else {
                A = points[i];
                B = point_of_intersection;
                flag = true;
            }
            Energy.emplace_back(alpha);
            interaction_points.at(i).emplace_back(B);
        } while (alpha > alpha_min || type == 2);
        Energies.at(i) = Energy;
    }
    return Energies;
}

void plot(std::array<std::vector<coord>, N>& points) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term wxt",
                                      "set multiplot",
                                      "set grid xtics ytics ztics",
                                      "set xrange [-5:5]",
                                      "set yrange [-5:5]",
                                      "set zrange [0:2]",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "splot '-' u 1:2:3 w lines"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    for(unsigned i = 0; i < N; i++) {
        for(unsigned j = 0; j < points[i].size(); j++) {
            double x = std::get<0>(points[i][j]);
            double y = std::get<1>(points[i][j]);
            double z = std::get<2>(points[i][j]);
            fprintf(gp, "%f\t%f\t%f'n", x, y, z);
        }
        fprintf(gp, "%c\n%s\n", 'e', "splot '-' u 1:2:3 w lines");
    }
    pclose(gp);
}