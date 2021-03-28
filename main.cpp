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
#include <algorithm>

const unsigned N = 1.0e3; //Number of points.
constexpr double eps = 1.0 / N;
const double R = 1;
const double pi = 3.14159265359;

bool flag = false;

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


std::vector <coord> polar ();

std::vector <coord> coordinate_transformation(std::vector<coord>& coords);

void data_file_creation (std::string& DataType, std::vector<coord>& xx);

void default_distribution_plot(std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title);

std::vector<longDoubleTuple> beam_direction (double sigma);

coord coordinates_of_the_interaction(longDoubleTuple& beams);

std::string exec(std::string str_obj);

void plot_of_the_1st_interaction(std::string& name);

std::vector<coord> definition_of_intersection_points (std::vector<coord>& initials, std::vector<longDoubleTuple>& beams);

std::tuple<double, double, double> statistical_weight(std::tuple<double, double, double> sigma);

unsigned interaction_type(std::tuple<double, double, double>& Sigma);

std::vector<unsigned> interaction ();

int main() {
    srand(time(nullptr));
    std::vector<coord> borning = std::move(polar());
    std::vector<coord> points = coordinate_transformation(borning);
    std::string name1 = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name1, points);
    default_distribution_plot(name1, name1, "x", "y", name1);

    std::vector<longDoubleTuple> beams = beam_direction(Sigma_air_sum);
    std::string name2 = "1st interaction";
    std::vector<coord> points2 = definition_of_intersection_points(points, beams);
    data_file_creation(name2, points2);
    exec("paste '" + name1 + "' '" + name2 + "' > test1");
    plot_of_the_1st_interaction(name1);
    flag = true; //The first interaction passed.



    /*if abs(x_i) == 1 ... sigma_...*/

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

//Function creates directions and lengths of flight of a particle before the first interaction.
std::vector<longDoubleTuple> beam_direction (double sigma) {
    std::vector <longDoubleTuple> beams;
    for(unsigned i = 0; i < N; i++) {
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
        beams.emplace_back(mu, cosPsi, sinPsi, L);
    }
    return beams;
}

coord coordinates_of_the_interaction (longDoubleTuple& beam) {
    double x = std::get<2>(beam) * std::get<3>(beam);
    double y = std::get<1>(beam) * std::get<3>(beam);
    double z = std::get<0>(beam) * std::get<3>(beam);
    if (flag == 0)
        z = std::abs(z);
    return std::make_tuple(x, y, z);
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
    const char* cmd = &str_obj[0];
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
        throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}

void plot_of_the_1st_interaction(std::string& name) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term svg",
                                      "set out \'" + name + ".svg\'",
                                      "set grid xtics ytics ztics",
                                      "set xrange [-1:1]",
                                      "set yrange [-1:1]",
                                      "set zrange [0:1]",
                                      //"set xlabel 'x'",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "splot \'" + name + "\',\
                                      \'cap\' w l,\
                                      \'test1\' w vectors nohead",
                                      "set terminal wxt",
                                      "set output",
                                      "replot"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

std::vector<coord> definition_of_intersection_points (std::vector<coord>& initials, std::vector<longDoubleTuple>& beams) {
    std::vector<coord> intersections;
    double x, y, z;
    for(unsigned i = 0; i < initials.size(); i++) {
        double x_init = std::get<0>(initials[i]);
        double y_init = std::get<1>(initials[i]);
        double z_init = std::get<2>(initials[i]);
        double cos1 = std::get<2>(beams[i]);
        double cos2 = std::get<1>(beams[i]);
        double cos3 = std::get<0>(beams[i]);
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
        if (std::get<3> (beams[i]) < std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2)))
            intersections.emplace_back(coordinates_of_the_interaction(beams[i]));
        else
            intersections.emplace_back(std::make_tuple(x, y, z));
    }
    return intersections;
}

std::tuple<double, double, double> statistical_weight (std::tuple<double, double, double>& sigma) {
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

std::vector<unsigned> interactions (std::vector<coord>& intersections) {
    /* Hint: There are some features for std::tuple that can make the function bellow shorter
     * in c++17 and c++20 like apply and for_each_in_tuple.*/
    std::tuple<double, double, double> p_air = statistical_weight(Sigma_air_2);
    std::tuple<double, double, double> p_Pb = statistical_weight(Sigma_Pb_2);
    for(unsigned i = 0; i < intersections.size(); i++) {
        std::tuple<double, double, double> Sigma;
        double x = std::get<0>(intersections[i]);
        double y = std::get<1>(intersections[i]);
        double z = std::get<2>(intersections[i]);
        do {
            if (std::abs(x) < 1 && std::abs(y) < 1 && z < 1)
                unsigned type = interaction_type(p_air);
            else
                unsigned type = interaction_type(p_Pb);

        } while (E > E_0);
    }
}

