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
#include <algorithm>


const unsigned N = 0.5e3; //Number of points. //Do not use more than 0.5e4 on old computers!
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

//std::array<std::vector<coord>, N> interaction_points;

std::array<std::vector<double>, N> Energies;

std::vector <coord> polar ();

std::vector <coord> coordinate_transformation (std::vector<coord>& coords);

void data_file_creation (std::string DataType, std::vector<coord>& xx);

void default_distribution_plot (std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title);

longDoubleTuple beam_direction (double sigma);

coord coordinates_of_the_interaction (longDoubleTuple& beams);

std::vector<std::pair<double, std::string>> statistical_weight(std::tuple<double, double, double> sigma);

std::array<std::vector<coord>, N> interactions (std::vector<coord>& points);

double cos_t (coord& A, coord& B, coord& C);

void cap();

void plot(std::array<std::vector<coord>, N>& points);

std::vector<coord> detector_coordinates ();

int main() {
    srand(time(nullptr));
    std::vector<coord> borning = std::move(polar());
    std::vector<coord> born_points = std::move(coordinate_transformation(borning));
    std::string name1 = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name1, born_points);
    default_distribution_plot(name1, name1, "x", "y", name1);
    std::array<std::vector<coord>, N> interaction_points = std::move(interactions(born_points));
    std::cout << (Energies[1][0]) << std::endl;
    cap();
    std::vector<coord> detectors = detector_coordinates();
    data_file_creation("detectors", detectors);
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

std::vector<coord> coordinate_transformation (std::vector<coord>& coords) {
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

void data_file_creation (std::string DataType, std::vector<coord>& xx) {
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << std::get<0>(xx[i]) << '\t' << std::get<1>(xx[i]) << '\t' << std::get<2>(xx[i]) << '\t' << std::endl;
    fout.close();
}

void default_distribution_plot (std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title) {
    //if you have problems with ".svg", you have to change ".svg" to ".pdf" in strings bellow.
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error ("Error opening pipe to GNUplot.");
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

longDoubleTuple beam_direction (double sigma) {
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
    data_file_creation("cap", cage);
}

void coordinates_from_tuple (double& x, double& y, double& z, coord& point) {
    x = std::get<0>(point);
    y = std::get<1>(point);
    z = std::get<2>(point);
}

coord definition_of_intersection_points(coord& initial_point, longDoubleTuple& beam) {
    double x, x_init, y, y_init, z, z_init;
    coordinates_from_tuple(x_init, y_init, z_init, initial_point);
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

std::vector<std::pair<double, std::string>> statistical_weight (std::tuple<double, double, double> sigma) {
    std::vector<std::pair<double, std::string>> ans;
    double sum = std::get<0>(sigma) + std::get<1>(sigma) + std::get<2>(sigma);
    double p_Compton = std::get<0>(sigma) / sum;
    double p_ph = std::get<1>(sigma) / sum;
    double p_pp = std::get<2>(sigma) / sum;
    ans = {std::make_pair(p_Compton, "p_Compton"), std::make_pair(p_ph, "p_ph"), std::make_pair(p_pp, "p_pp")};
    return ans;
}

std::string interaction_type (std::vector<std::pair<double, std::string>>& p) {
    std::sort(p.begin(), p.end());
    double gamma = eps*(rand()%(N+1));
    return (gamma <= p[0].first) ? p[0].second : (gamma <= p[1].first) ? p[1].second : p[2].second;
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
double cos_t(coord& A, coord& B, coord& C) {
    coord a = std::move(vector_creation(A, B));
    coord b = std::move(vector_creation(B, C));
    return scalar_prod_components(a, b) / (abs_components(a) * abs_components(b));
}

//The function returns energy steps for every particle.
std::array<std::vector<coord>, N> interactions (std::vector<coord>& points) {
    std::vector<std::pair<double, std::string>> p_air = statistical_weight(Sigma_air_2);
    std::vector<std::pair<double, std::string>> p_Pb = statistical_weight(Sigma_Pb_2);
    std::vector<double> Energy;
    //std::array<std::vector<double>, N> Energies;
    std::array<std::vector<coord>, N> interaction_points;
    double x, y, z;
    for(unsigned i = 0; i < points.size(); i++) {
        //double alpha_min = E_min / E_0;
        interaction_points.at(i).emplace_back(points[i]);
        double alpha = E_0 / E_e;
        std::string type;
        bool flag = false;
        double sigma_sum;
        coordinates_from_tuple(x, y, z, points[i]);
        coord A, C, B;
        longDoubleTuple direction;
        coord point_of_intersection;
        do {
            if (flag == 1) {
                if (std::abs(x) <= 1 && std::abs(y) <= 1 && z <= 1) {
                    sigma_sum = Sigma_air_sum;
                    type = interaction_type(p_air);
                } else {
                    sigma_sum = Sigma_Pb_sum;
                    type = interaction_type(p_Pb);
                }
                direction = std::move(beam_direction(sigma_sum));
                C = definition_of_intersection_points(B, direction);
                double cos_ab = cos_t(A, B, C);
                A = B;
                B = C;
                x = std::get<0>(B);
                y = std::get<1>(B);
                z = std::get<2>(B);
                alpha /= 1 + (1 - cos_ab)*alpha;
            } else {
                A = points[i];
                direction = beam_direction(Sigma_air_sum);
                B = definition_of_intersection_points(A, direction);
                flag = true;
            }
            Energy.emplace_back(alpha);
            interaction_points.at(i).emplace_back(B);
        } while (type.empty() == 1 || type == "p_Compton"); // && alpha > alpha_min);
        Energies.at(i) = Energy;
    }
    return interaction_points;
}

void plot(std::array<std::vector<coord>, N>& points) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term pop",
                                      "set multiplot",
                                      "set grid xtics ytics ztics",
                                      "set xrange [-2:2]",
                                      "set yrange [-2:2]",
                                      "set zrange [0:5]",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "splot \'detectors\' u 1:2:3 lw 3 lt rgb 'black'",
                                      "splot \'cap\' u 1:2:3 w lines lw 2 lt rgb 'black'",
                                      "splot \'cap\' u 1:2:3 w boxes lw 2 lt rgb 'black'",
                                      "splot '-' u 1:2:3 w lines"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    double x, y, z;
    for(unsigned i = 0; i < N; i++) {
        for(unsigned j = 0; j < points[i].size(); j++) {
            coordinates_from_tuple(x, y, z, points[i][j]);
            fprintf(gp, "%f\t%f\t%f\n", x, y, z);
        }
        fprintf(gp, "%c\n%s\n", 'e', "splot '-' u 1:2:3 w lines");
    }
    pclose(gp);
}

std::vector<coord> detector_coordinates () {
    std::vector<coord> detectors = {std::make_tuple(0, 0, 0.5),
                                    std::make_tuple(0, 0, 1),
                                    std::make_tuple(0, 0, 1.5),
                                    std::make_tuple(0, -1, 0.5),
                                    std::make_tuple(1, 0, 0.5)};
    return detectors;
}

//Function returns coordinates and energies of particles inside the box.
std::vector<std::pair<coord, double>> inside(std::array<std::vector<coord>, N>& points) {
    std::vector<std::pair<coord, double>> inside_the_box;
    double x, y, z;
    for(unsigned i = 0; i < N; i++)
        for(unsigned j = 0; j < points[i].size(); j++) {
            coordinates_from_tuple(x, y, z, points[i][j]);
            if (std::abs(x) <= 1 && std::abs(y) <= 1 && z < 1) {
                inside_the_box.emplace_back(std::make_pair(std::make_tuple(x, y, z), Energies[i][j]));
            }
        }
    return inside_the_box;
}

//function returns an array of 20 energy groups in range from 1.0e5 eV to 2.0e6 eV.
//std::array<std::vector<double>, 20> energy_groups() {
//std::array<std::vector<double>, 20> E;
/*unsigned energy_group(double& E) {
    unsigned number_of_groups = 20;
    double E_init = 1.0e5;
    double E_final = 2.0e6;
    double group_range = (E_final - E_init) / number_of_groups;
    std::vector<double> borders (number_of_groups);
    std::generate(borders.begin(), borders.end(), [&] {return E_init += group_range;});
    /*for (unsigned i = 1; i < number_of_groups; i++)
        if (E > borders[i-1] && E < borders[i])
            return i;
        else continue;*/



//Try to make Energies global and Points local.


//std::vector<coord> inside (std::array<std::vector<coord>, N>& Energies) {

//}