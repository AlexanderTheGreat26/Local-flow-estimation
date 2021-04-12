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
#include <iterator>


const unsigned N = 1e2; //Number of points. //Do not use more than 0.5e4 on old computers!
constexpr double eps = 1.0 / N;
const double R = 1;
const double pi = 3.14159265359;

double E_min = 2.0e4;
const double E_0 = 2.0e6;
const double E_e = 0.51099895e6;


typedef std::tuple<double, double, double> coord;

typedef std::tuple<double, double, double, double> longDoubleTuple;

const std::vector<longDoubleTuple> planes = {std::make_tuple(0, 0, 1, -1), //There's a plane equation presented like Ax+By+Cz+D=0.
                                             std::make_tuple(0, 1, 0, -1),
                                             std::make_tuple(1, 0, 0, 1),
                                             std::make_tuple(0, 1, 0, 1),
                                             std::make_tuple(1, 0, 0, -1)};

const double group_range = 1.0e4;

unsigned number_of_energy_groups = (E_0 - E_min) / group_range + 1; //Well, it's not safe I know.

std::vector<double> borders_of_groups (number_of_energy_groups);

std::array<std::vector<double>, N> Energies;

std::vector <coord> polar ();

std::vector <coord> coordinate_transformation (std::vector<coord>& coords);

void data_file_creation (std::string DataType, std::vector<coord>& xx);

void data_file_creation (std::string DataType, std::vector<double>& xx, std::vector<coord>& yy);

void default_distribution_plot (std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title);

longDoubleTuple beam_direction (double sigma);

coord coordinates_of_the_interaction (longDoubleTuple& beams);

std::vector<std::pair<double, std::string>> statistical_weight(std::tuple<double, double, double> sigma);

unsigned energy_group(double& E);

std::array<std::vector<coord>, N> interactions (std::vector<coord>& points);

double cos_t (coord& A, coord& B, coord& C);

void cap();

void plot(std::array<std::vector<coord>, N>& points);

std::vector<std::tuple<coord, double, unsigned >> inside (std::array<std::vector<coord>, N>& points);

std::vector<longDoubleTuple> database_read (std::string name);

double linear_interpolation (double& x_0, double& y_0, double& x_1, double& y_1, double& x);

std::vector<std::tuple<double, double, double>> interpolated_database (std::vector<longDoubleTuple>& default_database);

//void flow_detection (std::vector<std::tuple<coord, double, unsigned>>& inside_the_box, std::vector<coord>& detectors);

std::vector<std::tuple<double, double, double>> sigmas_air; // = std::move(database_read("air_sigmas_database"));
std::vector<std::tuple<double, double, double>> sigmas_Pb; // = std::move(database_read("Pb_sigmas_database"));

std::tuple<double, double, double> interpolation_for_single_particle (unsigned& group, double& E,
                                                                      std::vector<std::tuple<double, double, double>>& sigmas);

std::vector<coord> detectors = {std::make_tuple(0, 0, 0.5),
                                std::make_tuple(0, 0, 1),
                                std::make_tuple(0, 0, 1.5),
                                std::make_tuple(0, -1, 0.5),
                                std::make_tuple(1, 0, 0.5)};

void interpolation_plot (std::string matter, std::vector<double>& E, std::vector<std::tuple<double, double, double>>& sigmas);

int main() {
    std::cout << "Databases collecting...\t";

    borders_of_groups.at(0) = (E_min);
    double E = E_min;
    std::generate(borders_of_groups.begin()+1, borders_of_groups.end(), [&] { return E += group_range; });
    std::vector<longDoubleTuple> air_data = std::move(database_read("air_sigmas_database"));
    sigmas_air = std::move(interpolated_database(air_data));
    data_file_creation("air", borders_of_groups, sigmas_air);
    std::vector<longDoubleTuple> Pb_data = std::move(database_read("Pb_sigmas_database"));
    sigmas_Pb = std::move(interpolated_database(air_data));
    data_file_creation("Pb", borders_of_groups, sigmas_Pb);

    std::cout << "Done!" << std::endl;


    std::cout << "Computing...\t";

    srand(time(nullptr));
    std::vector<coord> borning = std::move(polar());
    std::vector<coord> born_points = std::move(coordinate_transformation(borning));
    std::string name1 = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name1, born_points);
    default_distribution_plot(name1, name1, "x", "y", name1);
    std::array<std::vector<coord>, N> interaction_points = std::move(interactions(born_points));

    std::cout << "Done!" << std::endl;


    std::cout << "Plotting...\t";

    interpolation_plot("air", borders_of_groups, sigmas_air);
    interpolation_plot("Pb", borders_of_groups, sigmas_air);
    cap();
    data_file_creation("detectors", detectors);
    plot(interaction_points);

    std::cout << "Done!" << std::endl;


    return 0;
}


//Function returns random points generated in polar coordinate system.
std::vector<coord> polar () {
    std::vector <coord> coordinates;
    for (unsigned i = 0; i < N; i++) {
        double rho = R * std::sqrt(eps * (rand() % (N + 1)));
        double phi = 2 * pi * eps * (rand() % (N + 1));
        double z = 0;
        coordinates.emplace_back(std::move(std::make_tuple(phi, rho, z)));
    }
    return coordinates;
}

//Function transform coordinates of points in polar system to Descart coordinate system.
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

//Function creates a data-file with coordinates. It's useful for plotting.
void data_file_creation (std::string DataType, std::vector<coord>& xx) {
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << std::get<0>(xx[i]) << '\t' << std::get<1>(xx[i]) << '\t' << std::get<2>(xx[i]) << '\t' << std::endl;
    fout.close();
}

//Function plots points of borning. It shows the initial distribution.
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
                                      "replot", "q"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

//Function returns the beam direction and free run length for particle.
longDoubleTuple beam_direction (double sigma) {
    double mu = 2 * eps * (rand() % (N + 1)) - 1; //cos(\phi);
    double L, a, b, d = 10;
    do {
        a = 2 * eps * (rand() % (N + 1)) - 1;
        b = 2 * eps * (rand() % (N + 1)) - 1;
        d = std::pow(a, 2) + std::pow(b, 2);
        L = - log(eps * (rand() % (N+1))) / sigma;
    } while (d > 1 || std::isfinite(L) == 0);
    double cosPsi = a / std::sqrt(d);
    double sinPsi = b / std::sqrt(d);
    return std::make_tuple(mu, cosPsi, sinPsi, L);
}

//Function returns vector of beam direction.
coord coordinates_of_the_interaction (longDoubleTuple& beam) {
    double x = std::get<2>(beam) * std::get<3>(beam);
    double y = std::get<1>(beam) * std::get<3>(beam);
    double z = std::get<0>(beam) * std::get<3>(beam);
    return std::make_tuple(x, y, std::abs(z));
}

//Just a function for sarcophagus creation.
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

//The templates below will be used for defining the angle between the vectors.
template<typename T, size_t... Is>
auto abs_components_impl(T const& t, std::index_sequence<Is...>) {
    return std::sqrt((std::pow(std::get<Is>(t), 2) + ...));
}

template <class Tuple>
double abs_components(const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return abs_components_impl(t, std::make_index_sequence<size>{});
}

coord definition_of_intersection_points (coord& initial_point, longDoubleTuple& beam) {
    double x, x_init, y, y_init, z, z_init, cos1, cos2, cos3, k, a, h, b, A, B, C, D;
    coordinates_from_tuple(x_init, y_init, z_init, initial_point);
    cos1 = std::get<2>(beam);
    cos2 = std::get<1>(beam);
    cos3 = std::get<0>(beam);
    k = cos2 / cos1;
    a = y_init - k * x_init;
    h = cos3 / cos1;
    b = z_init - h * x_init;
    unsigned j = 0;
    do {
        A = std::get<0>(planes[j]);
        B = std::get<1>(planes[j]);
        C = std::get<2>(planes[j]);
        D = std::get<3>(planes[j]);
        x = -(B*a + C*b - D) / (A + B*k + C*h);
        y = k*x + a;
        z = h*x + b;
        j++;
    } while (!(x >= -1 && x <= 1 && y >= -1 && y <= 1 && z >= 0 && z <= 1) && !std::isnan(x));
    return std::make_tuple(x, y, z);
}

/* Function returns coordinate of interaction for particles.
 * If it interaction outside the sarcophagus functions returns the point of intersection with one of the planes,
 * otherwise it returns the interaction point in air. */
coord definition_of_interaction_points (coord& initial_point, longDoubleTuple& beam) {
    double x, x_init, y, y_init, z, z_init;
    coordinates_from_tuple(x_init, y_init, z_init, initial_point);
    if (x_init > -1 && x_init < 1 && y_init > -1 && y_init < 1 && z_init >= 0 && z_init < 1) {
        coord intersection_point = std::move(definition_of_intersection_points(initial_point, beam));;
        if (std::get<3>(beam) < abs_components(intersection_point))
            return std::move(coordinates_of_the_interaction(beam));
        else
            return intersection_point;
    } else {
        coord interaction_point = std::move(coordinates_of_the_interaction(beam));
        coordinates_from_tuple(x, y, z, interaction_point);
        if (x >= -1 && x <= 1 && y >= -1 && y <= 1 && z <= 1) {
            return std::move(definition_of_intersection_points(interaction_point, beam));
        } else
            return interaction_point;
    }
}

//Function returns the probability for types of interaction for environment (which defines with argument).
std::vector<std::pair<double, std::string>> statistical_weight (std::tuple<double, double, double>& sigma,
                                                                double& sum) {
    std::vector<std::pair<double, std::string>> ans;
    double p_Compton = std::get<0>(sigma) / sum;
    double p_ph = std::get<1>(sigma) / sum;
    double p_pp = std::get<2>(sigma) / sum;
    ans = {std::make_pair(p_Compton, "Compton"), std::make_pair(p_ph, "ph"), std::make_pair(p_pp, "pp")};
    return ans;
}

//Function return random interaction type.
std::string interaction_type (std::vector<std::pair<double, std::string>>& p) {
    std::sort(p.begin(), p.end());
    double gamma = eps*(rand()%(N+1));
    return (gamma <= p[0].first) ? p[0].second : (gamma <= p[1].first) ? p[1].second : p[2].second;
}

template<typename T, size_t... Is>
auto sum_components_impl(T const& t, std::index_sequence<Is...>) {
    return (std::get<Is>(t) + ...);
}

template <class Tuple>
double sum_components(const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return sum_components_impl(t, std::make_index_sequence<size>{});
}

//Just a function for creating vector with two points.
coord vector_creation (coord& A, coord& B) {
    return std::make_tuple(std::get<0>(B) - std::get<0>(A),
                           std::get<1>(B) - std::get<1>(A),
                           std::get<2>(B) - std::get<2>(A));
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

//Function returns linear interpolation of interaction cross section for particles among neighboring borders.
double linear_interpolation (double& x_0, double& y_0, double& x_1, double& y_1, double& x) {
    return y_0 + (y_1 - y_0)/(x_1 - x_0) * (x - x_0);
}

void flow_detection (double& sigma_sum, std::vector<std::pair<double, std::string>>& p, std::string& type,
                     double& E, std::string environment) {
    unsigned group = energy_group(E);
    std::vector<std::tuple<double, double, double>> sigma;
    if (environment == "air")
        sigma = sigmas_air;
    else
        sigma = sigmas_Pb;
    std::tuple<double, double, double> particle_sigma = std::move(interpolation_for_single_particle(group, E, sigma));
    sigma_sum = sum_components(particle_sigma);
    p = statistical_weight(particle_sigma, sigma_sum);
    type = interaction_type(p);
    /*    for (unsigned i = 0; i < detectors.size(); i++) {
            coordinates_from_tuple(x_d, y_d, z_d, detectors[i]);
            if (z_d <= 1) {
                coord tau = vector_creation(particle_coordinate, detectors[i]);
                doubledistance = abs_components(tau);
                eta[group] += p[0].first * std::exp(-distance) / std::pow(distance, 2) * (particle_sigma)
                / std::get<0>(sigmas_air[group]);
            }
        }
    } else {

    }*/
}

//The function returns energy steps for every particle.
std::array<std::vector<coord>, N> interactions (std::vector<coord>& points) {
    std::vector<std::pair<double, std::string>> p_air, p_Pb;
    std::vector<std::tuple<double, double, double>> sigmas;
    std::vector<double> Energy;
    std::array<std::vector<coord>, N> interaction_points;
    double sigma_2_air_sum = sum_components(sigmas_air[sigmas_air.size()-1]);
    double x, y, z, alpha, E, sigma_sum, cos_ab;
    for (unsigned i = 0; i < points.size(); i++) {
        interaction_points.at(i).emplace_back(points[i]);
        alpha = E_0 / E_e;
        std::string type;
        bool flag = false;
        coord A, C, B;
        longDoubleTuple direction;
        do {
            if (flag == 1) {
                E = alpha * E_e;
                Energy.emplace_back(E);
                interaction_points.at(i).emplace_back(B);
                if (std::abs(x) <= 1 && std::abs(y) <= 1 && z <= 1)
                    flow_detection(sigma_sum, p_air, type, E, "air");
                else
                    flow_detection(sigma_sum, p_Pb, type, E, "Pb");
                direction = std::move(beam_direction(sigma_sum));
                C = definition_of_interaction_points(B, direction);
                cos_ab = cos_t(A, B, C);
                A = B;
                B = C;
                alpha /= 1 + (1 - cos_ab)*alpha;
            } else {
                A = points[i];
                direction = std::move(beam_direction(sigma_2_air_sum));
                B = definition_of_interaction_points(A, direction);
                flag = true;
            }
            coordinates_from_tuple(x, y, z, B);
            if (std::isnan(alpha)) break;
        } while (type.empty() == 1 || type == "Compton" && E > E_min);
        Energies.at(i) = Energy;
    }
    return interaction_points;
}

//The function plots the trajectories of particles. So it not fast, so you can comment it in main().
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
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < points[i].size(); j++) {
            coordinates_from_tuple(x, y, z, points[i][j]);
            fprintf(gp, "%f\t%f\t%f\n", x, y, z);
        }
        fprintf(gp, "%c\n%s\n", 'e', "splot '-' u 1:2:3 w lines");
    }
    fprintf(gp, "%c\n", 'q');
    pclose(gp);
}


//Function return the number of energy group for particle.
unsigned energy_group(double& E) {
    for (unsigned i = borders_of_groups.size() - 1; i > 0; i--)
        if (E >= borders_of_groups[i - 1] && E <= borders_of_groups[i])
            return i - 1;
}

namespace std {
    istream& operator >> (istream& in, longDoubleTuple& data) {
        double first, second, third, fourth;
        in >> first >> second >> third >> fourth;
        data = {first, second, third, fourth};
        return in;
    }

    ostream& operator << (ostream& out, const longDoubleTuple& data) {
        auto [first, second, third, fourth] = data;
        out << first << ' ' << second << ' ' << third << ' ' << fourth << ' ';
        return out;
    }
}

/* Function returns data in std::vector of std::tuple from text file include the interaction cross section for particle
 * with determined energy in environment. File looks like matrix 3xN. */
std::vector<longDoubleTuple> database_read (std::string name) {
    std::ifstream inFile(name);
    std::vector<longDoubleTuple> tuples_vector;
    copy(std::istream_iterator<longDoubleTuple> {inFile},
         std::istream_iterator<longDoubleTuple> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<longDoubleTuple>(std::cout, "\n"));
    return tuples_vector;
}

//Function returns database received via linear interpolation. It will be useful for presentation of results.
std::vector<std::tuple<double, double, double>> interpolated_database (std::vector<longDoubleTuple>& default_database) {
    std::vector<std::tuple<double, double, double>> new_data;
    unsigned i, j;
    i = j = 0;
    double lb_E, lb_Compton, lb_ph, lb_pp, rb_E, rb_Compton, rb_ph, rb_pp, Compton, ph, pp;
    while (i < borders_of_groups.size()) {
        if (std::get<0>(default_database[j + 1]) - std::get<0>(default_database[j]) == group_range) {
            new_data.emplace_back(std::make_tuple(std::get<1>(default_database[j]),
                                                  std::get<2>(default_database[j]),
                                                  std::get<3>(default_database[j])));
            j++;
            i = j;
            continue;
        }
        lb_E = std::get<0>(default_database[j]);
        lb_Compton = std::get<1>(default_database[j]);
        lb_ph = std::get<2>(default_database[j]);
        lb_pp = std::get<3>(default_database[j]);
        rb_E = std::get<0>(default_database[j+1]);
        rb_Compton = std::get<1>(default_database[j+1]);
        rb_ph = std::get<2>(default_database[j+1]);
        rb_pp = std::get<3>(default_database[j+1]);
        do {
            Compton = linear_interpolation(lb_E, lb_Compton, rb_E, rb_Compton, borders_of_groups[i]);
            ph = linear_interpolation(lb_E, lb_ph, rb_E, rb_ph, borders_of_groups[i]);
            pp = linear_interpolation(lb_E, lb_pp, rb_E, rb_pp, borders_of_groups[i]);
            new_data.emplace_back(std::make_tuple(Compton, ph, pp));
            i++;
        } while (borders_of_groups[i] < rb_E);
        j++;
    }
    return new_data;
}

//The group determines by left border. The function returns interpolated interaction cross sections for single particle.
std::tuple<double, double, double> interpolation_for_single_particle (unsigned& group, double& E,
                                                                      std::vector<std::tuple<double, double, double>>& sigmas) {
    double lb_E = borders_of_groups[group];
    double lb_Compton = std::get<0>(sigmas[group]);
    double lb_ph = std::get<1>(sigmas[group]);
    double lb_pp = std::get<2>(sigmas[group]);
    double rb_E = borders_of_groups[group+1];
    double rb_Compton = std::get<0>(sigmas[group+1]);
    double rb_ph = std::get<1>(sigmas[group+1]);
    double rb_pp = std::get<2>(sigmas[group+1]);
    double Compton = Compton = linear_interpolation(lb_E, lb_Compton, rb_E, rb_Compton, E);
    double ph = linear_interpolation(lb_E, lb_ph, rb_E, rb_ph, E);
    double pp = linear_interpolation(lb_E, lb_pp, rb_E, rb_pp, E);
    return std::make_tuple(Compton, ph, pp);
}

void interpolation_plot(std::string matter, std::vector<double>& E, std::vector<std::tuple<double, double, double>>& sigmas) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term svg",
                                      "set out \'Photon Cross Sections for " + matter + ".svg\'",
                                      "set grid xtics ytics",
                                      "set title \'Photon Cross Sections for " + matter + "\'",
                                      "plot \'" + matter + "\' u 1:2 w lines ti \'Compton Scattering\',\'"
                                      + matter + "\' u 1:2 lw 1 lt rgb 'black' ti \'Compton Scattering nodes\',\'"
                                      + matter + "\' u 1:3 w lines ti \'Photoelectric\',\'"
                                      + matter + "\' u 1:3 lw 1 lt rgb 'black' ti \'Photoelectric nodes\',\'"
                                      + matter + "\' u 1:4 w lines ti \'Pair Production\',\'"
                                      + matter + "\' u 1:4 lw 1 lt rgb 'black' ti \'Pair Production nodes\'",
                                      "set xlabel \'Energy, MeV\'",
                                      "set ylabel \'Cross sections, cm^2/g\'",
                                      "set terminal wxt",
                                      "set output",
                                      "replot", "q"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}

void data_file_creation (std::string DataType, std::vector<double>& xx, std::vector<coord>& yy) {
    //For reading created files via Matlab use command: M = dlmread('/PATH/file'); xi = M(:,i);
    std::ofstream fout;
    fout.open(DataType);
    for(unsigned i = 0; i < xx.size(); i++)
        fout << xx[i] << '\t' << std::get<0>(yy[i]) << '\t' << std::get<1>(yy[i])
             << '\t' << std::get<2>(yy[i]) << '\t' << std::endl;
    fout.close();
}