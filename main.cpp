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
#include <random>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <iterator>
#include <memory>
#include <stdexcept>


const int N = 2; //Number of points. //Do not use more than 0.5e4 on old computers!
constexpr double eps = 1.0 / N;
const double R = 1;
const double pi = 3.14159265359;

double E_min = 1.0e5;
const double E_0 = 2.0e6;
const double E_e = 0.51099895e6;

const double group_range = 1.0e5;

const int number_of_energy_groups = (E_0 - E_min) / group_range + 1; //Well, it's not safe I know.

std::vector<double> borders_of_groups (number_of_energy_groups);

std::array<std::vector<double>, N> Energies;

typedef std::tuple<double, double, double> coord;

typedef std::tuple<double, double, double, double> longDoubleTuple;

const std::vector<longDoubleTuple> planes = {std::make_tuple(0, 0, 1, -1), //There's a plane equation presented like Ax+By+Cz+D=0.
                                             std::make_tuple(0, 1, 0, -1),
                                             std::make_tuple(1, 0, 0, 1),
                                             std::make_tuple(0, 1, 0, 1),
                                             std::make_tuple(1, 0, 0, -1)};

std::vector <coord> polar ();

std::vector <coord> coordinate_transformation (std::vector<coord>& coords);

void data_file_creation (std::string DataType, std::vector<coord>& xx);

void data_file_creation (std::string DataType, std::vector<double>& xx, std::vector<coord>& yy);

void default_distribution_plot (std::string& name, std::string& data, std::string xlabel,
                                std::string ylabel, std::string title);

longDoubleTuple beam_direction (double sigma);

coord coordinates_of_the_interaction (longDoubleTuple& beams);

int energy_group(double& E);

std::array<std::vector<coord>, N> interactions (std::vector<coord>& points);

double cos_t (coord& A, coord& B, coord& C);

void cap();

void plot(std::array<std::vector<coord>, N>& points);

std::vector<longDoubleTuple> database_read (std::string name);

double linear_interpolation (double& x_0, double& y_0, double& x_1, double& y_1, double& x);

std::vector<std::tuple<double, double, double>> interpolated_database (std::vector<longDoubleTuple>& default_database);

std::string exec(std::string str);
//void flow_detection (std::vector<std::tuple<coord, double, int>>& inside_the_box, std::vector<coord>& detectors);

std::string PATH = exec("echo $PWD");

std::vector<std::tuple<double, double, double>> sigmas_air;
std::vector<std::tuple<double, double, double>> sigmas_Pb;

std::tuple<double, double, double> interpolation_for_single_particle (int& group, double& E,
                                                                      std::vector<std::tuple<double, double, double>>& sigmas);

std::vector<coord> detectors = {std::make_tuple(0, 0, 0.5),
                                std::make_tuple(0, 0, 1),
                                std::make_tuple(0, 0, 1.5),
                                std::make_tuple(0, -1, 0.5),
                                std::make_tuple(1, 0, 0.5)};

std::vector<std::vector<double>> eta;

void interpolation_plot (std::string matter, std::vector<double>& E,
                         std::vector<std::tuple<double, double, double>>& sigmas);


void path_def (std::string& path);

/*template<typename T, std::enable_if_t<std::is_floating_point<T>>* = nullptr>
T my_rand(T min, T max);*/

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

int main() {
    std::cout << "Databases collecting...\t";

    path_def(PATH);
    borders_of_groups.at(0) = (E_min);
    double E = E_min;
    std::generate(borders_of_groups.begin()+1, borders_of_groups.end(), [&] { return E += group_range; });
    std::vector<longDoubleTuple> air_data = std::move(database_read(PATH + "air_sigmas_database"));
    sigmas_air = std::move(interpolated_database(air_data));
    data_file_creation("air", borders_of_groups, sigmas_air);
    std::vector<longDoubleTuple> Pb_data = std::move(database_read(PATH + "Pb_sigmas_database"));
    sigmas_Pb = std::move(interpolated_database(air_data));
    data_file_creation("Pb", borders_of_groups, sigmas_Pb);

    std::cout << "Done!" << std::endl;


    std::cout << "Computing...\t";

    std::vector<coord> borning = std::move(polar());
    std::vector<coord> born_points = std::move(coordinate_transformation(borning));
    std::string name = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name, born_points);
    std::array<std::vector<coord>, N> interaction_points = std::move(interactions(born_points));

    std::cout << "Done!" << std::endl;


    std::cout << "Plotting...\t";

    default_distribution_plot(name, name, "x", "y", name);
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
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < N; i++) {
        double rho = R * std::sqrt(dis(gen));
        double phi = 2 * pi * dis(gen);
        double z = 0;
        coordinates.emplace_back(std::move(std::make_tuple(phi, rho, z)));
    }
    return coordinates;
}

//Function transform coordinates of points in polar system to Descart coordinate system.
std::vector<coord> coordinate_transformation (std::vector<coord>& coords) {
    std::vector<coord> xOy;
    for(int i = 0; i < coords.size(); i++) {
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
    fout.open(PATH + DataType);
    for(int i = 0; i < xx.size(); i++)
        fout << std::get<0>(xx[i]) << '\t' << std::get<1>(xx[i]) << '\t'
                                << std::get<2>(xx[i]) << '\t' << std::endl;
    fout.close();
}

//Function plots points of borning. It shows the initial distribution.
void default_distribution_plot (std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title) {
    //if you have problems with ".svg", you have to change ".svg" to ".pdf" in strings bellow.
    //std::cout << PATH + data << std::endl;
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error ("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term svg",
                                      "set out \'" + PATH + name + ".svg\'",
                                      "set xrange [-1:1]",
                                      "set yrange [-1:1]",
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set grid xtics ytics",
                                      "set title \'" + title + "\'",
                                      "plot \'" + PATH + data + "\' using 1:2 lw 1 lt rgb 'orange' ti \'Nodes\'",
                                      "set key box top right",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
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

//Function returns the beam direction and free run length for particle.
/* Note: always cos_gamma > 0, but it's n*/
longDoubleTuple beam_direction (double sigma) {
    //std::cout << sigma << std::endl;
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double L, mu, a, b, cos_psi, cos_gamma, d = 10;
    do {
        mu = 2 * dis(gen) - 1;
        do {
            a = 2 * dis(gen) - 1;
            b = 2 * dis(gen) - 1;
            d = std::pow(a, 2) + std::pow(b, 2);
            L = -log(dis(gen)) / sigma;
        } while (d > 1 || !std::isfinite(L));
        cos_psi = a / std::sqrt(d);
        cos_gamma = std::sqrt(1.0 - (std::pow(mu, 2) + std::pow(cos_psi, 2)));
    } while (std::pow(mu, 2) + std::pow(cos_psi, 2) > 1);
    //std::cout << cos_gamma << std::endl;
    return std::make_tuple(cos_gamma, mu, cos_psi, L);
}

//Function returns vector of beam direction.
coord coordinates_of_the_interaction (longDoubleTuple& beam) {
    double x = std::get<2>(beam) * std::get<3>(beam);
    double y = std::get<1>(beam) * std::get<3>(beam);
    double z = std::get<0>(beam) * std::get<3>(beam);
    //Разыгрывай по флагу знак.
    return std::make_tuple(x, y, z);
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


//The functions bellow using for solving systems of linear equations via LU-decomposition.
void LU_decomposition (std::vector<std::array<double, 3>>& default_matrix,
                       std::vector<std::array<double, 3>>& L,
                       std::vector<std::array<double, 3>>& U) {
    U = default_matrix;
    /*for (int i = 0; i < default_matrix.size(); i++) {
        for (int j = 0; j < default_matrix[i].size(); j++) {
            std::cout << U[i][j] << '\t';
        }
        std::cout << std::endl;
    }*/
    std::cout << std::endl;
    for (int i = 0; i < default_matrix.size(); i++) {
        for (int j = 0; j < default_matrix[i].size(); j++) {
            U.at(i).at(j) = 0;
            L.at(i).at(j) = 0;
        }
        L.at(i).at(i) = 1;
    }
    for (int i = 0; i < default_matrix.size(); i++) {
        for (int j = 0; j < default_matrix[i].size(); j++) {
            double sum = 0;
            for (int k = 0; k < i; k++)
                sum += L[i][k]*U[k][j];
            if (i <= j)
                U.at(i).at(j) = default_matrix[i][j] - sum;
            else
                L.at(i).at(j) = (default_matrix[i][j] - sum) / U[j][j];

        }
    }
}

std::vector<double> direct (std::vector<std::array<double, 3>>& L, std::vector<double>& b) {
    std::vector<double> y;
    for (int i = 0; i < L.size(); i++) {
        double sum = 0;
        for (int j = 0; j < i; j++)
            sum += L[i][j] * y[j];
        y.emplace_back(b[i] - sum);
    }
    return y;
}

std::vector<double> reverse (std::vector<std::array<double, 3>>& U, std::vector<double>& y) {
    std::vector<double> x = y;
    for (int i = y.size()-1; i >= 0; i--) {
        for (int j = i + 1; j < y.size(); j++)
            x.at(i) -= U[i][j] * x[j];
        x.at(i) /= U[i][i];
    }
    return x;
}

coord solve(std::vector<std::array<double, 3>> matrix, std::vector<double>& free_numbers_column) {
    std::vector<std::array<double, 3>> L(3), U(3);
    LU_decomposition(matrix, L, U);
    std::vector<double> straight_run_results = std::move(direct(L, free_numbers_column));
    std::vector<double> x = std::move(reverse(U, straight_run_results));
    return std::make_tuple(x[0], x[1], x[2]);
}

//First of all you need to analise is cos's and axes...
coord definition_of_intersection_points (coord& initial_point, longDoubleTuple& beam) {
    double x, x_init, y, y_init, z, z_init, cos_alpha, cos_beta, cos_gamma, A, B, C, D;
    cos_alpha = std::get<2>(beam);
    cos_beta = std::get<1>(beam);
    cos_gamma = std::get<0>(beam);
    coordinates_from_tuple(x_init, y_init, z_init, initial_point);
    coord intersection_coordinate;
    //std::cout << x_init << '\t' << y_init <<'\t' <<z_init << std::endl;
    int i = 0;
    do {
        A = std::get<0>(planes[i]);
        B = std::get<1>(planes[i]);
        C = std::get<2>(planes[i]);
        D = std::get<3>(planes[i]);
        std::cout << A << '\t' << B << '\t' << C << '\t' << D << std::endl;
        std::vector<std::array<double, 3>> matrix = {{cos_beta, -cos_alpha, 0},
                                                     {0, cos_gamma,  -cos_beta},
                                                     {A,        B,          C}};
        std::vector<double> right_part = {x_init*cos_beta - y_init*cos_alpha,
                                          y_init*cos_gamma - z_init*cos_beta,
                                          -D};
        intersection_coordinate = std::move(solve(matrix, right_part));
        coordinates_from_tuple(x, y, z, intersection_coordinate);
        i++;
        if  (x == x_init || y == y_init || z == z_init) continue;
    } while (i < planes.size() && (std::abs(x) <= 1 && std::abs(y) <= 1 && z > 0 && z <= 1));

    return intersection_coordinate;
}

//Just a function for creating vector with two points.
coord vector_creation (coord& A, coord& B) {
    return std::make_tuple(std::get<0>(B) - std::get<0>(A),
                           std::get<1>(B) - std::get<1>(A),
                           std::get<2>(B) - std::get<2>(A));
}

/* Function returns coordinate of interaction for particles.
 * If it interaction outside the sarcophagus functions returns the point of intersection with one of the planes,
 * otherwise it returns the interaction point in air. */
coord definition_of_interaction_points (coord& initial_point, longDoubleTuple& beam) {
    double x, x_init, y, y_init, z, z_init;
    coordinates_from_tuple(x_init, y_init, z_init, initial_point);
    coord trajectory;
    //std::cout << x_init << '\t' << y_init << '\t' << z_init << std::endl;
    // First of all we check if current frame of reference is inside the sarcophagus.
    if (std::abs(x_init) <= 1 && std::abs(y_init) <= 1 && z_init >= 0 && z_init <= 1) {
        coord intersection_point = std::move(definition_of_intersection_points(initial_point, beam));
        trajectory = std::move(vector_creation(initial_point, intersection_point));
        // Then, if the distance from one of the planes in this direction is bigger the free run length,
        // we will use last.
        if (std::get<3>(beam) < abs_components(trajectory))
            return std::move(coordinates_of_the_interaction(beam));
        else
            return intersection_point;
    } else { // Otherwise, when the current frame of reference inside the Pb, we have two ways too:
        coord interaction_point = std::move(coordinates_of_the_interaction(beam));
        trajectory = std::move(vector_creation(initial_point, interaction_point));
        coordinates_from_tuple(x, y, z, interaction_point);
        // if particle moves towards the inner walls of the sarcophagus
        // and the free run length lets her get to the sarcophagus we define the interaction point
        // like a point of sarcophagus and trajectory intersection.
        if (std::abs(x) <= 1 && std::abs(y) <= 1 && z <= 1 &&
        std::get<3>(beam) > abs_components(trajectory))
            return std::move(definition_of_intersection_points(interaction_point, beam));
        else // Otherwise we just define new frame of reference in Pb.
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
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double gamma = dis(gen);
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
                     double& E, std::string environment, coord& particle_coordinate) {
    int group = energy_group(E);
    std::vector<std::tuple<double, double, double>> sigma;
    if (environment == "air")
        sigma = sigmas_air;
    else
        sigma = sigmas_Pb;
    std::tuple<double, double, double> particle_sigma = std::move(interpolation_for_single_particle(group, E, sigma));
    sigma_sum = sum_components(particle_sigma);
    p = statistical_weight(particle_sigma, sigma_sum);
    type = interaction_type(p);
    double x, y, z;
    coordinates_from_tuple(x, y, z, particle_coordinate);
    /*for (int i = 0; i < detectors.size(); i++) {
        coordinates_from_tuple(x, y, z, detectors[i]);
        coord tau = vector_creation(particle_coordinate, detectors[i]);
        double distance = abs_components(tau);
        //Averadge, sum or max?
        /*if (z <= 1) //inside the box
            eta[i].emplace_back(p[0].first * std::exp(-distance)/std::pow(distance, 2) *
            std::get<0>(particle_sigma)/((std::get<0>(sigmas_air[group])+std::get<0>(sigmas_air[group+1]))/2.0));
        else
            eta[i].emplace_back(p[0].first * std::exp(-distance)/std::pow(distance, 2) *
            std::get<0>(particle_sigma)/((std::get<0>(sigmas_Pb[group])+std::get<0>(sigmas_air[group+1]))/2.0));
    }
    //Well, let's take average cross sections in group.
    /*    for (int i = 0; i < detectors.size(); i++) {
            coordinates_from_tuple(x_d, y_d, z_d, detectors[i]);
            if (z_d <= 1) {
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
    for (int i = 0; i < points.size(); i++) {
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
                    flow_detection(sigma_sum, p_air, type, E, "air", B);
                else
                    flow_detection(sigma_sum, p_Pb, type, E, "Pb", B);
                std::cout << x << '\t' << y << '\t' << z << std::endl;
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
            if (std::isnan(alpha)) break;
            coordinates_from_tuple(x, y, z, B);
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
    std::vector<std::string> stuff = {//"set term pdf",
                                      //"set output \'" + PATH + "test.pdf\'",
                                      "set term wxt",
                                      "set multiplot",
                                      "set grid xtics ytics ztics",
                                      "set xrange [-3:3]",
                                      "set yrange [-3:3]",
                                      "set zrange [0:3]",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "splot \'" + PATH + "Distribution of " + std::to_string(N) +" points\' u " +
                                                                                                  "1:2:3 lw 1 lt rgb 'red' ti \'Nodes\'",
                                      "splot \'" + PATH + "detectors\' u 1:2:3 lw 3 lt rgb 'black'",
                                      "splot \'" + PATH + "cap\' u 1:2:3 w lines lw 2 lt rgb 'black'",
                                      "splot \'" + PATH + "cap\' u 1:2:3 w boxes lw 2 lt rgb 'black'",
                                      "splot '-' u 1:2:3 w lines"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    double x, y, z;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < points[i].size(); j++) {
            coordinates_from_tuple(x, y, z, points[i][j]);
            fprintf(gp, "%f\t%f\t%f\n", x, y, z);
            //std::cout << x << '\t' << y << '\t' << z << std::endl;
        }
        fprintf(gp, "%c\n%s\n", 'e', "splot '-' u 1:2:3 w lines");
        std::cout << std::endl;
    }
    fprintf(gp, "%c\n", 'q');
    pclose(gp);
}


//Function return the number of energy group for particle.
int energy_group(double& E) {
    for (int i = borders_of_groups.size() - 1; i > 0; i--)
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
    int i, j, k;
    for (j = 0; j < default_database.size(); j++)
        if (std::get<0>(default_database[j]) == borders_of_groups[0])
            k = j;
    i = 0;
    j = k;
    double lb_E, lb_Compton, lb_ph, lb_pp, rb_E, rb_Compton, rb_ph, rb_pp, Compton, ph, pp;
    while (j < default_database.size()) {
        if (std::get<0>(default_database[j + 1]) - std::get<0>(default_database[j]) == group_range) {
            new_data.emplace_back(std::make_tuple(std::get<1>(default_database[j]),
                                                  std::get<2>(default_database[j]),
                                                  std::get<3>(default_database[j])));
            j++;
            i = (k == 0) ? j : j - k;
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
        } while (borders_of_groups[i] < rb_E && i < borders_of_groups.size());
        j++;
    }
    return new_data;
}

//The group determines by left border. The function returns interpolated interaction cross sections for single particle.
std::tuple<double, double, double> interpolation_for_single_particle (int& group, double& E,
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
                                      "set out \'" + PATH + "Photon Cross Sections for " + matter + ".svg\'",
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
    DataType += PATH;
    fout.open(DataType);
    for(int i = 0; i < xx.size(); i++)
        fout << xx[i] << '\t' << std::get<0>(yy[i]) << '\t' << std::get<1>(yy[i])
             << '\t' << std::get<2>(yy[i]) << '\t' << std::endl;
    fout.close();
}

//Just a function for returning the terminal output.
std::string exec(std::string str) {
    const char* cmd = str.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}

/* When you use cmake, files located in "cmake-build-debug",
 * when you use smth like "g++ -o %project_name main.cpp",
 * files located in root-directory of project.
 * So this function allows store files in root-directory of project
 * in both cases described. */
void path_def (std::string& path) {
    if (path.find("cmake-build-debug") != std::string::npos)
        path += "/../";
}



//std::array<std::array<double>,

/* We need to know number of virtual particles which detector registers.
 * So we have an array with coordinates and array with Energies for every particle.
 * let's collect their contribution in every detector. */

//TEST THE INTERSECTIONS ALGORITHM!