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

const unsigned N = 1.0e3; //Number of points.
constexpr long double eps = 1.0 / N;
const long double R = 1;
const long double pi = 3.14159265359;
constexpr long double Sigma = 0.0349 + 0.00523 + 0.005;
bool flag = false;


typedef std::tuple<long double, long double, long double> coord;

typedef std::tuple<long double, long double, long double, long double> longDoubleTuple;


std::vector <coord> polar ();

std::vector <coord> coordinate_transformation(std::vector<coord>& coords);

void data_file_creation (std::string& DataType, std::vector<coord>& xx);

void default_distribution_plot(std::string& name, std::string& data, std::string xlabel, std::string ylabel, std::string title);

std::vector<longDoubleTuple> direction_beam();

std::vector<coord> coordinates_of_the_1st_interaction(std::vector<longDoubleTuple>& beams);

void cap();

std::string exec(const char* cmd);

void plot_of_the_1st_interaction(std::string& name);


int main() {
    srand(time(nullptr));
    std::vector<coord> borning = std::move(polar());
    std::vector<coord> points = std::move(coordinate_transformation(borning));
    std::string name1 = "Distribution of " + std::to_string(N) + " points";
    data_file_creation(name1, points);
    default_distribution_plot(name1, name1, "x", "y", name1);

    
    std::vector<longDoubleTuple>beams = std::move(direction_beam());
    std::vector<coord> first_interaction = std::move(coordinates_of_the_1st_interaction(beams));
    std::string name2 = "1st interaction";
    data_file_creation(name2, first_interaction);

    cap();
    char* char_arr;
    std::string str_obj("paste '" + name1 + "' '" + name2 + "' > test1");
    char_arr = &str_obj[0];
    exec(char_arr);
    plot_of_the_1st_interaction(name1);

    
    return 0;
}

std::vector <coord> polar () {
    std::vector <coord> coords;
    for (unsigned i = 0; i < N; i++) {
        long double rho = R * std::sqrt(eps * (rand() % (N + 1)));
        long double phi = 2 * pi * eps * (rand() % (N + 1));
        long double z = 0;
        coords.emplace_back(std::move(std::make_tuple(phi, rho, z)));
    }
    return coords;
}

std::vector<coord> coordinate_transformation(std::vector<coord>& coords) {
    std::vector<coord> xOy;
    for(unsigned i = 0; i < coords.size(); i++) {
        long double phi = std::get<0>(coords[i]);
        long double rho = std::get<1>(coords[i]);
        long double x = rho * cos(phi);
        long double y = rho * sin(phi);
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
    else {
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
}

//Function creates directions and lengths of flight of a particle before the first interaction.
std::vector<longDoubleTuple> direction_beam() {
    std::vector <longDoubleTuple> beams;
    for(unsigned i = 0; i < N; i++) {
        long double mu = 2 * eps * (rand() % (N + 1)) - 1; //cos(\phi);
        long double L, a, b, d = 10;
        do {
            a = 2 * eps * (rand() % (N + 1)) - 1;
            b = 2 * eps * (rand() % (N + 1)) - 1;
            d = std::pow(a, 2) + std::pow(b, 2);
            L = - log(eps * (rand() % (N + 1))) / Sigma;
        } while (d > 1 || std::isfinite(L) == 0);
        long double cosPsi = a / std::sqrt(d);
        long double sinPsi = b / std::sqrt(d);
        beams.emplace_back(mu, cosPsi, sinPsi, L);
    }
    return beams;
}

std::vector<coord> coordinates_of_the_1st_interaction(std::vector<longDoubleTuple>& beams) {
    std::vector<coord> coordinates;
    for(unsigned i = 0; i < beams.size(); i++) {
        long double x = std::get<0>(beams[i]) * std::get<3>(beams[i]);
        long double y = std::get<1>(beams[i]) * std::get<3>(beams[i]);
        long double z = std::get<2>(beams[i]) * std::get<3>(beams[i]);
        if (flag == 0)
            z = std::abs(z);
        coordinates.emplace_back(x, y, z);
    }
    flag = true;
    return coordinates;
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

void plot_of_the_1st_interaction(std::string& name) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp)
        throw std::runtime_error("Error opening pipe to GNUplot.");
    else {
        std::vector<std::string> stuff = {"set term svg",
                                          "set out \'" + name + ".svg\'",
                                          "set grid xtics ytics ztics",
                                          "set xrange [-1:1]",
                                          "set yrange [-1:1]",
                                          "set zrange [0:1]",
                                          "set ticslevel 0",
                                          "set border 4095",
                                          "splot \'" + name + "\',\
                                          \'cap\' w l,\
                                          \'test1\' w vectors",
                                          "set terminal pop",
                                          "set output",
                                          "replot"};
        for (const auto& it : stuff)
            fprintf(gp, "%s\n", it.c_str());
        pclose(gp);
    }
}

//This function can return the terminal output. We will use it just for input.
std::string exec(const char* cmd) {
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