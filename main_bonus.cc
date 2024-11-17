#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <optional>
#include <assert.h>
#include "platform_key.h"

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <string_view>

using std::cerr;
using std::cout;
using std::endl;
using std::nullopt;
using std::optional;
using std::string;
using std::string_view;
using std::unordered_map;
using std::vector;

using adjacency_matrix = std::vector<std::vector<size_t>>;

// Definition for simulate
void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
              const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
              const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
              size_t total_processes, const unordered_map<Platform_Key, vector<int>, Platform_Key_Hasher> &line_failures,
              const unordered_map<int, vector<int>> &train_failures);

enum LineColor {
    GREEN = 'g',
    YELLOW = 'y',
    BLUE = 'b',
    RED = 'r',
    BROWN = 'w',
    PURPLE = 'p',
    TURQUOISE = 't',
    PINK = 'k',
    LIME = 'l',
    GREY = 'e'
};

const LineColor colors[] = {GREEN, YELLOW, BLUE, RED, BROWN, PURPLE, TURQUOISE, PINK, LIME, GREY};

vector<string> extract_station_names(string &line) {
    constexpr char space_delimiter = ' ';
    vector<string> stations{};
    line += ' ';
    size_t pos;
    while ((pos = line.find(space_delimiter)) != string::npos) {
        stations.push_back(line.substr(0, pos));
        line.erase(0, pos + 1);
    }
    return stations;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << argv[0] << " <input_file>\n";
        std::exit(1);
    }

    int rank, tp;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tp);

    std::ifstream ifs(argv[1], std::ios_base::in);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open " << argv[1] << '\n';
        std::exit(2);
    }

    // Read S & V
    size_t S;
    size_t V;
    ifs >> S;
    ifs >> V;

    // Read station names.
    string station;
    std::vector<string> station_names{};
    station_names.reserve(S);
    for (size_t i = 0; i < S; ++i) {
        ifs >> station;
        station_names.emplace_back(station);
    }

    // Read P popularity
    size_t p;
    std::vector<size_t> popularities{};
    popularities.reserve(S);
    for (size_t i = 0; i < S; ++i) {
        ifs >> p;
        popularities.emplace_back(p);
    }

    // Form adjacency mat
    adjacency_matrix mat(S, std::vector<size_t>(S));
    for (size_t src{}; src < S; ++src) {
        for (size_t dst{}; dst < S; ++dst) {
            ifs >> mat[src][dst];
        }
    }

    ifs.ignore();

    string stations_buf;
    unordered_map<char, vector<string>> station_lines;

    for (size_t i = 0; i < V; ++i) {
        std::getline(ifs, stations_buf);
        vector<string> station_names = extract_station_names(stations_buf);
        station_lines[colors[i]] = std::move(station_names);
    }

    // N time ticks
    size_t N;
    ifs >> N;

    // number of trains per line
    unordered_map<char, size_t> num_trains;
    for (size_t i = 0; i < V; ++i) {
        size_t line_num_trains;
        ifs >> line_num_trains;
        num_trains[colors[i]] = line_num_trains;
    }

    size_t num_ticks_to_print;
    ifs >> num_ticks_to_print;

    // Read Bonus inputs
    size_t NLF, NTF;
    ifs >> NLF;
    ifs >> NTF;

    int tick;
    string station1, station2, type;
    unordered_map<Platform_Key, vector<int>, Platform_Key_Hasher> line_failures; // {platform_1: {tick_1, tick_2, ...}, ...}
    for (size_t i = 0; i < NLF; ++i) {
        ifs >> type;
        ifs >> station1;
        ifs >> station2;
        ifs >> tick;
        Platform_Key platform_key(station1, station2, station_names);
        auto iter = line_failures.find(platform_key);
        if (iter == line_failures.end()) {
            line_failures.emplace(std::piecewise_construct, std::forward_as_tuple(platform_key),
                                  std::forward_as_tuple(1, tick));
        } else {
            (iter->second).emplace_back(tick);
        }
    }

    int train_id;
    unordered_map<int, vector<int>> train_failures; // {tick: {train_1, train_2, ...}, ...}
    for (size_t i = 0; i < NTF; ++i) {
        ifs >> type;
        ifs >> train_id;
        ifs >> tick;
        auto iter = train_failures.find(tick);
        if (iter == train_failures.end()) {
            train_failures.emplace(std::piecewise_construct, std::forward_as_tuple(tick),
                                  std::forward_as_tuple(1, train_id));
        } else {
            (iter->second).emplace_back(train_id);
        }
    }

    // Start timing with MPI_Wtime
    double start_time = MPI_Wtime();

    // Call student implementation
    simulate(S, station_names, popularities, mat, station_lines, N, num_trains, num_ticks_to_print, rank, tp, line_failures, train_failures);

    // Barrier to make sure all processes are finished before timing
    MPI_Barrier(MPI_COMM_WORLD);

    // End timing with MPI_Wtime
    double end_time = MPI_Wtime();

    // We will use rank 0's timing in the end
    double total_time;
    if (rank == 0) {
        total_time = end_time - start_time;
        cerr << std::fixed << "mpi_time:" << total_time << "s" << endl;
    }

    MPI_Finalize();

    return 0;
}
