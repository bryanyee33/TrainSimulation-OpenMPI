#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <deque>
#include <iostream>

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include "platform_load_time_gen.hpp"

using std::string;
using std::unordered_map;
using std::vector;
using std::deque;
using adjacency_matrix = std::vector<std::vector<size_t>>;

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

struct Platform {
    int platform;
    int rank;
    Platform(int platform, int rank) : platform(platform), rank(rank) {}
    Platform() : platform(-1), rank(-1) {}
};

inline int pair_idx(vector<int> &keys, vector<int> &vals, int key, int val) {
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] == key && vals[i] == val) {
            return i;
        }
    }
    return -1;
}

inline void assign_platforms(int platforms_per_process, int platform_idx, vector<int> &platform_from, vector<int> &platform_to,
        vector<int> &platform_dist, vector<vector<int>> &platform_colors, vector<vector<int>> &platform_spawns, vector<vector<bool>> &platform_spawns_is_last,
        vector<int> &all_platform_from, vector<int> &all_platform_to, const adjacency_matrix &mat,
        vector<vector<int>> &station_lines_vec, vector<int> &num_trains_vec) {

    int temp;
    for (size_t color = 0; color < station_lines_vec.size(); ++color) {
        vector<int> &station_line = station_lines_vec[color];

        // Forwards direction
        for (size_t j = 0; j < station_line.size() - 1; ++j) {
            if ((temp = pair_idx(platform_from, platform_to, station_line[j], station_line[j + 1])) != -1) { // Platform already added & incharge
                platform_colors[temp].emplace_back(color);

                if (j == 0) { // Starting terminal platform
                    platform_spawns[temp][color] = (num_trains_vec[color] + 1) / 2; // Round up
                }
                continue;
            } else if (pair_idx(all_platform_from, all_platform_to, station_line[j], station_line[j + 1]) != -1) { // Platform already added
                continue;
            } else if (-platforms_per_process < platform_idx && platform_idx <= 0) { // Found unadded platform for this process
                platform_from.emplace_back(station_line[j]);
                platform_to.emplace_back(station_line[j + 1]);
                platform_dist.emplace_back(mat[station_line[j]][station_line[j + 1]]);
                platform_colors[-platform_idx].emplace_back(color);

                if (j == 0) { // Starting terminal platform
                    platform_spawns[-platform_idx][color] = (num_trains_vec[color] + 1) / 2; // Round up
                }
            }
            all_platform_from.emplace_back(station_line[j]);
            all_platform_to.emplace_back(station_line[j + 1]);
            --platform_idx;
        }

        // Backwards direction
        for (size_t j = station_line.size() - 1; j > 0; --j) {
            if ((temp = pair_idx(platform_from, platform_to, station_line[j], station_line[j - 1])) != -1) { // Platform already added & incharge
                platform_colors[temp].emplace_back(color);

                if (j == station_line.size() - 1) { // Ending terminal platform
                    platform_spawns[temp][color] = num_trains_vec[color] / 2;
                    platform_spawns_is_last[temp][color] = true;
                }
                continue;
            } else if (pair_idx(all_platform_from, all_platform_to, station_line[j], station_line[j - 1]) != -1) { // Platform already added
                continue;
            } else if (-platforms_per_process < platform_idx && platform_idx <= 0) { // Found unadded platform for this process
                platform_from.emplace_back(station_line[j]);
                platform_to.emplace_back(station_line[j - 1]);
                platform_dist.emplace_back(mat[station_line[j]][station_line[j - 1]]);
                platform_colors[-platform_idx].emplace_back(color);

                if (j == station_line.size() - 1) { // Ending terminal platform
                    platform_spawns[-platform_idx][color] = num_trains_vec[color] / 2;
                    platform_spawns_is_last[-platform_idx][color] = true;
                }
            }
            all_platform_from.emplace_back(station_line[j]);
            all_platform_to.emplace_back(station_line[j - 1]);
            --platform_idx;
        }
    }
}

inline void find_connected_platforms(vector<int> &platform_from, vector<int> &platform_to, vector<unordered_map<int, Platform>> &outgoing,
        vector<unordered_map<int, Platform>> &incoming, int platforms_per_process, vector<vector<int>> &platform_colors,
        vector<int> &all_platform_from, vector<int> &all_platform_to, vector<vector<int>> &station_lines_vec) {

    for (size_t platform = 0; platform < platform_to.size(); ++platform) {
        for (int color : platform_colors[platform]) {
            vector<int> &station_line = station_lines_vec[color];

            size_t from_idx = std::find(station_line.begin(), station_line.end(), platform_from[platform]) - station_line.begin();
            size_t to_idx = std::find(station_line.begin(), station_line.end(), platform_to[platform]) - station_line.begin();
            int before_from, after_to, result;

            if (from_idx < to_idx) {
                if (from_idx == 0) {
                    before_from = station_line[1];
                } else {
                    before_from = station_line[from_idx - 1];
                }
                if (to_idx == station_line.size() - 1) {
                    after_to = station_line[station_line.size() - 2];
                } else {
                    after_to = station_line[to_idx + 1];
                }
                
            } else {
                if (from_idx == station_line.size() - 1) {
                    before_from = station_line[station_line.size() - 2];
                } else {
                    before_from = station_line[from_idx + 1];
                }
                if (to_idx == 0) {
                    after_to = station_line[1];
                } else {
                    after_to = station_line[to_idx - 1];
                }
            }

            result = pair_idx(all_platform_from, all_platform_to, platform_to[platform], after_to);
            outgoing[platform].emplace(std::piecewise_construct, std::forward_as_tuple(color),
                                       std::forward_as_tuple(result, result / platforms_per_process));

            result = pair_idx(all_platform_from, all_platform_to, before_from, platform_from[platform]);
            incoming[platform].emplace(std::piecewise_construct, std::forward_as_tuple(color),
                                       std::forward_as_tuple(result, result / platforms_per_process));
        }
    }
}

void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
              const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
              const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
              size_t total_processes) {

    vector<int> num_trains_vec;
    for (size_t color = 0; color < num_trains.size(); ++color) {
        num_trains_vec.emplace_back(num_trains.at(colors[color]));
    }

    int num_platforms = 0;
    for (size_t i = 0; i < num_stations; ++i) {
        for (size_t j = 0; j < num_stations; ++j) {
            if (mat[i][j] != 0) {
                ++num_platforms;
            }
        }
    }

    // Station lines in terms of indices
    vector<vector<int>> station_lines_vec(station_lines.size());
    for (size_t color = 0; color < station_lines.size(); ++color) {
        for (const string &station_name: station_lines.at(colors[color])) {
            station_lines_vec[color].emplace_back(std::find(station_names.begin(), station_names.end(), station_name) - station_names.begin());
        }
    }

    int num_colors = station_lines.size();
    int platforms_per_process = (num_platforms + total_processes - 1) / total_processes; // Round up
    int platform_startidx = platforms_per_process * mpi_rank;
    vector<int> platform_from, platform_to, platform_dist; // Only this process' platforms
    vector<vector<int>> platform_colors(platforms_per_process); // Colors for this process' platforms
    vector<vector<int>> platform_spawns(platforms_per_process, vector<int>(num_colors, 0));
    vector<vector<bool>> platform_spawns_is_last(platforms_per_process, vector<bool>(num_colors, false));
    vector<int> all_platform_from, all_platform_to; // Sorted based on color, all process' platforms
    all_platform_from.reserve(num_platforms);
    all_platform_to.reserve(num_platforms);

    assign_platforms(platforms_per_process, platform_startidx, platform_from, platform_to, platform_dist, platform_colors, platform_spawns,
            platform_spawns_is_last, all_platform_from, all_platform_to, mat, station_lines_vec, num_trains_vec);
    
    
    int num_platforms_rank = platform_to.size(); // No. of platforms this rank is in charge of
    // [{pform_1_color_idx_1: pform_1_rank_1, pform_1_color_idx_2: pform_1_rank_2, ...}, {pform_2_...}]
    vector<unordered_map<int, Platform>> outgoing(num_platforms_rank);
    vector<unordered_map<int, Platform>> incoming(num_platforms_rank);

    find_connected_platforms(platform_from, platform_to, outgoing, incoming, platforms_per_process, platform_colors, all_platform_from, all_platform_to, station_lines_vec);

    vector<PlatformLoadTimeGen> pltg;
    for (int platform = 0; platform < num_platforms_rank; ++platform) {
        pltg.emplace_back(popularities[platform_from[platform]]);
    }

    int NEG_ONE = -1, out_color;
    bool spawns_left = true;
    vector<int> travelling(num_platforms_rank);
    vector<int> travel_time_left(num_platforms_rank, -1); // -1 for unused
    vector<int> unload_time_left(num_platforms_rank, -1);
    vector<int> train_color; // Sorted by train id
    train_color.reserve(ticks << 1); // Max of 2 * ticks trains
    vector<deque<int>> trains(num_platforms_rank); // First element is in platform. The rest are in holding area.

    vector<int> sendbuff(num_platforms_rank); // Each platform sends at most 1 train
    vector<int> recvbuff(num_colors * num_platforms_rank, -1); // Each platform can receive trains up to the number of colors available
    vector<MPI_Request> recvstats(num_colors * num_platforms_rank, MPI_REQUEST_NULL);
    MPI_Request req; // For unneeded output

    for (int tick = 0; tick < (int)ticks; ++tick) {
        // Update and send trains
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            if (--travel_time_left[platform] == 0) {
                sendbuff[platform] = travelling[platform];
                out_color = train_color[travelling[platform]];
            } else {
                out_color = -1;
            }

            for (int color : platform_colors[platform]) {
                if (color == out_color) {
                    MPI_Isend(&sendbuff[platform], 1, MPI_INT, outgoing[platform][color].rank, outgoing[platform][color].platform,
                              MPI_COMM_WORLD, &req);
                } else {
                    MPI_Isend(&NEG_ONE, 1, MPI_INT, outgoing[platform][color].rank, outgoing[platform][color].platform,
                              MPI_COMM_WORLD, &req);
                }
            }
        }
        // Receive trains from connected platforms
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            for (int color : platform_colors[platform]) {
                MPI_Irecv(&recvbuff[platform * num_colors + color], 1, MPI_INT, incoming[platform][color].rank,
                          platform_startidx + platform, MPI_COMM_WORLD, &recvstats[platform * num_colors + color]);
            }
        }

        // Move trains from platform to travelling
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            if (unload_time_left[platform] > 0) { // Remain at 0 to indicate still need to travel
                --unload_time_left[platform];
            }
            if (unload_time_left[platform] == 0 && travel_time_left[platform] <= 0) { // No train on link
                travelling[platform] = trains[platform].front();
                trains[platform].pop_front();
                travel_time_left[platform] = platform_dist[platform];
                unload_time_left[platform] = -1; // -1 to indicate the train already left platform
            }
        }

        // Update trains received from connected platforms
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            for (int color : platform_colors[platform]) {
                MPI_Wait(&recvstats[platform * num_colors + color], MPI_STATUS_IGNORE); // Wait till all trains received
            }
            // Sort received trains by id
            std::sort(recvbuff.begin() + (platform * num_colors), recvbuff.begin() + ((platform + 1) * num_colors));

            for (int i = 0; i < num_colors; ++i) {
                if (recvbuff[platform * num_colors + i] != -1) { // -1 for no train
                    trains[platform].emplace_back(recvbuff[platform * num_colors + i]);
                    recvbuff[platform * num_colors + i] = -1; // Reset for next recv
                }
            }
        }

        // Spawn more trains. Should be after receiving from connected platforms since new train ids will always be higher.
        if (spawns_left) {
            spawns_left = false;
            for (int color = 0; color < num_colors; ++color) {
                for (int platform = 0; platform < num_platforms_rank; ++platform) {
                    if (platform_spawns[platform][color]-- > 0) { // Spawn for current platform
                        trains[platform].emplace_back(train_color.size() + platform_spawns_is_last[platform][color]);
                    }
                }
                int num_spawned = std::min(2, num_trains_vec[color]);
                for (int i = 0; i < num_spawned; ++i) {
                    train_color.emplace_back(color);
                    spawns_left = true;
                }
                num_trains_vec[color] -= num_spawned; // Update for next tick
            }
        }

        // Unload trains
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            if (unload_time_left[platform] >= 0 || trains[platform].size() == 0) {
                continue;
            }
            // Update with new unloading time
            unload_time_left[platform] = pltg[platform].next(trains[platform].front());
        }

        MPI_Barrier(MPI_COMM_WORLD);

        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            if (travel_time_left[platform] > 0) {
                int train = travelling[platform];
                std::cout << tick << ": " <<  (char)colors[train_color[train]] << train << "-" <<
                        station_names[platform_from[platform]] << "->" << station_names[platform_to[platform]] << std::endl;
            }
            for (size_t i = 0; i < trains[platform].size(); ++i) {
                int train = trains[platform][i];
                if (i == 0) {
                    std::cout << tick << ": " << (char)colors[train_color[train]] << train << "-" <<
                        station_names[platform_from[platform]] << "%" << std::endl;
                } else {
                    std::cout << tick << ": " << (char)colors[train_color[train]] << train << "-" <<
                        station_names[platform_from[platform]] << "#" << std::endl;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi_rank == 0) {
            printf("%d, %d\nps:\n", platforms_per_process, platform_startidx);

            for (size_t i = 0; i < platform_from.size(); ++i) {
                printf("from: %d, to: %d\n",platform_from[i], platform_to[i]);
            }
            printf("\npo:\n");
            for (size_t i = 0; i < all_platform_from.size(); ++i) {
                printf("from: %d, to: %d\n",all_platform_from[i], all_platform_to[i]);
            }
            printf("\npc:\n");
            for (auto i : platform_colors) {
                printf("nxt color:\n");
                for (auto j : i) {
                    printf("%d,",j);
                }
                printf("\n");
            }
            printf("\npspawns:\n");
            for (auto i : platform_spawns) {
                printf("nxt spawn:\n");
                for (auto j : i) {
                    printf("%d,",j);
                }
                printf("\n");
            }
            printf("\noutgoing:\n");
            for (auto i : outgoing) {
                printf("nxt outgoing:\n");
                for (auto j : i) {
                    printf("key: %d, platform: %d, rank: %d\n", j.first, j.second.platform, j.second.rank);
                }
            }
            printf("\nincoming:\n");
            for (auto i : incoming) {
                printf("nxt incoming:\n");
                for (auto j : i) {
                    printf("key: %d, platform: %d, rank: %d\n", j.first, j.second.platform, j.second.rank);
                }
            }
            printf("\n");
        }
}