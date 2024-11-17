#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <deque>
#include "platform_key.h"

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

vector<LineColor> colors = {GREEN, YELLOW, BLUE, RED, BROWN, PURPLE, TURQUOISE, PINK, LIME, GREY};

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

// // Find lexicographical ordering of numbers, starting from 0. n must be at least 1.
inline void lexico(unsigned int n, vector<int> &result) {
    unsigned int curr = 1;
    result.emplace_back(0);

    for (unsigned int i = 1; i < n; ++i) { // Get n numbers (excluding first 0)
        result.emplace_back(curr);

        if (curr * 10 < n) { // Endings with 0 are ordered first
            curr *= 10;
        } else{
            if (++curr >= n) { // Exceeded limit
                curr = (curr + 9) / 10; // (curr + 10) / 10 == (curr / 10) + 1
            }
            while (curr % 10 == 0) { // Remove trailing 0s for the correct ordering
                curr /= 10;
            }
        }
    }
}

// Takes the original vector's index -> value, and reverses it into value -> index in reverse vector
inline void reverse_map(vector<int> &original, vector<int> &reverse) {
    for (size_t i = 0; i < original.size(); ++i) {
        reverse[original[i]] = i;
    }
}

// Assign platforms in color order, to maximise self-communication rather than communicating with other processes.
// Each process calculates this to prevent communication of this giant list.
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

// Finds all platforms that precedes this platform, and all platforms that is after this platform.
// These are the platforms this process needs to communicate with (outgoing or incoming).
inline void find_connected_platforms(vector<int> &platform_from, vector<int> &platform_to, vector<unordered_map<int, Platform>> &outgoing,
        vector<unordered_map<int, Platform>> &incoming, int platforms_per_process, vector<vector<int>> &platform_colors,
        vector<int> &all_platform_from, vector<int> &all_platform_to, vector<vector<int>> &station_lines_vec) {

    for (size_t platform = 0; platform < platform_to.size(); ++platform) {
        for (int color : platform_colors[platform]) {
            vector<int> &station_line = station_lines_vec[color];

            // Find index within this current color
            size_t from_idx = std::find(station_line.begin(), station_line.end(), platform_from[platform]) - station_line.begin();
            size_t to_idx = std::find(station_line.begin(), station_line.end(), platform_to[platform]) - station_line.begin();
            int before_from, after_to, result;

            if (from_idx < to_idx) { // Travelling forwards
                if (from_idx == 0) { // Terminal station
                    before_from = station_line[1];
                } else {
                    before_from = station_line[from_idx - 1];
                }
                if (to_idx == station_line.size() - 1) { // Terminal station
                    after_to = station_line[station_line.size() - 2];
                } else {
                    after_to = station_line[to_idx + 1];
                }
                
            } else { // Travelling backwards
                if (from_idx == station_line.size() - 1) {
                    before_from = station_line[station_line.size() - 2]; // Terminal station
                } else {
                    before_from = station_line[from_idx + 1];
                }
                if (to_idx == 0) { // Terminal station
                    after_to = station_line[1];
                } else {
                    after_to = station_line[to_idx - 1];
                }
            }

            result = pair_idx(all_platform_from, all_platform_to, platform_to[platform], after_to); // Find index based on all platforms
            outgoing[platform].emplace(std::piecewise_construct, std::forward_as_tuple(color),
                                       std::forward_as_tuple(result, result / platforms_per_process));

            result = pair_idx(all_platform_from, all_platform_to, before_from, platform_from[platform]); // Find index based on all platforms
            incoming[platform].emplace(std::piecewise_construct, std::forward_as_tuple(color),
                                       std::forward_as_tuple(result, result / platforms_per_process));
        }
    }
}

void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
              const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
              const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
              size_t total_processes, const unordered_map<Platform_Key, vector<int>, Platform_Key_Hasher> &line_failures,
              const unordered_map<int, vector<int>> &train_failures) {

    int num_colors = station_lines.size();
    int num_platforms = 0;
    for (size_t i = 0; i < num_stations; ++i) {
        for (size_t j = 0; j < num_stations; ++j) {
            if (mat[i][j] != 0) {
                ++num_platforms;
            }
        }
    }

    // Number of trains in terms of color indices
    vector<int> num_trains_vec;
    int total_num_trains = 0;
    for (int color = 0; color < num_colors; ++color) {
        num_trains_vec.emplace_back(num_trains.at(colors[color]));
        total_num_trains += num_trains.at(colors[color]);
    }

    // Station lines in terms of color and station_name indices
    vector<vector<int>> station_lines_vec(num_colors);
    for (int color = 0; color < num_colors; ++color) {
        for (const string &station_name: station_lines.at(colors[color])) {
            station_lines_vec[color].emplace_back(std::find(station_names.begin(), station_names.end(), station_name) - station_names.begin());
        }
    }

    int platforms_per_process = (num_platforms + total_processes - 1) / total_processes; // Round up
    int platform_startidx = platforms_per_process * mpi_rank;

    vector<int> platform_from, platform_to, platform_dist; // Only this process' platforms
    vector<vector<int>> platform_colors(platforms_per_process); // Colors for this process' platforms
    vector<vector<int>> platform_spawns(platforms_per_process, vector<int>(num_colors, 0)); // Amount needed to spawn for this process' platforms
    vector<vector<bool>> platform_spawns_is_last(platforms_per_process, vector<bool>(num_colors, false)); // Whether this spawn is the last terminal or first
    vector<int> all_platform_from, all_platform_to; // Sorted based on color, all process' platforms
    all_platform_from.reserve(num_platforms);
    all_platform_to.reserve(num_platforms);

    assign_platforms(platforms_per_process, platform_startidx, platform_from, platform_to, platform_dist, platform_colors, platform_spawns,
            platform_spawns_is_last, all_platform_from, all_platform_to, mat, station_lines_vec, num_trains_vec);
    
    int num_platforms_rank = platform_to.size(); // No. of platforms this rank is in charge of

    vector<vector<int>> platform_failures(num_platforms_rank);
    for (int i = 0; i < num_platforms_rank; ++i) {
        platform_failures[i] = std::move(line_failures.at(Platform_Key(platform_from[i], platform_to[i])));
        // Sort in descending order, so that we can pop from back when done with that tick
        std::sort(platform_failures[i].begin(), platform_failures[i].end(), greater<>());
    }

    // [{pform_0_color_idx_0: target_pform_0, pform_0_color_idx_1: target_pform_1, ...}, {pform_1_...}]
    vector<unordered_map<int, Platform>> outgoing(num_platforms_rank);
    vector<unordered_map<int, Platform>> incoming(num_platforms_rank);

    find_connected_platforms(platform_from, platform_to, outgoing, incoming, platforms_per_process, platform_colors, all_platform_from, all_platform_to, station_lines_vec);

    // Output ordering
    bool spawns_left = true;
    vector<int> output_order; // Each value is a train id, sorted based on color first, then by lexicographical order of the train id
    vector<int> reverse_output_order; // Each index is a train id, and the value is its index in the output_order vector
    
    if (mpi_rank == 0) { // Only rank 0 is printing
        vector<int> sorted_colors(num_colors); // colors (represented by indices) sorted lexicographically
        vector<int> start_free_by_color(num_colors); // Starting index for each color group
        vector<int> next_free_by_color(num_colors); // Next index for a train of this color to be placed at
        std::iota(sorted_colors.begin(), sorted_colors.end(), 0); // Range from 0 to num_colors - 1
        std::sort(sorted_colors.begin(), sorted_colors.end(), [&](int i, int j) {return colors[i] < colors[j];}); // Sort lexicographically

        int cum_sum = 0;
        for (int color : sorted_colors) {
            // Starting index of this color depends on previous colors
            start_free_by_color[color] = cum_sum;
            next_free_by_color[color] = cum_sum;
            cum_sum += num_trains.at(colors[color]);
        }

        reverse_output_order.assign(total_num_trains, 0);
        output_order.reserve(total_num_trains);

        // The vectors here are reused so the vector names do not correspond to the actual items. The comments show what is stored.
        lexico(total_num_trains, output_order); // The first total_num_trains numbers sorted lexicographically (starting from 0)
        reverse_map(output_order, reverse_output_order); // The reverse mapping of the above (used as precomputation for speed)

        int train_id = 0, num_spawned;

        // Spawn all trains from all processes and store the ordering
        while (spawns_left) {
            spawns_left = false;
            for (int color = 0; color < num_colors; ++color) {
                num_spawned = std::min(2, num_trains_vec[color]);
                for (int spawn = 0; spawn < num_spawned; ++spawn) {
                    output_order[next_free_by_color[color]++] = train_id++;
                    spawns_left = true;
                }
                num_trains_vec[color] -= num_spawned;
            }
        }

        // Reset num_trains_vec for use within main loop
        for (int color = 0; color < num_colors; ++color) {
            num_trains_vec[color] = num_trains.at(colors[color]);
        }

        // Sort each group of colors in lexicographical order. next_free_by_color now has the ending index.
        auto lexi_comp = [&](int i, int j) {return reverse_output_order[i] < reverse_output_order[j];};
        for (int i = 0; i < num_colors; ++i) {
            std::sort(output_order.begin() + start_free_by_color[i], output_order.begin() + next_free_by_color[i], lexi_comp);
        }
        reverse_map(output_order, reverse_output_order); // Store the actual reverse mapping now that output_order has been calculated
    }


    vector<PlatformLoadTimeGen> pltg;
    for (int platform = 0; platform < num_platforms_rank; ++platform) {
        pltg.emplace_back(popularities[platform_from[platform]]);
    }

    int NEG_ONE = -1, out_color, idx;
    int int_ticks = (int)ticks, int_num_ticks_to_print = (int)num_ticks_to_print;
    spawns_left = true;

    // Status of trains
    vector<int> travelling(num_platforms_rank);
    vector<int> travel_time_left(num_platforms_rank, -1); // -1 for unused
    vector<int> unload_time_left(num_platforms_rank, -1);
    vector<int> train_color; // Sorted by train id, for all process' trains
    train_color.reserve(ticks << 1); // Max of 2 * ticks trains
    vector<deque<int>> trains(num_platforms_rank); // First element is in platform. The rest are in holding area.

    // MPI communication
    vector<int> sendbuff(num_platforms_rank); // Each platform sends at most 1 train
    vector<int> recvbuff(num_colors * num_platforms_rank, -1); // Each platform can receive trains up to the number of colors available
    vector<MPI_Request> recvstats(num_colors * num_platforms_rank, MPI_REQUEST_NULL);
    MPI_Request unneeded_req; // For unneeded output.
    
    // Each train status will be represented by 3 ints:
    // Tag: Train idx (so only 2 ints sent in buffer)
    // int_1: From station
    // int_2: To station if travelling, -1 if on platform, -2 if in holding
    vector<int> output_sendbuff(total_num_trains << 1);
    vector<int> output_recvbuff;
    if (mpi_rank == 0) {
        output_recvbuff.assign(total_num_trains << 1, -1);
    }

    // 0 for no failure, 1 for link failure, 2 for train failure, 3 for both failure
    vector<int> failure_status(num_platforms_rank);

    for (int tick = 0; tick < int_ticks; ++tick) {
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
                     // Use tags to send to the correct platform, since each rank may have multiple platforms
                    MPI_Isend(&sendbuff[platform], 1, MPI_INT, outgoing[platform][color].rank, outgoing[platform][color].platform,
                              MPI_COMM_WORLD, &unneeded_req); // Request not needed since waiting at receiving end
                } else {
                    MPI_Isend(&NEG_ONE, 1, MPI_INT, outgoing[platform][color].rank, outgoing[platform][color].platform,
                              MPI_COMM_WORLD, &unneeded_req);
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
            MPI_Waitall(num_colors, &recvstats[platform * num_colors], MPI_STATUSES_IGNORE);

            // Sort received trains by id
            std::sort(recvbuff.begin() + (platform * num_colors), recvbuff.begin() + ((platform + 1) * num_colors));

            for (int i = 0; i < num_colors; ++i) {
                if (recvbuff[platform * num_colors + i] != -1) { // -1 for no train
                    trains[platform].emplace_back(recvbuff[platform * num_colors + i]); // Add train to queue
                    recvbuff[platform * num_colors + i] = -1; // Reset for next recv
                }
            }
        }

        // Spawn more trains. Should be added to queue after receiving from connected platforms since new train ids will always be higher.
        if (spawns_left) {
            spawns_left = false;
            for (int color = 0; color < num_colors; ++color) {
                for (int platform = 0; platform < num_platforms_rank; ++platform) {
                    if (platform_spawns[platform][color]-- > 0) { // Spawn for current platform
                        // train_color has all trains spawned till now (its size is next train idx)
                        trains[platform].emplace_back(train_color.size() + platform_spawns_is_last[platform][color]);
                    }
                }
                int num_spawned = std::min(2, num_trains_vec[color]);
                for (int i = 0; i < num_spawned; ++i) {
                    train_color.emplace_back(color); // Store color information for all process' trains
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

        // Link failures
        for (int platform = 0; platform < num_platforms_rank; ++platform) {
            if (platform_failures[platform].back() == tick) {
                platform_failures[platform].pop_back();
                if (failure_status[platform] == 1) { // Same failure
                    continue;
                } else if (failure_status[platform] < 3) {
                    ++failure_status[platform];
                    travel_time_left[platform] << 1; // Will not affect if <= 0, since it will still be <= 0 after *2
                }
            }
        }

        if (tick >= int_ticks - int_num_ticks_to_print) {
            idx = 0; // Reset
            MPI_Barrier(MPI_COMM_WORLD); // Needed since rank 0 could still be receiving MPI msgs

            // Send train details for output
            for (int platform = 0; platform < num_platforms_rank; ++platform) {
                if (travel_time_left[platform] > 0) {
                    output_sendbuff[idx] = platform_from[platform];
                    output_sendbuff[idx + 1] = platform_to[platform];
                    MPI_Isend(&output_sendbuff[idx], 2, MPI_INT, 0, travelling[platform], MPI_COMM_WORLD, &unneeded_req);
                    idx += 2; // Random order on sender's side in output_sendbuff
                }
                
                for (size_t i = 0; i < trains[platform].size(); ++i) {
                    output_sendbuff[idx] = platform_from[platform];
                    if (i == 0) {
                        output_sendbuff[idx + 1] = -1; // Platform
                    } else {
                        output_sendbuff[idx + 1] = -2; // Holding area
                    }
                    MPI_Isend(&output_sendbuff[idx], 2, MPI_INT, 0, trains[platform][i], MPI_COMM_WORLD, &unneeded_req);
                    idx += 2; // Random order on sender's side in output_sendbuff
                }
            }

            // Output train details
            if (mpi_rank == 0) {
                // train_color has the exact number of trains we have in this tick
                for (size_t train = 0; train < train_color.size(); ++train) {
                    // Place them in specified order so don't need to sort
                    MPI_Recv(&output_recvbuff[reverse_output_order[train] << 1], 2, MPI_INT, MPI_ANY_SOURCE, train, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                printf("%d:", tick);

                // No need to clear output_recvbuff since all trains in the next tick will be in the same position.
                // All -1 will remain until that train spawns.
                for (int train = 0; train < total_num_trains; ++train) {
                    idx = train << 1;
                    if (output_recvbuff[idx] == -1) { // From station empty. Train not yet spawned.
                        continue;
                    }

                    if (output_recvbuff[idx + 1] == -1) { // Platform
                        printf(" %c%d-%s%%", colors[train_color[output_order[train]]], output_order[train],
                            station_names[output_recvbuff[idx]].c_str());
                    } else if (output_recvbuff[idx + 1] == -2) { // Holding area
                        printf(" %c%d-%s#", colors[train_color[output_order[train]]], output_order[train],
                            station_names[output_recvbuff[idx]].c_str());
                    } else { // Travelling
                        printf(" %c%d-%s->%s", colors[train_color[output_order[train]]], output_order[train],
                            station_names[output_recvbuff[idx]].c_str(), station_names[output_recvbuff[idx + 1]].c_str());
                    }
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronise all ranks
    }
}