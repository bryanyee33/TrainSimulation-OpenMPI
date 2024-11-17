#include <string>
#include <vector>
#include <algorithm>

struct Platform_Key {
    int from;
    int to;

    Platform_Key(std::string from_str, std::string to_str, const std::vector<std::string> &station_names) {
        from = std::find(station_names.begin(), station_names.end(), from_str) - station_names.begin();
        to = std::find(station_names.begin(), station_names.end(), to_str) - station_names.begin();
    }
    Platform_Key(int from, int to) : from(from), to(to) {}

    bool operator==(const Platform_Key &other) const {
        return from == other.from && to == other.to;
    }
};

struct Platform_Key_Hasher {
    size_t operator()(const Platform_Key &p) const {
        return p.from ^ p.to;
    }
};