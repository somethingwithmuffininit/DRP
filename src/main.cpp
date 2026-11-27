#include "mixing.hpp"
#include <charconv>
#include <iostream>
#include <thread>

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout << "mixbench run CONFIG.json [--output FRESH_DIR] [--jobs N]\n";
            return 0;
        }
        mix::require(argc >= 3 && std::string(argv[1]) == "run",
                     "usage: mixbench run CONFIG.json [--output DIR] [--jobs N]");
        std::filesystem::path config = argv[2],
                              output = std::filesystem::path("results") / config.stem();
        int jobs = std::max(1u, std::thread::hardware_concurrency());
        for (int i = 3; i < argc; ++i) {
            std::string flag = argv[i];
            mix::require(++i < argc, "missing value for " + flag);
            std::string value = argv[i];
            if (flag == "--output")
                output = value;
            else if (flag == "--jobs") {
                auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), jobs);
                mix::require(ec == std::errc{} && end == value.data() + value.size() && jobs > 0,
                             "invalid job count");
            } else
                throw std::invalid_argument("unknown flag: " + flag);
        }
        mix::run(config, output, jobs);
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
