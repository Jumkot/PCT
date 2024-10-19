#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main() {
    // midpoints
    std::vector<std::string> midpoints = {
        "midpoint_4.txt", "midpoint_8.txt", "midpoint_16.txt",
        "midpoint_20.txt", "midpoint_28.txt", "midpoint_32.txt"
    };

    std::ofstream output_mp("midpoints.dat");
    if (!output_mp.is_open()) {
        std::cerr << "Не удалось открыть файл midpoints.dat для записи!\n";
        return 1;
    }

    for (const auto& filename : midpoints) {
        std::ifstream input_mp(filename);
        if (!input_mp.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        std::string line;
        if (std::getline(input_mp, line)) {
            output_mp << line << "\n";
        }

        input_mp.close();
    }

    output_mp.close();
    
    // montecarlos7
    std::vector<std::string> montecarlos7 = {
        "montecarlo7_4.txt", "montecarlo7_8.txt", "montecarlo7_16.txt",
        "montecarlo7_20.txt", "montecarlo7_28.txt", "montecarlo7_32.txt"
    };

    std::ofstream output_mc7("montecarlos7.dat");
    if (!output_mc7.is_open()) {
        std::cerr << "Не удалось открыть файл montecarlos7.dat для записи!\n";
        return 1;
    }

    for (const auto& filename : montecarlos7) {
        std::ifstream input_mc7(filename);
        if (!input_mc7.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        std::string line;
        if (std::getline(input_mc7, line)) {
            output_mc7 << line << "\n";
        }

        input_mc7.close();
    }

    output_mc7.close();

    // montecarlos8
    std::vector<std::string> montecarlos8 = {
        "montecarlo8_4.txt", "montecarlo8_8.txt", "montecarlo8_16.txt",
        "montecarlo8_20.txt", "montecarlo8_28.txt", "montecarlo8_32.txt"
    };

    std::ofstream output_mc8("montecarlos8.dat");
    if (!output_mc8.is_open()) {
        std::cerr << "Не удалось открыть файл montecarlos8.dat для записи!\n";
        return 1;
    }

    for (const auto& filename : montecarlos8) {
        std::ifstream input_mc8(filename);
        if (!input_mc8.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        std::string line;
        if (std::getline(input_mc8, line)) {
            output_mc8 << line << "\n";
        }

        output_mc8.close();
    }

    output_mc8.close();

    return 0;
}
