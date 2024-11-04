#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main()
{
    // gauss28
    std::vector<std::string> gauss28
            = {"gauss28_4.txt",
               "gauss28_8.txt",
               "gauss28_16.txt",
               "gauss28_32.txt"};

    std::ofstream output_g28("gauss28.dat");
    if (!output_g28.is_open()) {
        std::cerr << "Не удалось открыть файл gauss28.dat для записи!\n";
        return 1;
    }

    for (const auto& filename : gauss28) {
        std::ifstream input_g28(filename);
        if (!input_g28.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        std::string line;
        if (std::getline(input_g28, line)) {
            output_g28 << line << "\n";
        }

        input_g28.close();
    }

    output_g28.close();

    // gauss45
    std::vector<std::string> gauss45
            = {"gauss45_4.txt",
               "gauss45_8.txt",
               "gauss45_16.txt",
               "gauss45_32.txt"};

    std::ofstream output_g45("gauss45.dat");
    if (!output_g45.is_open()) {
        std::cerr << "Не удалось открыть файл gauss45.dat для записи!\n";
        return 1;
    }

    for (const auto& filename : gauss45) {
        std::ifstream input_g45(filename);
        if (!input_g45.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        std::string line;
        if (std::getline(input_g45, line)) {
            output_g45 << line << "\n";
        }

        output_g45.close();
    }

    output_g45.close();

    return 0;
}
