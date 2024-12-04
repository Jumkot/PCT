#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double read_serial_time(const std::string& filename) {
    std::ifstream serial_file(filename);
    if (!serial_file.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << "!\n";
        return -1.0;
    }

    double serial_time;
    serial_file >> serial_time;
    serial_file.close();
    
    return serial_time;
}

void process_files(const std::vector<std::string>& files, const std::string& output_filename, double serial_time) {
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Не удалось открыть файл " << output_filename << " для записи!\n";
        return;
    }

    for (const auto& filename : files) {
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            std::cerr << "Не удалось открыть файл " << filename << "!\n";
            continue;
        }

        int num_processes;
        double parallel_time;
        while (input_file >> num_processes >> parallel_time) {
            // Вычисление ускорения
            double speedup = serial_time / parallel_time;
            output_file << num_processes << "\t" << speedup << "\n";
        }

        input_file.close();
    }

    output_file.close();
}

int main() {
    double serial_time1 = read_serial_time("results/serial1.txt");
    double serial_time3 = read_serial_time("results/serial3.txt");

    if (serial_time1 < 0 || serial_time3 < 0) {
        return 1;
    }

    std::vector<std::string> inverse1 = {
        "results/inverse1_4.txt",
        "results/inverse1_8.txt",
        "results/inverse1_16.txt",
        "results/inverse1_32.txt"
    };
    process_files(inverse1, "results/inverse1.dat", serial_time1);

    // std::vector<std::string> laplace4 = {
    //     "results/laplace4_8.txt",
    //     "results/laplace4_16.txt",
    //     "results/laplace4_32.txt"
    // };
    // process_files(laplace4, "results/laplace4.dat", serial_time3);

    return 0;
}
