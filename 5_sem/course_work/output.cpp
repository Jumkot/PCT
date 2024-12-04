#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double read_serial_time(const std::string& filename)
{
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

void process_files(
        const std::vector<std::string>& files,
        const std::string& output_filename,
        double serial_time)
{
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Не удалось открыть файл " << output_filename
                  << " для записи!\n";
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

int main()
{
    double serial_time1 = read_serial_time("results/serial/serial1.txt");
    double serial_time2 = read_serial_time("results/serial/serial2.txt");
    double serial_time4 = read_serial_time("results/serial/serial4.txt");

    if (serial_time1 < 0 || serial_time2 < 0 || serial_time4 < 0) {
        return 1;
    }

    std::vector<std::string> inverse1
            = {"results/1000/inverse1_4.txt",
               "results/1000/inverse1_8.txt",
               "results/1000/inverse1_12.txt",
               "results/1000/inverse1_16.txt",
               "results/1000/inverse1_20.txt",
               "results/1000/inverse1_24.txt",
               "results/1000/inverse1_28.txt",
               "results/1000/inverse1_32.txt"};
    process_files(inverse1, "results/inverse1.dat", serial_time1);

    std::vector<std::string> inverse2
            = {"results/2500/inverse2_4.txt",
                "results/2500/inverse2_8.txt",
                "results/2500/inverse2_12.txt",
                "results/2500/inverse2_16.txt",
                "results/2500/inverse2_20.txt",
                "results/2500/inverse2_24.txt",
                "results/2500/inverse2_28.txt",
                "results/2500/inverse2_32.txt"};
    process_files(inverse2, "results/inverse2.dat", serial_time2);

    std::vector<std::string> inverse4
            = {"results/4000/inverse4_4.txt",
               "results/4000/inverse4_8.txt",
               "results/4000/inverse4_12.txt",
               "results/4000/inverse4_16.txt",
               "results/4000/inverse4_20.txt",
               "results/4000/inverse4_24.txt",
               "results/4000/inverse4_28.txt",
               "results/4000/inverse4_32.txt"};
    process_files(inverse4, "results/inverse4.dat", serial_time4);

    return 0;
}
