#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
// #include "gtest/gtest.h"
#include <filesystem>
#include <iostream>
#include <string>

#include <ctime>

int main(){
    const std::string filePath = std::string(PROJECT_SOURCE_DIR)+"/Benchmarks/test.qasm";

    // Open the file
    std::ifstream fileStream(filePath);
    if (!fileStream.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return -1; // or handle error appropriately
    }

    // Read the file content into a string
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    std::string fileContent = buffer.str();

    // Close the file (optional here since the ifstream will close itself on destruction)
    fileStream.close();

    // Use the file content with QuantumComputation::fromQASM
    const auto qc = qc::QuantumComputation::fromQASM(fileContent);
    std::shared_ptr<qc::QuantumComputation> QC = std::make_shared<qc::QuantumComputation>(std::move(qc));
    auto ddPack = std::make_shared<dd::Package<>>(3*QC->getNqubits());
    auto ts = cir_2_tn(QC,ddPack);
    std::cout << ts.infor() << std::endl;
    std::clock_t start = std::clock();  // 获取开始时间
    auto tdd = ts.cont(ddPack.get());
    std::cout<<"Key: "<<tdd.e.p->v<<std::endl;
    std::cout<<"final node: " << ddPack->size(tdd.e) <<std::endl;

    std::clock_t end = std::clock();  // 获取结束时间

    // 计算运行时间（单位：秒）
    double duration = double(end - start) / CLOCKS_PER_SEC;

    std::cout << "Time: " << duration << "s" << std::endl;

    ddPack->statistics();
    // double aa = 3.141592653589793238462643383279502884197169399375105820974;
    // double bb = std::cos(aa*2/32768);
    // double cc = std::sin(aa*2/32768);
    // std::cout << bb << std::endl;
    // std::cout << cc << std::endl;
    // std::cout << std::cos(aa*2/32768)*std::cos(aa*2/32768)+std::sin(aa*2/32768)*std::sin(aa*2/32768)<<std::endl;
    // std::cout << bb*bb+cc*cc<<std::endl;
    dd::export2Dot(tdd.e, "test");
    int number;
    std::cout << "Enter an integer: ";
    std::cin >> number;
    system("pause");
    system("pause");
    system("pause");
    return 0;
}