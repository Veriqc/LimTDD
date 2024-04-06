#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include <filesystem>
#include <iostream>
#include <string>

using namespace dd;
TDD makezero(int n, dd::Package<>* ddpackage) {
    TensorNetwork tn;
    for(int i=0; i < n; i++){
        xt::xarray<dd::ComplexValue> zero = {complex_one,complex_zero};
        Tensor temp = Tensor(zero,{{"x"+std::to_string(i)+"_0",0}});
        tn.add_ts(temp);
    }
    return tn.cont(ddpackage);
}
TDD cont(dd::TensorNetwork* tn,dd::Package<>* ddpackage, int n,bool release = true) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tn->tensors.size() == 0) {
        throw std::runtime_error("null tensor network");
    }

    clock_t start,end;
    TDD res_dd = makezero(n,ddpackage);

    ddpackage->incRef(res_dd.e);
    TDD temp_dd;
    unsigned int MAX_NODE = ddpackage->size(res_dd.e);
    start = clock();
    for (int i=0; i < tn->tensors.size(); ++i) {
        auto tensor = tn->tensors[i];
        try{
            temp_dd = ddpackage->cont(res_dd, tensor.to_tdd(ddpackage));
            if (release) {
                ddpackage->incRef(temp_dd.e);
                ddpackage->decRef(res_dd.e);
                ddpackage->garbageCollect();
            }
            res_dd = temp_dd;
            MAX_NODE = std::max(MAX_NODE, ddpackage->size(res_dd.e));
        }
        catch(...){
            std::exception_ptr p = std::current_exception();
            std::clog <<(p ? p.__cxa_exception_type()->name() : "null ") << std::endl;
        }
    }
    end = clock();
    std::cout<<"time: " << double(end-start)/CLOCKS_PER_SEC << "s" <<std::endl;
    std::cout<<"MAX node: " << MAX_NODE  <<std::endl;

    return res_dd;
};
int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <number>\n";
        return 1;
    }

    std::string path = std::string(PROJECT_SOURCE_DIR)+"/Benchmarks/clifford_s/";
    std::string file_name = argv[1]+std::string("time") + ".qasm";
	std::cout << path+file_name << std::endl;
    
    std::ifstream fileStream(path+file_name);
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    std::string fileContent = buffer.str();
    fileStream.close();

    // Use the file content with QuantumComputation::fromQASM
    const auto qc = qc::QuantumComputation::fromQASM(fileContent);
    std::shared_ptr<qc::QuantumComputation> QC = std::make_shared<qc::QuantumComputation>(std::move(qc));
    auto ddPack = std::make_shared<dd::Package<>>(3*QC->getNqubits());
    auto tn = cir_2_tn(QC,ddPack);

	dd::TDD tdd = cont(&tn,ddPack.get(),QC->getNqubits());
    
    std::cout<<"final node: " << ddPack->size(tdd.e) <<std::endl;
    return 0;
}
