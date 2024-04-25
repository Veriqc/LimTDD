#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include <filesystem>
#include <iostream>
#include <string>

using namespace dd;
xt::xarray<dd::ComplexValue> stateToArray(const BasisStates& state){
    switch (state) {
                    case BasisStates::zero:
                        return {complex_one,complex_zero};
                    case BasisStates::one:
                        return {complex_zero,complex_one};
                    case BasisStates::plus:
                        return {complex_SQRT2_2,complex_SQRT2_2};
                    case BasisStates::minus:
                        return {complex_SQRT2_2,complex_mSQRT2_2};
                    case BasisStates::right:
                        return {complex_SQRT2_2,complex_iSQRT2_2};
                    case BasisStates::left:
                        return {complex_SQRT2_2,complex_miSQRT2_2};
    }
}
TDD makezero(int n, dd::Package<>* ddpackage, std::vector<BasisStates> states) {
    TensorNetwork tn;
    if(n> states.size()){
        throw std::invalid_argument("wrong qubit number");
    }
    for(int i=0; i < n; i++){
        xt::xarray<dd::ComplexValue> array = stateToArray(states[i]);
        Tensor temp = Tensor(array,{{"x"+std::to_string(i)+"_0",0}});
        tn.add_ts(temp);
    }
    return tn.cont(ddpackage);
}
TDD cont(dd::TensorNetwork* tn,dd::Package<>* ddpackage, int n,bool simulate,const std::vector<BasisStates>& states,bool release = true) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tn->tensors.size() == 0) {
        throw std::runtime_error("null tensor network");
    }

    clock_t start,end;
    start = clock();
    TDD res_dd = simulate ? makezero(n, ddpackage, states) : tn->tensors[0].to_tdd(ddpackage);
    ddpackage->incRef(res_dd.e);
    unsigned int MAX_NODE = ddpackage->size(res_dd.e);

    // The loop starts from 0 if simulating, 1 otherwise.
    for (size_t i = simulate ? 0 : 1; i < tn->tensors.size(); ++i) {
        try {
            TDD temp_dd = ddpackage->cont(res_dd, tn->tensors[i].to_tdd(ddpackage));
            if (release) {
                ddpackage->incRef(temp_dd.e);
                ddpackage->decRef(res_dd.e);
                ddpackage->garbageCollect();
            }
            res_dd = temp_dd;
            MAX_NODE = std::max(MAX_NODE, ddpackage->size(res_dd.e));
        } catch (...) {
            std::exception_ptr p = std::current_exception();
            std::clog << (p ? p.__cxa_exception_type()->name() : "null ") << std::endl;
        }
    }
    end = clock();
    std::cout<<"time: " << double(end-start)/CLOCKS_PER_SEC << "s" <<std::endl;
    std::cout<<"MAX node: " << MAX_NODE  <<std::endl;

    return res_dd;
};


BasisStates charToBasisState(char c) {
    static const std::map<char, BasisStates> stateMap = {
        {'0', BasisStates::zero},
        {'1', BasisStates::one},
        {'+', BasisStates::plus},
        {'-', BasisStates::minus},
        {'>', BasisStates::right},
        {'<', BasisStates::left}
    };

    auto it = stateMap.find(c);
    if (it != stateMap.end()) {
        return it->second;
    } else {
        throw std::invalid_argument("Invalid basis state character: " + std::string(1, c));
    }
}
std::vector<BasisStates> stringToBasisStates(const std::string& states) {
    std::vector<BasisStates> basisStates;
    for (char c : states) {
        basisStates.push_back(charToBasisState(c));
    }
    return basisStates;
}
int main(int argc, char *argv[]) {
    // filename, initial state
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <number>\n";
        return 1;
    }
    
    
    std::string filename = argv[1];
	std::cout << filename << std::endl;
    
    std::ifstream fileStream(filename);
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    std::string fileContent = buffer.str();
    fileStream.close();

    // Use the file content with QuantumComputation::fromQASM
    const auto qc = qc::QuantumComputation::fromQASM(fileContent);
    std::shared_ptr<qc::QuantumComputation> QC = std::make_shared<qc::QuantumComputation>(std::move(qc));
    auto ddPack = std::make_shared<dd::Package<>>(3*QC->getNqubits());
    auto tn = cir_2_tn(QC,ddPack);

    bool simulate = false;
    std::vector<BasisStates> initialStates;
    if(argc > 2){
        simulate = true;
         try {
            initialStates = stringToBasisStates(argv[2]);
            std::cout << "Initial states vector size: " << initialStates.size() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
    std::cout <<"simulate:" << simulate << std::endl;

	dd::TDD tdd = cont(&tn,ddPack.get(),QC->getNqubits(),simulate, initialStates);
    
    std::cout<<"final node: " << ddPack->size(tdd.e) <<std::endl;
    return 0;
}
