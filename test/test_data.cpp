#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include "dd/Tensor.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

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
            // std::clog << (p ? p.__cxa_exception_type()->name() : "null ") << std::endl;
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

bool envFlagEnabled(const char* name) {
    const auto* value = std::getenv(name);
    return value != nullptr && std::string_view(value) == "1";
}

int runXarraySelftest() {
    dd::ComplexValue one = {1, 0};
    dd::ComplexValue zero = {0, 0};
    dd::ComplexValue two = {2, 0};
    dd::ComplexValue three = {3, 0};

    auto ddPack = std::make_shared<dd::Package<>>(10);
    ddPack->varOrder = {{"x0", 0}, {"y0", 1}, {"x1", 2}, {"y1", 3}};

    xt::xarray<dd::ComplexValue> tensorCnot = {
        {{{one, zero}, {zero, one}}, {{zero, zero}, {zero, zero}}},
        {{{zero, zero}, {zero, zero}}, {{zero, one}, {one, zero}}},
    };
    xt::xarray<dd::ComplexValue> permutedTensorCnot = {
        {{{zero, zero}, {zero, zero}}, {{zero, zero}, {zero, zero}}},
        {{{zero, zero}, {zero, zero}}, {{zero, zero}, {zero, zero}}},
    };
    for (std::size_t a = 0; a < 2; ++a) {
        for (std::size_t b = 0; b < 2; ++b) {
            for (std::size_t c = 0; c < 2; ++c) {
                for (std::size_t d = 0; d < 2; ++d) {
                    permutedTensorCnot(a, b, c, d) = tensorCnot(c, d, a, b);
                }
            }
        }
    }

    std::vector<dd::Index> tensorIndices = {{"x0", 0}, {"y0", 0}, {"x1", 0}, {"y1", 0}};
    std::vector<dd::Index> permutedIndices = {{"x1", 0}, {"y1", 0}, {"x0", 0}, {"y0", 0}};

    auto tensorTdd = dd::Tensor(tensorCnot, tensorIndices, "tensor_cnot").to_tdd(ddPack.get());
    auto permutedTdd = dd::Tensor(permutedTensorCnot, permutedIndices, "permuted_tensor_cnot").to_tdd(ddPack.get());

    const bool orderEqual = tensorTdd.e == permutedTdd.e;
    std::cout << "xarray_selftest.equal: " << orderEqual << std::endl;
    std::cout << "xarray_selftest.tensor_nodes: " << ddPack->size(tensorTdd.e) << std::endl;
    std::cout << "xarray_selftest.permuted_nodes: " << ddPack->size(permutedTdd.e) << std::endl;

    auto hyperPack = std::make_shared<dd::Package<>>(10);
    hyperPack->varOrder = {{"x0", 0}, {"x1", 2}, {"y1", 3}};
    xt::xarray<dd::ComplexValue> repeatedIndexTensor = tensorCnot;
    xt::xarray<dd::ComplexValue> offDiagonalPerturbed = tensorCnot;
    offDiagonalPerturbed(0, 1, 0, 0) = two;
    offDiagonalPerturbed(0, 1, 1, 1) = three;
    offDiagonalPerturbed(1, 0, 0, 1) = three;
    offDiagonalPerturbed(1, 0, 1, 0) = two;

    std::vector<dd::Index> repeatedIndices = {{"x0", 0}, {"x0", 1}, {"x1", 0}, {"y1", 0}};
    auto repeatedTdd = dd::Tensor(repeatedIndexTensor, repeatedIndices, "repeated_index_base").to_tdd(hyperPack.get());
    auto perturbedTdd = dd::Tensor(offDiagonalPerturbed, repeatedIndices, "repeated_index_perturbed").to_tdd(hyperPack.get());
    const bool repeatedIndexEqual = repeatedTdd.e == perturbedTdd.e;
    std::cout << "xarray_selftest.repeated_index_equal: " << repeatedIndexEqual << std::endl;
    std::cout << "xarray_selftest.repeated_index_nodes: " << hyperPack->size(repeatedTdd.e) << std::endl;
    std::cout << "xarray_selftest.repeated_index_perturbed_nodes: " << hyperPack->size(perturbedTdd.e) << std::endl;

    return (orderEqual && repeatedIndexEqual) ? 0 : 2;
}

int main(int argc, char *argv[]) {
    if (envFlagEnabled("LIMTDD_XARRAY_SELFTEST")) {
        return runXarraySelftest();
    }

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
    ddPack->enableRegressionDiagnostics = envFlagEnabled("LIMTDD_REGRESSION_DIAG");
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
    // dd::export2Dot(tdd.e,"test",true,true);
    
    std::cout<<"final node: " << ddPack->size(tdd.e) <<std::endl;
    if (ddPack->enableRegressionDiagnostics) {
        ddPack->printRegressionDiagnostics();
    }
    return 0;
}
