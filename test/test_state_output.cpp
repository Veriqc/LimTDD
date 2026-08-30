#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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
    throw std::invalid_argument("Unsupported basis state");
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

Edge<mNode> sliceStateEdge(const Edge<mNode>& e, const int x, const int c, dd::Package<>* ddpackage) {
    assert(e.w != Complex::zero);
    if (e.p->v == -1 || e.p->v < x) {
        return e;
    }
    if (e.p->v != x) {
        throw std::runtime_error("sliceStateEdge only supports direct slicing on the current variable");
    }

    if (e.p->v != e.map->level) {
        auto temp = e.p->e[c];
        if (temp.w != Complex::zero) {
            temp.w = ddpackage->cn.mulCached(temp.w, e.w);
            auto mr = ddpackage->mapmul(e.map, temp.map); temp.map = mr.map;
            ddpackage->cn.mul(temp.w, temp.w, ddpackage->cn.getTemporary(cos(mr.phase * rotate_angle), sin(mr.phase * rotate_angle)));
        }
        return temp;
    }

    if (e.map->x == 0) {
        auto temp = e.p->e[c];
        if (temp.w != Complex::zero) {
            temp.w = ddpackage->cn.mulCached(temp.w, e.w);
            auto mr = ddpackage->mapmul(e.map->father, temp.map); temp.map = mr.map;
            ddpackage->cn.mul(temp.w, temp.w, ddpackage->cn.getTemporary(cos(mr.phase * rotate_angle), sin(mr.phase * rotate_angle)));
            if (c == 1) {
                ddpackage->cn.mul(temp.w, temp.w, ddpackage->cn.getTemporary(cos(e.map->rotate * rotate_angle), sin(e.map->rotate * rotate_angle)));
            }
        }
        return temp;
    }

    auto temp = e.p->e[1 - c];
    if (temp.w != Complex::zero) {
        temp.w = ddpackage->cn.mulCached(temp.w, e.w);
        auto mr = ddpackage->mapmul(e.map->father, temp.map); temp.map = mr.map;
        ddpackage->cn.mul(temp.w, temp.w, ddpackage->cn.getTemporary(cos(mr.phase * rotate_angle), sin(mr.phase * rotate_angle)));
        if (c == 0) {
            ddpackage->cn.mul(temp.w, temp.w, ddpackage->cn.getTemporary(cos(e.map->rotate * rotate_angle), sin(e.map->rotate * rotate_angle)));
        }
    }
    return temp;
}

Complex amplitudeForBitstring(const TDD& tdd, const std::string& basisState, dd::Package<>* ddpackage) {
    auto edge = tdd.e;
    for (const auto bitChar : basisState) {
        if (edge.p->v == -1) {
            break;
        }
        const auto bit = bitChar == '1' ? 1 : 0;
        edge = sliceStateEdge(edge, edge.p->v, bit, ddpackage);
    }
    return edge.w;
}

std::vector<std::pair<std::string, Complex>> tddToStateVector(const TDD& tdd, dd::Package<>* ddpackage, const std::size_t qubitCount) {
    std::vector<std::pair<std::string, Complex>> stateVector;
    const auto basisCount = static_cast<std::size_t>(1ULL << qubitCount);
    stateVector.reserve(basisCount);
    for (std::size_t basisIndex = 0; basisIndex < basisCount; ++basisIndex) {
        auto basisState = std::bitset<64>(basisIndex).to_string();
        basisState = basisState.substr(64 - qubitCount);
        stateVector.emplace_back(basisState, amplitudeForBitstring(tdd, basisState, ddpackage));
    }
    return stateVector;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <qasm-file> [initial-state] [basis-state]" << std::endl;
        return 1;
    }

    const std::string filename = argv[1];
    std::ifstream fileStream(filename);
    if (!fileStream.is_open()) {
        std::cerr << "Unable to open QASM file: " << filename << std::endl;
        return 2;
    }
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    const auto qc = qc::QuantumComputation::fromQASM(buffer.str());
    auto QC = std::make_shared<qc::QuantumComputation>(std::move(qc));
    auto ddPack = std::make_shared<dd::Package<>>(3 * QC->getNqubits());
    auto tn = cir_2_tn(QC, ddPack);

    bool simulate = false;
    std::vector<BasisStates> initialStates;
    std::string basisStateQuery;
    if (argc > 2) {
        simulate = true;
        try {
            initialStates = stringToBasisStates(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Invalid initial state: " << e.what() << std::endl;
            return 3;
        }
    }
    if (argc > 3) {
        basisStateQuery = argv[3];
        if (basisStateQuery.size() != QC->getNqubits()) {
            std::cerr << "Invalid basis-state length: expected " << QC->getNqubits()
                      << ", got " << basisStateQuery.size() << std::endl;
            return 4;
        }
        for (const auto bit : basisStateQuery) {
            if (bit != '0' && bit != '1') {
                std::cerr << "Invalid basis-state bit: " << bit << std::endl;
                return 5;
            }
        }
    }

    auto tdd = cont(&tn, ddPack.get(), QC->getNqubits(), simulate, initialStates, true);

    if (!basisStateQuery.empty()) {
        const auto amplitude = amplitudeForBitstring(tdd, basisStateQuery, ddPack.get());
        std::cout << "CPP_LIMTDD_STATE_BEGIN" << std::endl;
        std::cout << "qubits\t" << QC->getNqubits() << std::endl;
        std::cout << std::setprecision(17);
        std::cout << "STATE\t" << basisStateQuery << "\t"
                  << CTEntry::val(amplitude.r) << "\t"
                  << CTEntry::val(amplitude.i) << std::endl;
        std::cout << "CPP_LIMTDD_STATE_END" << std::endl;
        return 0;
    }

    const auto stateVector = tddToStateVector(tdd, ddPack.get(), QC->getNqubits());

    std::cout << "CPP_LIMTDD_STATE_BEGIN" << std::endl;
    std::cout << "qubits\t" << QC->getNqubits() << std::endl;
    std::cout << std::setprecision(17);
    for (const auto& [basisState, amplitude] : stateVector) {
        std::cout << "STATE\t" << basisState << "\t"
                  << CTEntry::val(amplitude.r) << "\t"
                  << CTEntry::val(amplitude.i) << std::endl;
    }
    std::cout << "CPP_LIMTDD_STATE_END" << std::endl;
    return 0;
}