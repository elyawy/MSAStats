#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "MsaStatsCalculator.h"
#include "summary_stats_defs.h"
namespace py = pybind11;

PYBIND11_MODULE(msastats, m) {

    m.def("calculate_msa_stats", [](const std::vector<std::string> &msa) {
        MsaStatsCalculator msacalc(msa);
        msacalc.recomputeStats();
        return msacalc.getStatVec();
    });

    m.def("calculate_fasta_stats", [](std::string fasta_path) {
        MsaStatsCalculator msacalc(fasta_path);
        msacalc.recomputeStats();
        return msacalc.getStatVec();
    });

    m.def("stats_names", []() {

        std::vector<std::string> names;
        for (auto name: stat_type::_names()){
            names.push_back(name);
        }

        return names;
    });

}