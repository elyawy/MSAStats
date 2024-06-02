#include <iostream>
#include <cassert>

#include "../src/MSA.h"

int main() {

    
    MsaStatsCalculator msa("msas/dummy_msa.fasta");
    msa.recomputeStats();

    assert(abs(msa.getStatValByType(stat_type::AVG_GAP_SIZE) - 2.71) < 1e-2);

    assert(msa.getStatValByType(stat_type::MSA_LEN) == 37);
    assert(msa.getStatValByType(stat_type::MSA_MAX_LEN) == 28);
    assert(msa.getStatValByType(stat_type::MSA_MIN_LEN) == 15);

    assert(msa.getStatValByType(stat_type::TOT_NUM_GAPS) == 24);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_ONE) == 6);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_TWO) == 6);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_THREE) == 6);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_AT_LEAST_FOUR) == 6);

    assert(abs(msa.getStatValByType(stat_type::AVG_UNIQUE_GAP_SIZE) - 2.67) < 1e-2);
    assert(msa.getStatValByType(stat_type::TOT_NUM_UNIQUE_GAPS) == 12);

    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_ONE_POS_1_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_ONE_POS_2_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_ONE_POS_N_MINUS_1_GAPS) == 1);

    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_TWO_POS_1_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_TWO_POS_2_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_TWO_POS_N_MINUS_1_GAPS) == 1);

    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_THREE_POS_1_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_THREE_POS_2_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_THREE_POS_N_MINUS_1_GAPS) == 1);

    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_AT_LEAST_FOUR_POS_1_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_AT_LEAST_FOUR_POS_2_GAPS) == 1);
    assert(msa.getStatValByType(stat_type::NUM_GAPS_LEN_AT_LEAST_FOUR_POS_N_MINUS_1_GAPS) == 1);

    assert(msa.getStatValByType(stat_type::MSA_POSITION_WITH_0_GAPS) == 10);
    assert(msa.getStatValByType(stat_type::MSA_POSITION_WITH_1_GAPS) == 2);
    assert(msa.getStatValByType(stat_type::MSA_POSITION_WITH_2_GAPS) == 12);
    assert(msa.getStatValByType(stat_type::MSA_POSITION_WITH_N_MINUS_1_GAPS) == 13);


    return 0;

}