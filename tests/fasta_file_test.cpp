#include <iostream>

#include "../src/MSA.h"

int main() {

    MSA msa("msas/dummy_msa.fasta");
    msa.printMSA();

    return 0;

}