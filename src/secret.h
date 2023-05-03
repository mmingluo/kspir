#ifndef SECRET_H
#define SECRET_H

#include <stdint.h>
#include <vector>

#include "params.h"


#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

enum LWE_TYPE {LWE, RLWE};

class Secret
{
private:
    /* data */
    LWE_TYPE lwe_type = RLWE;

    uint64_t length = N;
    uint64_t modulus;

    std::vector<uint64_t> data;

    bool nttform = false;

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
#endif

public:

    Secret(uint64_t module = bigMod);

    Secret(LWE_TYPE type, uint64_t module = mod); // LWE construct

    std::vector<uint64_t> getData();

    uint64_t getData(int32_t index);

    uint64_t getModulus();

    int32_t getLength();

    LWE_TYPE getLweType();

    bool isNttForm();

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT();
#endif

    void toCoeffForm();

    void toNttForm();

    ~Secret();
};

#endif
