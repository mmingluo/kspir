#include <cassert>

#include "secret.h"
#include "ntt.h"

#include "crt.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

// default construct: set Hanmming weight as 2/3*n
Secret::Secret(uint64_t module, bool ntt): modulus(module), nttform(ntt)
{
    this->length = N;
    this->lwe_type = RLWE;
#ifndef INTEL_HEXL
    int64_t* temp = new int64_t[this->length];
    for (size_t i = 0; i < length; i++)
    {
        // TODO: add modulus, add time random seed
        temp[i] = rand() % 3 - 1;
    }
    // save in ntt form
    ntt(temp);
    this->nttform = true;
    
    // save secret in ntt form
    for (size_t i = 0; i < length; i++)
    {
        // this->data.push_back(temp[i] < 0 ? (uint64_t)(temp[i] + mod) : (uint64_t)temp[i]);
        this->data.push_back(temp[i]);
    }
   
    delete[] temp;
#endif

#ifdef INTEL_HEXL
    // int64_t* temp = new int64_t[this->length];

    if (modulus == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, modulus);
        this->ntts = ntts;
    }

    for (size_t i = 0; i < length; i++)
    {
        int32_t temp = rand() % 3 - 1;
        this->data.push_back(temp < 0 ? (uint64_t)(temp + modulus) : (uint64_t)temp);
    }

    if (ntt)
    {
        this->ntts.ComputeForward(data.data(), data.data(), 1, 1);
    }
#endif
}

Secret::Secret(LWE_TYPE type, uint64_t module) : lwe_type(type), modulus(module)
{
    assert(type == LWE); // this type is lwe construct
    this->length = N;
    this->nttform = false;

#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }
#endif

    // int32_t* temp = new int32_t[this->length];
    for (size_t i = 0; i < length; i++)
    {
        // TODO: add modulus, add time random seed
        int32_t temp = rand() % 3 - 1;
        this->data.push_back(temp < 0 ? (uint64_t)(temp + modulus) : (uint64_t)temp);
    }
}

uint64_t Secret::getModulus()
{
    return this->modulus;
}

int32_t Secret::getLength()
{
    return this->length;
}

std::vector<uint64_t> Secret::getData()
{
    return this->data;
}

uint64_t Secret::getData(int32_t index)
{
    return this->data[index];
}

LWE_TYPE Secret::getLweType()
{
    return this->lwe_type;
}

bool Secret::isNttForm()
{
    return this->nttform;
}

#ifdef INTEL_HEXL
intel::hexl::NTT Secret::getNTT()
{
    return ntts;
}
#endif

void Secret::toCoeffForm()
{
     // TODO
#ifdef INTEL_HEXL
    ntts.ComputeInverse(data.data(), data.data(), 1, 1);
#endif
    this->nttform = false;  
}

void Secret::toNttForm()
{
    // TODO
#ifdef INTEL_HEXL
    ntts.ComputeForward(data.data(), data.data(), 1, 1);
#endif
    this->nttform = true;
}

Secret::~Secret()
{
}
