#include <stddef.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <stdint.h>

#include "pir.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

void test_ntt_crt()
{
    size_t length = 4096;
    uint64_t modulus = crtMod;

    std::vector<uint64_t> input1(length, 0);

    // sample_random8_vector(input1.data(), length);
    sample_random(input1, modulus);

#ifdef INTEL_HEXL
    // uint64_t root_of_unity = 10297991595; //minimum root of unity
    intel::hexl::NTT ntts(length, modulus, root_of_unity_crt);

int ntimes = 1000;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    ntts.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts.ComputeInverse(input1.data(), input1.data(), 1, 1);
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt costs " << glapsed.count() << " us." << std::endl;
#endif
}

void test_ntt_crt_miniroot()
{
    size_t length = 4096;
    uint64_t modulus = crtMod;

    std::vector<uint64_t> input1(length, 0);

    // sample_random8_vector(input1.data(), length);
    sample_random(input1, modulus);

#ifdef INTEL_HEXL
    uint64_t root_of_unity = 10297991595; //minimum root of unity
    intel::hexl::NTT ntts(length, modulus, root_of_unity);

int ntimes = 1000;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    ntts.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts.ComputeInverse(input1.data(), input1.data(), 1, 1);
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (mini root) costs " << glapsed.count() << " us." << std::endl;
#endif
}

void test_ntt50()
{
    size_t length = 4096;
    uint64_t modulus = bigMod;

    std::vector<std::vector<uint64_t> > input1(length, std::vector<uint64_t>(length, 0));

    // sample_random8_vector(input1.data(), length);
    for (size_t i = 0; i < length; i++)
        sample_random(input1[i], modulus);

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);

int ntimes = 10;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t nt = 0; nt < ntimes; nt++)
{
    for (size_t i = 0; i < length; i++)
    {
        ntts.ComputeForward(input1[i].data(), input1[i].data(), 1, 1);
        ntts.ComputeInverse(input1[i].data(), input1[i].data(), 1, 1);
    }
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (56 bits) costs " << glapsed.count() << " us." << std::endl;
#endif
}

void test_ntt56()
{
    size_t length = 4096;
    uint64_t modulus = intel::hexl::GeneratePrimes(1, 55, false, 4096)[0];
    modulus = 72057594037641217;

    std::vector<std::vector<uint64_t> > input1(length, std::vector<uint64_t>(length, 0));

    // sample_random8_vector(input1.data(), length);
    for (size_t i = 0; i < length; i++)
        sample_random(input1[i], modulus);

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);

int ntimes = 10;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t nt = 0; nt < ntimes; nt++)
{
    for (size_t i = 0; i < length; i++)
    {
        ntts.ComputeForward(input1[i].data(), input1[i].data(), 1, 1);
        ntts.ComputeInverse(input1[i].data(), input1[i].data(), 1, 1);
    }
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (56 bits) costs " << glapsed.count() << " us." << std::endl;
#endif
}


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    // test_ntt_crt();
    // test_ntt_crt_miniroot();
    test_ntt50();
    // test_ntt56();

    return 0;
}
