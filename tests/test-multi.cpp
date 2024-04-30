
#include <stdint.h>
#include <stddef.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <chrono>
#include <assert.h>

#include "pir.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

inline uint64_t barrett_raw_u64(uint64_t input, uint64_t const_ratio_1, uint64_t modulus)
{
    unsigned long long tmp[2];
    tmp[1] = static_cast<unsigned long long>(((static_cast<__uint128_t>(input) * static_cast<__uint128_t>(const_ratio_1)) >> 64));

    // Barrett subtraction
    tmp[0] = input - tmp[1] * modulus;

    // One more subtraction is enough
    return (tmp[0] >= modulus ? (tmp[0] - modulus) : (tmp[0]));
}

inline uint64_t barrett_coeff(uint64_t val, size_t n)
{
    // return val % moduli[n];
    return (n == 0) ? barrett_raw_u64(val, cr1_p, crtq1) : barrett_raw_u64(val, cr1_b, crtq2);
    return val;
}

void reorientCipher(uint64_t* cipherbuf, std::vector<std::vector<uint64_t> >& inputAcrt,
                    std::vector<std::vector<uint64_t> >& inputBcrt, int32_t N1, int32_t ringdim)
{
    assert(ringdim == inputAcrt[0].size() / 2);
    assert(ringdim == inputBcrt[0].size() / 2);

    for (size_t i = 0; i < ringdim; i++)
    {
        for (size_t j = 0; j < N1; j++)
        {
            // uint64_t temp = inputA[j][i];
            // cipherbuf[i * 2 * N1 + 2 * j] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
            // temp = inputB[j][i];
            // cipherbuf[i * 2 * N1 + 2 * j + 1] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
            cipherbuf[i * 2 * N1 + 2 * j] = inputAcrt[j][i] | inputAcrt[j][i + ringdim] << 32;
            cipherbuf[i * 2 * N1 + 2 * j + 1] = inputBcrt[j][i] | inputBcrt[j][i + ringdim] << 32;
        }
    }
}

// encode to crt form
void database_tocrt(uint64_t* datacrt, std::vector<std::vector<uint64_t> >& data_ntt,
                    const int32_t firstdim, const int32_t seconddim, const int32_t ringdim)
{    
    /** datacrt = [a, e, i, m, b, f, j, \cdots] 
    data_ntt = [a, b, c, d
                e, f, g, h
                i, j, k, l
                m, n, o, p]
    **/

    for (size_t i = 0; i < ringdim; i++)
    {
        for (size_t j = 0; j < seconddim; j++)
        {
            for (size_t k = 0; k < firstdim; k++)
            {
                uint64_t temp = data_ntt[j * firstdim + k][i];
                uint64_t temp2 = data_ntt[j * firstdim + k][i + ringdim];
                // TODO: check
                // datacrt[i * firstdim * seconddim + j * firstdim + k] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
                datacrt[i * firstdim * seconddim + j * firstdim + k] = temp | temp2 << 32; 
            }
        }
    }
}

void test_mul256MB()
{
    int32_t firstdim = 0x01 << 6; // 8;
    int32_t seconddim = 0x01 << 5; // 7;
    int32_t packing_num = 1; // 4;
    uint64_t ringdim = 4096;

    assert(ringdim == N);
    assert(ringdim == poly_len);

    uint64_t modulus = bigMod;

    // total_vec * ringdim * 8 bits = 256 MB
    int32_t total_vec = firstdim * seconddim * packing_num;

    std::cout << "total_vec = " << total_vec << std::endl;

    std::vector<std::vector<uint64_t> > data(total_vec, std::vector<uint64_t>(ringdim, 0));
    std::vector<std::vector<uint64_t> > data_ntt(total_vec, std::vector<uint64_t>(ringdim, 0));
    // sample database
    for (size_t i = 0; i < total_vec; i++)
    {
        for (size_t j = 0; j < ringdim; j++)
        {
            data[i][j] = rand() & 0x0ff;
        }
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(ringdim, modulus);
    for (size_t i = 0; i < total_vec; i++)
    {
        ntts.ComputeForward(data_ntt[i].data(), data[i].data(), 1, 1); // convert to ntt form
    }

    // dummy ciphertext
    std::vector<std::vector<uint64_t> > inputA(firstdim, std::vector<uint64_t>(ringdim, 0));
    std::vector<std::vector<uint64_t> > inputB(firstdim, std::vector<uint64_t>(ringdim, 0));
    for (size_t i = 0; i < firstdim; i++)
    {
        sample_random8_vector(inputA[i].data(), ringdim);
        ntts.ComputeForward(inputA[i].data(), inputA[i].data(), 1, 1);
        sample_random8_vector(inputB[i].data(), ringdim);
        ntts.ComputeForward(inputB[i].data(), inputB[i].data(), 1, 1);
    }

    // perform multiplication
    std::vector<std::vector<uint64_t> > resultA(seconddim, std::vector<uint64_t>(ringdim, 0));
    std::vector<std::vector<uint64_t> > resultB(seconddim, std::vector<uint64_t>(ringdim, 0));
    std::vector<uint64_t> tempA(N);
    std::vector<uint64_t> tempB(N);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();
for (size_t nt = 0; nt < ntimes; nt++)
{
    for (size_t i = 0; i < seconddim; i++)
    {  
        for (size_t j = 0; j < firstdim; j++)
        {
            intel::hexl::EltwiseMultMod(tempA.data(), inputA[j].data(), data_ntt[i * firstdim + j].data(),
                        ringdim, modulus, 1);
            intel::hexl::EltwiseMultMod(tempB.data(), inputB[j].data(), data_ntt[i * firstdim + j].data(),
                        ringdim, modulus, 1);
            intel::hexl::EltwiseAddMod(resultA[i].data(), resultA[i].data(), tempA.data(), ringdim,
                        modulus);
            intel::hexl::EltwiseAddMod(resultB[i].data(), resultB[i].data(), tempB.data(), ringdim,
                        modulus);
        }
    }
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "dummy multiplication costs " << glapsed.count() << " us." << std::endl;

    for (size_t i = 0; i < seconddim; i++)
    {
        ntts.ComputeInverse(resultA[i].data(), resultA[i].data(), 1, 1);
        ntts.ComputeInverse(resultB[i].data(), resultB[i].data(), 1, 1);
    }
    for (size_t i = 0; i < firstdim; i++)
    {
        ntts.ComputeInverse(inputA[i].data(), inputA[i].data(), 1, 1);
        ntts.ComputeInverse(inputB[i].data(), inputB[i].data(), 1, 1);
    }
#endif

    // crt version
    std::vector<std::vector<uint64_t> > inputAcrt(firstdim, std::vector<uint64_t>(2 * ringdim, 0));
    std::vector<std::vector<uint64_t> > inputBcrt(firstdim, std::vector<uint64_t>(2 * ringdim, 0));
    for (size_t i = 0; i < firstdim; i++)
    {
        computeForward(inputAcrt[i].data(), inputA[i].data());
        computeForward(inputBcrt[i].data(), inputB[i].data());
    }
    std::vector<std::vector<uint64_t> > data_nttcrt(total_vec, std::vector<uint64_t>(2 * ringdim, 0));
    for (size_t i = 0; i < total_vec; i++)
    {
        computeForward(data_nttcrt[i].data(), data[i].data());
    }

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * firstdim * 2 * N);
    reorientCipher(cipherbuf, inputAcrt, inputBcrt, firstdim, ringdim);
    size_t num_words = firstdim * seconddim * N;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    database_tocrt(datacrt, data_nttcrt, firstdim, seconddim, ringdim);

    std::vector<std::vector<uint64_t> > output(seconddim, std::vector<uint64_t>(4 * N, 0));

    auto start_crt = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    fastMultiplyQueryByDatabaseDim1(output, datacrt, cipherbuf, firstdim, seconddim);
}
    auto stop_crt = std::chrono::high_resolution_clock::now();

    auto glapsed_crt = std::chrono::duration_cast<std::chrono::microseconds>(stop_crt - start_crt);
    std::cout << "crt dummy multiplication costs " << glapsed_crt.count() << " us." << std::endl;

    // check
    std::vector<uint64_t> result2A(N, 0);
    std::vector<uint64_t> result2B(N, 0);
    computeInverse(result2A.data(), output[2].data());
    computeInverse(result2B.data(), output[2].data() + 2 * N);

    std::cout << "checked result" << std::endl;

    showLargeVector(resultA[2], "resultA[2] = ");
    showLargeVector(result2A, "result2A = ");
    std::cout << std::endl;
    showLargeVector(resultB[2], "resultB[2] = ");
    showLargeVector(result2B, "result2B = ");
}

int main(int argc, char** argv)
{
    // srand(time(NULL));

    test_mul256MB();

    return 0;
}
