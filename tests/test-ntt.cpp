#include <stdint.h>
#include <stddef.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "ntt.h"
#include "pir.h"

// #ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
// #endif

int64_t barret_reduce_hexl(int128_t a) {
    int64_t t;
    
    // v = round(2^(40+64) /bigQ.);
    const int64_t v = 4611686044196143104;
    t = (__int128_t)v*a >> 104;

    t *= bigMod;
    return a - t;
}

uint64_t bit_inverse(uint64_t a, int loglen)
{
    uint64_t result = 0;

    for (size_t i = 0; i < loglen; i++)
    {
        uint64_t temp = (a >> i) & 0x01;
        result |= temp << (loglen - i - 1); 
    }

    return result;
}

void compute_bits_inverse(int loglen) // e.g. len = 1024
{
    int64_t len = 0x01 << loglen;
    std::vector<uint64_t> revs;

    for (size_t i = 0; i < len; i++)
    {
        revs.push_back(bit_inverse(i, loglen));
    }

    for (auto r : revs)
    {
        std::cout << r << ", ";
    }
}

/*
// refer to intel hexl [https://github.com/intel/hexl]
uint64_t InverseMod(uint64_t input, uint64_t modulus) {
  uint64_t a = input % modulus;
  // HEXL_CHECK(a != 0, input << " does not have a InverseMod");

  if (modulus == 1) {
    return 0;
  }

  int64_t m0 = static_cast<int64_t>(modulus);
  int64_t y = 0;
  int64_t x = 1;
  while (a > 1) {
    // q is quotient
    int64_t q = static_cast<int64_t>(a / modulus);

    int64_t t = static_cast<int64_t>(modulus);
    modulus = a % modulus;
    a = static_cast<uint64_t>(t);

    // Update y and x
    t = y;
    y = x - q * y;
    x = t;
  }

  // Make x positive
  if (x < 0) x += m0;

  return uint64_t(x);
}
*/

// test ntt when q = 8380417 (\approx 2^23), N = 2048
void test_small_ntt()
{
    int32_t a[dim] = {1, 6};
    int32_t s[dim] = {1, 1, 1};
    int32_t result[dim];

    ntt(a);
    ntt(s);
    hadamard_mult(result, a, s);
    invntt_tomont(result);
}

// test ntt when q = 4398046486529 (\approx 2^42), N = 2048
void test_big_ntt()
{
    int64_t a[dim] = {1, 6};
    int64_t s[dim] = {1, 1, 1};
    int64_t result[dim];
    // std::cout << "dim = " << dim << std::endl;
    ntt(a);
    ntt(s);
    hadamard_mult(result, a, s);
    invntt_tomont(result);
}

void test_barret_reduce()
{
    int64_t t = 4398046486531 + 4398046486531;
    std::cout << barret_reduce(t) << std::endl;
}

uint64_t powd(uint64_t x, uint64_t y, uint64_t modulus)
{
    uint64_t s = 1;
    for (size_t i = 0; i < y; i++)
    {
        s = (__int128_t)s * x % modulus;
    }

    return s;
}

void test_omega()
{
    uint64_t omega = 384399401;
    uint64_t modulus = 4398046486529;

    for (size_t i = 1; i < 2 * dim; i += 2)
    {
        std::cout << powd(omega, i, modulus) << ", ";
    }

    std::cout << std::endl;
    std::cout << "dim = " << dim << std::endl;
    std::cout << "omega^(2*dim) = " << powd(omega, 2 * dim, modulus) << std::endl;
}

inline int16_t barret_reduce_7681(int16_t a)
{
    int16_t t;
    // v = round(2^(11+16) /7681.);
    const int16_t v = 17474;
    t = (int32_t) v*a >> 27;

    t *= 7681;
    return a - t;
}

inline int16_t barret_reduce_7681_mult(int32_t a)
{
    int16_t t;
    // v = round(2^(11+32) /7681.);
    const int32_t v = 1145175501;
    t = (int64_t) v*a >> 43;

    t *= 7681;
    return a - t;
}

/*
// return corresponding [0, 1, 2,……, p-1]
int16_t barrett_reduce(int16_t a)
{
    int16_t t;
#if SCHEME_P == 3329
    // v = round(2^(10+16) /3329.);
    const int16_t v = 20159;
    t = (int32_t)v*a >> 26;
#elif SCHEME_P == 769
    // v = round(2^(8+16) /769.);
    const int16_t v = 21817;
    t = (int32_t)v*a >> 24;
#else
#error "Invalid SCHEME_P! SCHEME_P should be belong to {3329, 769}."
#endif
    t *= SCHEME_P;
    return a - t;
}
*/

// the returned value is not strictly belong to [0, q-1] 
inline uint16_t special_reduce_7681(int16_t a)
{
    uint16_t t = a >> 13;
    uint16_t u = a & 0x1fff;
    u -= t;
    t <<= 9;
    return u + t;
}

void test_reduce_7681()
{
    uint64_t length = 0x01 << 25;
    uint64_t modulus = 7681;
    
    std::vector<uint16_t> input1(length);
    std::vector<uint16_t> input2(length);
    std::vector<uint16_t> output(length);
    // sample_random(input1, modulus);
    // sample_random(input2, modulus);
    for (size_t i = 0; i < length; i++)
    {
        input1[i] = rand() % modulus;
        input2[i] = rand() % modulus;
    }

    /*****************************  specail reduce  **********************************/
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < length; i++)
    {
        output[i] = special_reduce_7681(input1[i] * input2[i]);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "special time costs " << glapsed.count() << " us" << std::endl;

    /*****************************  barret reduce  **********************************/
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < length; i++)
    {
        output[i] = barret_reduce_7681_mult((int32_t)input1[i] * input2[i]);
    }
    stop = std::chrono::high_resolution_clock::now();
    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "barret reduce time costs " << glapsed.count() << " us" << std::endl;

    /*****************************  native reduce  **********************************/
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < length; i++)
    {
        output[i] = ((int32_t)input1[i] * input2[i]) % 7681;
    }
    stop = std::chrono::high_resolution_clock::now();
    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "native reduce time costs " << glapsed.count() << " us" << std::endl;

    std::cout << "";
}

void test_reduce()
{
    std::vector<uint64_t> input1(N);
    std::vector<uint64_t> input2(N);
    std::vector<uint64_t> output(N);
    sample_random(input1, bigMod);
    sample_random(input2, bigMod);

    output[0] = (uint128_t)input1[0] * input2[0] % bigMod;
    output[0] = static_cast<uint128_t>(input1[0]) * input2[0] % bigMod;
    output[0] = barret_reduce((uint128_t)input1[0] * input2[0]); // note: error reduce
}

/*
void BM_EltwiseReduceModInPlace()
{  //  NOLINT
  size_t input_size = 1024;
  uint64_t modulus = 0xffffffffffc0001ULL;

  //auto input1 =
  //    intel::hexl::GenerateInsecureUniformIntRandomValues(input_size, 0, 100 * modulus);
  // const uint64_t input_mod_factor = modulus;
  // const uint64_t output_mod_factor = 1;
    
    // intel::hexl::EltwiseReduceMod(input1.data(), input1.data(), input_size, modulus,
    //                input_mod_factor, output_mod_factor);
}
*/

// #ifdef INTEL_HEXL
void test_intel_hexl()
{
    // uint64_t length = 256;
    // uint64_t modulus = 7681;
    uint64_t length = N;
    uint64_t modulus = bigMod;
    intel::hexl::NTT ntts(length, modulus);

    std::vector<uint64_t> a(length, 0), b(length, 0);
    
    for (size_t i = 0; i < length; i++)
        a[i] = rand() % modulus;

    std::cout << "a[0] = " << a[0] << std::endl;

    ntts.ComputeForward(a.data(), a.data(), 1, 1);
    ntts.ComputeInverse(b.data(), a.data(), 1, 1);

    uint64_t temp = 0;

    for (auto iter = a.begin(); iter != a.end(); iter++)
        temp += *iter;

    // temp = (uint128_t)temp * bNinv % modulus;
    uint128_t t = (uint128_t)temp * bNinv;
    barret_reduce_hexl(t);

    intel::hexl::InverseMod(length, modulus);

    // intel::hexl::BarrettReduce128()
    std::cout << "sum(ntt(a)) * N^-1 = " << temp << std::endl;
}
// #endif

int main(int argc, char** argv)
{
    srand(time(NULL));

    // test_small_ntt();
    // test_big_ntt();
    // test_barret_reduce();
    // compute_bits_inverse(11);
    // test_omega();

#ifdef INTEL_HEXL
    /**
     * @brief intel hexl test
     * 
     */
    // BM_EltwiseReduceModInPlace();
    // test_intel_hexl();
#endif

    // test_reduce_7681();
    test_reduce();

    return 0;
}

