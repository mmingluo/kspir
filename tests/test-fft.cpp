
#include <stdint.h>
#include <stddef.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <chrono>
#include <assert.h>

#include "pir.h"


void test_fft_multiply()
{
    int32_t length = N;
    uint64_t modulus = 0x01UL << 56;
    uint64_t pmod = 0x01 << 16;

    std::vector<uint64_t> input1(length, 0);
    sample_random(input1, modulus);

    std::vector<uint64_t> message(length, 0);
    sample_random(message, pmod);

    std::vector<double> input1_double(length);
    for (size_t i = 0; i < length; i++)
    {
        input1_double[i] = input1[i];
    }

    std::vector<long double> message_double(length);
    for (size_t i = 0; i < length; i++)
    {
        message_double[i] = message[i];
    }

    std::vector<long double> output_double(length);
    for (size_t i = 0; i < length; i++)
    {
        output_double[i] = input1_double[i] * message_double[i];
    }

    std::vector<uint64_t> output(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        output[i] =  static_cast<uint128_t>(output_double[i]) & ((0x01UL << 56) - 1);
    }

    std::cout << "finished" << std::endl;
}

int main(int argc, char** argv)
{
    // srand(time(NULL));

    test_fft_multiply();

    return 0;
}
