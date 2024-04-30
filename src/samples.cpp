#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <iostream>

#include "params.h"
#include "samples.h"

void sample_random(uint64_t* a, uint64_t modulus, int32_t size)
{
    
    for (size_t i = 0; i < size; i++)
    {
        // TODO: this sampling is not uniform
        uint64_t temp = ((uint64_t)rand() << 31) | rand();
        a[i] = temp % modulus;
    }
}

void sample_random(std::vector<uint64_t>& result, uint64_t modulus)
{
    /*
    for (auto re : result)
    {
        // TODO: this sampling is not uniform
        uint64_t temp = ((uint64_t)rand() << 31) | rand();
        re = temp % modulus;
    }
    */
    for (auto iter = result.begin(); iter != result.end(); iter++)
    {
        uint64_t temp = ((uint64_t)rand() << 31) | rand();
        *iter = temp % modulus;
    }

}

uint64_t sample_guass(uint64_t modulus)
{
    // TODO:    
    int32_t temp = rand() % 3 - 1;

    return temp < 0 ? modulus - 1: temp;
}

void sample_guass(uint64_t* result, uint64_t modulus)
{
    // TODO:
    for (size_t i = 0; i < N; i++)
    {
        int32_t temp = rand() %3 -1;
        result[i] = temp < 0 ? modulus - 1: temp;
    }
}

void guass_to_modulus(uint64_t* result, uint64_t modulus1, uint64_t modulus2)
{
    for (size_t i = 0; i < N; i++)
    {
        if (result[i] > modulus1 / 2)
        {
            result[i] = (int64_t)result[i] - modulus1 + modulus2;
        }
    }
}

void sample_database(uint64_t** data)
{
    // the database has N * N entries, and each entry has 8 bits
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            // TODO
            data[i][j] = rand() & 0x0ff;
        }
    }

}

void sample_database(std::vector<std::vector<uint64_t> >& data)
{
    // the database has N * N entries, and each entry has 16 bits    
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            data[i][j] = rand() & 0x0ffff;
        }
    }

}

uint64_t sample_random62()
{
    uint64_t temp = rand();
    return (temp << 31) | rand();
}

void sample_random62_vector(uint64_t* a, size_t length)
{
    for (size_t i = 0; i < length; i++)
    {
        uint64_t temp = rand();
        a[i] = (temp << 31) | rand();
    }
}

uint64_t sample_random8()
{
    return rand() & 0x0ff;
}

void sample_random8_vector(uint64_t* a, size_t length)
{
    for (size_t i = 0; i < length; i++)
    {
        a[i] = rand() & 0x0ff;
    }
}
