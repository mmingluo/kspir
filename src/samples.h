#ifndef SAMPLES_H
#define SAMPLES_H

#include <stdint.h>
#include <stddef.h>
#include <vector>

void sample_random(uint64_t* a, uint64_t modulus, int32_t size);

void sample_random(std::vector<uint64_t>& result, uint64_t modulus);

uint64_t sample_guass(uint64_t modulus);

/**
 * @brief 
 * 
 */
void sample_guass(uint64_t* result, uint64_t modulus);

/**
 * @brief evaluate automorphic transform
 * 
 * @param result 
 * @param modulus1 input modulus
 * @param modulus2 output modulus
 */
void guass_to_modulus(uint64_t* result, uint64_t modulus1, uint64_t modulus2);

void sample_database(uint64_t** data);

void sample_database(std::vector<std::vector<uint64_t> >& data);

// for dummy data test
uint64_t sample_random62();

void sample_random62_vector(uint64_t* a, size_t length);

uint64_t sample_random8();

void sample_random8_vector(uint64_t* a, size_t length);

#endif
