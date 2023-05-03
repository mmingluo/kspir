#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include <stddef.h>
#include <vector>
#include <string>

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

void copy_vector(int64_t* result, uint64_t* a);

void copy_vector(uint64_t* result, int64_t* a, uint64_t modulus);

template<typename T1, typename T2>
void copy_vector(T1* result, T2* a);

template<typename T>
void showLargeVector(std::vector<T>& vals);

void showLargeVector(std::vector<uint64_t>& vals, std::string ss = "");

void showLargeIntervalVector(std::vector<uint64_t>& vals, int32_t interval, std::string ss = "");

template<typename T>
void showVector(std::vector<T>& vals, std::string ss = "");

void showVector(std::vector<uint64_t>& vals, std::string ss = "");

uint64_t QInv(int lwenum);

void database_to_signed(std::vector<std::vector<uint64_t> >& data, uint64_t pbits, uint64_t modulus);

void database_to_rnssigned(std::vector<std::vector<uint64_t> >& data, uint64_t pbits,
                            uint64_t modulus1, uint64_t modulus2);

void database_tontt(std::vector<std::vector<uint64_t> >& data);

void database_to_rnsntt(std::vector<std::vector<uint64_t> >& data);

void data_to_setupdata(std::vector<std::vector<uint64_t> >& setup_data,
                        std::vector<std::vector<uint64_t> >& data,
                        uint64_t pbits, uint64_t modulus1, uint64_t modulus2);

void negate(uint64_t* result, uint64_t length, uint64_t modulus);

void multConst(std::vector<uint64_t>& result, uint64_t cosntNum);

void lweToRlwe(std::vector<uint64_t>& result);

void transpose(std::vector<std::vector<uint64_t> >& a);

std::vector<uint64_t> powerOfBg(uint64_t base, uint64_t BBg, int32_t ellnum);

std::vector<uint128_t> powerOfBg128(uint64_t base, uint64_t BBg, int32_t ellnum);

// inline uint64_t powerOfBg(uint64_t base, uint64_t BBg, int32_t num);

void crt_inv(std::vector<uint128_t>& result, const std::vector<uint64_t>& input1,
            const std::vector<uint64_t>& input2, uint64_t modulus1, uint64_t modulus2);

void decompose(uint64_t** result, const uint64_t* input,
                int32_t ellnum, uint64_t base, uint64_t BBg);

void decompose(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
                int32_t ellnum, uint64_t base, uint64_t BBg);

std::vector<uint64_t> recontruct(std::vector<std::vector<uint64_t> >& dec_a,
                            int32_t ellnum, uint64_t base, uint64_t BBg);

void check_recontruct(std::vector<std::vector<uint64_t> >& dec_a, const std::vector<uint64_t>& a,
                int32_t ellnum, uint64_t base, uint64_t BBg);

void decompose_variant(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
                int32_t ellnum, uint64_t base, uint64_t BBg);

void decompose_rlwe(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input1,
                const std::vector<uint64_t>& input2, int32_t ellnum, uint64_t base, uint64_t BBg);

void decompose_crt(std::vector<std::vector<uint64_t> >& result1, std::vector<std::vector<uint64_t> >& result2,
                const std::vector<uint64_t>& input1, const std::vector<uint64_t>& input2,
                int32_t ellnum, uint64_t base, uint64_t BBg);

inline void element_to_vector(std::vector<uint64_t>& result, uint64_t input)
{
    for (auto iter = result.begin(); iter != result.end(); iter++)
        *iter = input;
}

void automorphic(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,
                const int32_t index, const uint64_t modulus);

void encode_crt(std::vector<uint64_t>& result);

void compute_indicator(int32_t& data_index, bool& reverse, size_t j, int32_t s_index);

#endif
