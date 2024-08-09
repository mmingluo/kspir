#ifndef ANSWER_H
#define ANSWER_H

#include <stdint.h>
#include <stddef.h>
#include <vector>

#include "lwe.h"
#include "secret.h"

void decompose_a(uint64_t** dec_a, uint64_t* a);

#ifndef INTEL_HEXL
void ksKey_hint(uint64_t* b_rlwe, uint64_t* a_rlwe, uint64_t* a, uint64_t** ks);
#else
int ksKey_hint(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& ks);

int ksKey_hint_variant(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& ks);

int ksKey_hint_first_construction(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& pk,
                std::vector<std::vector<uint64_t> >& data,
                int32_t ellnum = 3, uint64_t base = 0x01 << 17, uint64_t BBg = 0x01 << 11);
#endif

#ifndef INTEL_HEXL
int answer(uint64_t* b_rlwe, uint64_t* a_rlwe, uint64_t* b, uint64_t* a,
           std::vector<std::vector<uint64_t> >& ks, std::vector<std::vector<uint64_t> >& data,
           bool isntt_b = 1, bool isntt_data = 1);
#else
int answer_first_dimension(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
            std::vector<uint64_t>& query_b,// std::vector<uint64_t>& a,
            std::vector<std::vector<uint64_t> >& data,// std::vector<std::vector<uint64_t> >& ks,
            bool isntt_b = 1, bool isntt_data = 1);

/**
 * @brief 
 * 
 * @param result 
 * @param b_rlwe in coefficient form
 * @param a_rlwe in coefficient form
 * @param b in ntt form
 * @param data in ntt form
 * @param query 
 * @param isntt_b 
 * @param isntt_data 
 * @return int 
 */
int answer_two_dimensions(RlweCiphertext& result,
            std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
            std::vector<uint64_t>& b,
            std::vector<std::vector<uint64_t> >& data,
            RGSWCiphertext& query,
            bool isntt_b = 1, bool isntt_data = 1);

void pkKeyGen(std::vector<std::vector<uint64_t> >& pk,
                 Secret& queryKey, Secret& answerKey,
                 int32_t ellnum = 3, uint64_t base = 0x01 << 17, uint64_t BBg = 0x01 << 11);

void dummy_ksKeyGen(std::vector<std::vector<uint64_t> >& ks,
                 Secret& queryKey, Secret& answerKey,
                 std::vector<std::vector<uint64_t> >& data);
#endif

/**
 * @brief evaluate automorphic transform
 * 
 * @param result 
 * @param index 
 * @param autokey 
 */
void evalAuto(RlweCiphertext& result, const int32_t index, const AutoKey& autokey);

void packingLWEs(RlweCiphertext& result, std::vector<RlweCiphertext>& lwes,
                 const AutoKey& autokey);

void packingRLWEs(RlweCiphertext& result, std::vector<RlweCiphertext>& rlwes,
                 const AutoKey& autokey);

void evalAutoRNS(RlweCiphertext& result1, RlweCiphertext& result2,
                 const int32_t index, const AutoKeyRNS& autokey);

void evalExpandRNS(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
                    const AutoKeyRNS& autokey);

void evalTrNnRNS(std::vector<RlweCiphertext>& result, int32_t log2Nn, const AutoKeyRNS& autokey);

void publicKeyGen(std::vector<std::vector<uint64_t> >& publicKey,
                Secret& queryKey, Secret& answerKey);

void evalKsKeyGen(std::vector<std::vector<uint64_t> >& ks,
                 std::vector<std::vector<uint64_t> >& pk,
                 std::vector<std::vector<uint64_t> >& data);

/**
 * @brief evaluate automorphic transform
 * 
 * @param result dim1: 2 * ell * N, dim2: N. store in ntt form
 * @param ks dim1: 2 * ell * N * 2(CRT), dim2: N. store in ntt form
 * @param 
 */
void modSwitch(std::vector<std::vector<uint64_t> >& result,
                std::vector<std::vector<uint64_t> >& ks);

/**
 * @brief rns modswitch
 * 
 * @param result store in ntt form
 * @param input1 rlwe ciphertext in modulus1
 * @param input2 rlwe ciphertext in modulus2
 */
void modSwitch(RlweCiphertext& result, RlweCiphertext& input1, RlweCiphertext& input2);

#endif
