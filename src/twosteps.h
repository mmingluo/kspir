#ifndef TWO_STEPS_H
#define TWO_STEPS_H

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "lwe.h"
#include "secret.h"

void sample_database_bsgs(std::vector<std::vector<uint64_t> >& data);

/**
 * @brief database to bsgs ntt form
 * @param data input database
 * @param modulus 
 *  
*/
void database_tobsgsntt(std::vector<std::vector<uint64_t> >& result,
                        std::vector<std::vector<uint64_t> >& data, uint64_t modulus, int32_t N1 = N / 2);

void inverse_encode(std::vector<uint64_t>& message);

void query_bsgs(RlweCiphertext& cipher, Secret& queryKey, uint64_t row);

void query_bsgs_rns(std::vector<RlweCiphertext>& cipher, Secret& queryKey, uint64_t row);

void decrypt_bsgs(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret);

void decrypt_bsgs_total(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret,
                        const int32_t rlwesNum);

void evalAutoRNSCRT(RlweCiphertext& result, std::vector<std::vector<uint64_t> >& input,
              const int32_t N2, const AutoKeyBSGSRNS& autokey,
              const int32_t ellnum, const uint64_t PP, const uint64_t BBg);

void reorientCipher(uint64_t* cipherbuf, RlweCiphertext& input,
                    std::vector<RlweCiphertext>& rotated_cipher, int32_t N1);

void matrix_vector_mul(RlweCiphertext& result, RlweCiphertext& input,
                        std::vector<std::vector<uint64_t> >& data,
                        const AutoKeyBSGS& autoKey);

void matrix_vector_mul_bsgs(RlweCiphertext& result, RlweCiphertext& input,
                        std::vector<std::vector<uint64_t> >& data,
                        const AutoKeyBSGS& autoKey, const int32_t N1);

void matrix_vector_mul_bsgs_crt(RlweCiphertext& result, RlweCiphertext& input,
                        uint64_t* data, const AutoKeyBSGS& autoKey, const int32_t N1);

void matrix_vector_mul_bsgs_rns_crt(RlweCiphertext& result, std::vector<RlweCiphertext>& input,
                        uint64_t* data, const AutoKeyBSGSRNS& autoKey, const int32_t N1,
                        std::vector<std::vector<int32_t> >& permutations);

void matrix_vector_mul_bsgs_rns_crt_large(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
                        uint64_t* data, const AutoKeyBSGSRNS& autoKey, const int32_t N1,
                        std::vector<std::vector<int32_t> >& permutations, int32_t r = 1);

void genAutoKeyFromOffline(AutoKeyBSGSRNS& result, AutoKeyBSGSRNS& result1, AutoKeyBSGSRNS& result2,
              const AutoKeyBSGSRNS& offlineKey1, const AutoKeyBSGSRNS& offlineKey2,
              const int32_t N1);

#endif
