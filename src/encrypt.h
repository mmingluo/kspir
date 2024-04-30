#include <stdint.h>

#include "secret.h"
#include "lwe.h"

/**
 * @brief encrypt an LWE cipertext
 * 
 * @param b 
 * @param a 
 * @param secret 
 * @param message 
 */
void encrypt(uint64_t* b, uint64_t* a, Secret& secret, uint64_t message);

/**
 * @brief encrypt an RLWE ciphertex
 * 
 * @param b return in ntt form
 * @param a return in ntt form
 * @param secret coefficient or ntt form
 * @param message note: shoule be ntt form
 */
void encrypt(uint64_t* b, uint64_t* a, Secret& secret, uint64_t* message);

/**
 * @brief encrypt an RLWE ciphertex
 * 
 * @param b1 return in ntt form
 * @param a1 return in ntt form
 * @param b1 another part of CRT, return in ntt form
 * @param a1 another part of CRT, return in ntt form
 * @param secret coefficient or ntt form
 * @param secret1 store in ntt form
 * @param message note: shoule be ntt form
*/
void encrypt(uint64_t* b1, uint64_t* a1, uint64_t* b2, uint64_t* a2,
            Secret& secret, uint64_t* message,
            uint64_t* secret1, uint64_t* secret2,
            intel::hexl::NTT& ntts1, intel::hexl::NTT& ntts2,
            uint64_t modulus1, uint64_t modulus2);


/**
 * @brief note that the message will be crypted as: message * bigMod
 * 
 * @param cipher1 
 * @param cipher2 
 * @param secret 
 * @param message 
 */
void encrypt_rns(RlweCiphertext& cipher1, RlweCiphertext& cipher2,
                Secret& secret, std::vector<uint64_t>& message);

void encrypt_rns(std::vector<uint64_t>& b1, std::vector<uint64_t>& a1,
                    std::vector<uint64_t>& b2, std::vector<uint64_t>& a2,
                    Secret& secret, std::vector<uint64_t>& message1, std::vector<uint64_t>& message2,
                    int32_t num = 0, uint64_t modulus2 = bigMod2);

void encrypt_rns_bsgs_autokey(std::vector<uint64_t>& b1, std::vector<uint64_t>& a1,
                    std::vector<uint64_t>& b2, std::vector<uint64_t>& a2,
                    Secret& secret, std::vector<uint64_t>& message1, std::vector<uint64_t>& message2,
                    uint64_t modulus2);

/**
 * @brief used to remove term 2^l. (lazy variant)
 * Refer to `4.1 Removing the Leading Term N' in [].
 *
 * 
 * @param b 
 * @param a 
 * @param secret 
 * @param message 
 * @param lwesnum 
 */
void encrypt_special_rlwe(uint64_t* b, uint64_t* a, Secret& secret, uint64_t* message, int32_t lwesnum);

void encrypt(RlweCiphertext& cipher, Secret& secret, std::vector<uint64_t>& message);

/**
 * @brief encrypt for rns bsgs algorithm
 * @param cipher encrypted cipher
 * @param secret secret
 * @param message stored in coefficients form
 *
*/
void encrypt_rns_bsgs(std::vector<RlweCiphertext>& cipher, Secret& secret, std::vector<uint64_t>& message);

/**
 * @brief 
 * 
 * @param b 
 * @param a 
 * @param secret 
 * @return uint64_t 
 */
uint64_t decrypt(uint64_t* b, uint64_t* a, Secret& secret);

/**
 * @brief decrypt an RLWE ciphertext
 * 
 * @param message 
 * @param b 
 * @param a 
 * @param secret 
 */
void decrypt(uint64_t* message, uint64_t* b, uint64_t* a, Secret& secret);

void decrypt(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret);

/**
 * @brief decrypt an RLWE ciphertext in another modulus
 * 
 * @param message 
 * @param b 
 * @param a 
 * @param secret 
 * @param modulus2 the second modulus
 */
void decrypt(uint64_t* message, uint64_t* b, uint64_t* a, Secret& secret, uint64_t modulus2);
