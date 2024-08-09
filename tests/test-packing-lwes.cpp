#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"

void test_packing_lwes()
{
    int32_t l = 6; // 6; // 5
    int32_t lwesNum = 0x01 << l;
    // uint64_t modulus = bigMod;
    uint64_t modulus = crtMod;

    // Secret answerKey(bigMod);
    Secret answerKey(LWE, modulus);

    uint64_t length = answerKey.getLength();

    std::vector<RlweCiphertext> lwes;
    for (size_t i = 0; i < lwesNum; i++)
    {
        lwes.push_back(RlweCiphertext(N, modulus));
    }

    std::vector<uint64_t> message(lwesNum), decryptd_message(length);
    sample_random(message, modulus);
    //for (auto iter = message.begin(); iter != message.end(); iter++)
    //    *iter = 1234567891234;
    // showVector(message, "The message is ");
    showLargeVector(message, "The message is ");

    auto start_b = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < lwesNum; i++)
    {
        encrypt(lwes[i].b.data(), lwes[i].a.data(), answerKey, message[i]);
        lweToRlwe(lwes[i].a, modulus);
    }
    auto stop_b = std::chrono::high_resolution_clock::now();

    auto glapsed_b = std::chrono::duration_cast<std::chrono::microseconds>(stop_b - start_b);
    std::cout << "encryption costs " << glapsed_b.count() << " us." << std::endl;

    AutoKey autokey(length, modulus, 4, 0, 0x01 << 14);
    // AutoKey autokey(length, modulus, 3, 0x01 << 2, 0x01 << 18); section 4.2 in eprint 2020-015
    // AutoKey autokey;
    autokey.keyGen(answerKey, lwesNum);

    // the results
    RlweCiphertext result(N, modulus);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    packingLWEs(result, lwes, autokey);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " packing lwes costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message, result, answerKey);

    multConst(decryptd_message, QInv(lwesNum, modulus), modulus); // remove term 2^l
    std::cout << std::endl;

    showLargeIntervalVector(decryptd_message, length / lwesNum, "The decrypted value is ");
}

/**
 * @brief packing/convert r RLWEs into a single RLWE
 * algorithms include external product + the subprocedure of r-LWE-to-RLWE algorithm
 * in [https://eprint.iacr.org/2020/015].
 * 
 * input a target_row \in [r]
 * output data[target_row], where the size of data[target_row] is around 8 KB.
 */
void test_packing_rlwes()
{    
    int32_t l = 4; // 4; // 6; // 5
    int32_t rlwesNum = 0x01 << l; // the number of RLWEs

    int32_t target_row = rand() % rlwesNum;
    std::cout << "rlweNum: " << rlwesNum << std::endl;
    std::cout << "target row: " << target_row << std::endl;

    uint64_t modulus = crtMod;

    Secret answerKey(crtMod, false);

    uint64_t length = answerKey.getLength();

    std::vector<RlweCiphertext> rlwes;
    std::vector<RlweCiphertext> ext_rlwes;
    for (size_t i = 0; i < rlwesNum; i++)
    {
        rlwes.push_back(RlweCiphertext(N, modulus));
        ext_rlwes.push_back(RlweCiphertext(N, modulus));
    }

    std::vector<std::vector<uint64_t> > message(rlwesNum, std::vector<uint64_t>(N));
    std::vector<uint64_t> decryptd_message(length);
    
    for (size_t i = 0; i < rlwesNum; i++)
       sample_random(message[i], modulus);

    // add dummy messages
    /** * for target_row == 0
     * the input RLWE ciphertexts are:
     * [1, ..., r + 1, ...]
     * [2, ..., r + 2, ...]
     * [3, ..., r + 3, ...]
     * ...
     * [r, ..., 2 * r, ...]
     * the output RLWE ciphertext is:
     * ---> [1, 2, 3,..., N]
    */
    for (size_t i = 0; i < rlwesNum; i++)
        for (size_t j = 0; j < N; j += rlwesNum)
        {
            uint64_t encode = (j + i) + 1;
            // message[i][j + target_row] = encode * 1000000000000;
            message[i][j + target_row] = encode * bsgsDelta;
        }

    intel::hexl::NTT ntts = answerKey.getNTT();

    for (size_t i = 0; i < rlwesNum; i++)
    {
        ntts.ComputeForward(message[i].data(), message[i].data(), 1, 1);
        encrypt(rlwes[i].b.data(), rlwes[i].a.data(), answerKey, message[i].data());
        rlwes[i].setIsNtt(true);
    }

    auto start_b = std::chrono::high_resolution_clock::now();
    RGSWCiphertext query(N, modulus, 2, 0x01 << 20, 0x01 << 18);
    // RGSWCiphertext query(N, modulus, 3, 0x01 << 11, 0x01 << 15);
    query.keyGen(answerKey, target_row, true);
    auto stop_b = std::chrono::high_resolution_clock::now();

    auto glapsed_b = std::chrono::duration_cast<std::chrono::microseconds>(stop_b - start_b);
    std::cout << "encryption costs " << glapsed_b.count() << " us." << std::endl;

    auto start_pk = std::chrono::high_resolution_clock::now();
    // AutoKey autokey(length, modulus, 4, 0, 0x01 << 14);
    AutoKey autokey(length, modulus, 3, 0, 0x01 << 19);
    // AutoKey autokey(length, modulus, 10, 0, 0x01 << 6);
    autokey.keyGen(answerKey, rlwesNum, true);
    auto stop_pk = std::chrono::high_resolution_clock::now();

    auto glapsed_pk = std::chrono::duration_cast<std::chrono::microseconds>(stop_pk - start_pk);
    std::cout << "building packing keys costs " << glapsed_pk.count() << " us." << std::endl;

    // the results
    RlweCiphertext result(N, modulus);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < rlwesNum; i++)
    {
        externalProduct(ext_rlwes[i], rlwes[i], query);
    }

for (size_t i = 0; i < ntimes; i++)
{
    // packingRLWEs(result, rlwes, autokey);
    packingRLWEs(result, ext_rlwes, autokey);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " packing rlwes costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message, result, answerKey);

    multConst(decryptd_message, QInv(rlwesNum, modulus), modulus); // remove term 2^l
    std::cout << std::endl;

    // convert to the original format
    for (auto & me : decryptd_message)
    {
        if (me > modulus / 2)
        {
            me = round(((long double)me - modulus) / bsgsDelta) + bsgsp;
        } else {
            me = round((long double)me / bsgsDelta);
        }
    }

    showLargeVector(decryptd_message, "The decrypted value is ");
}

int main(int argc, char** argv)
{
    srand(time(NULL));
    
    // test_packing_lwes();
    test_packing_rlwes();

    return 0;
}
