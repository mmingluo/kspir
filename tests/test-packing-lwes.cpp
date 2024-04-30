#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"

void test_packing_lwes()
{
    int32_t l = 4; // 6; // 5
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
    // AutoKey autokey(length, modulus, 3, 0x01 << 14, 0x01 << 12);
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


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_packing_lwes();

    return 0;
}
