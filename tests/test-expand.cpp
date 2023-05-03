#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <math.h>

#include "pir.h"


void test_crt()
{
    std::vector<uint64_t> a(N), b(N);
    std::vector<uint128_t> input(N), output(N);
    uint64_t modulus1 = bigMod;
    uint64_t modulus2 = bigMod2;

    // sample and crt
    sample_random(a, modulus1);
    sample_random(b, modulus2);

    for (size_t i = 0; i < N; i++)
    {
        input[i] = (static_cast<uint128_t>(a[i]) << 10 | b[i]);
        // input[i] = modulus2 + 1;
        a[i] = input[i] % modulus1;
        b[i] = input[i] % modulus2;
    }
    
    // crt inverse
    crt_inv(output, a, b, modulus1, modulus2);

    std::cout << "Done." << std::endl;
}


void copy_value(std::vector<uint64_t>& result, std::vector<uint64_t>& input)
{
    for (size_t i = 0; i < N; i++)
    {
        result[i] = input[i];
    }
}

void add(std::vector<uint64_t>& result, std::vector<uint64_t>& input1,
            std::vector<uint64_t>& input2, uint64_t modulus)
{
    for (size_t i = 0; i < N; i++)
    {
        result[i] = (input1[i] + input2[i]) % modulus;
    }   
}

void mult_factor(std::vector<uint64_t>& result, std::vector<uint64_t>& input,
                    uint64_t neg_factor, uint64_t length, uint64_t modulus)
{
    uint64_t number = length - neg_factor;

    for (size_t i = 0; i < length - number; i++)
    {
        result[i + number] = modulus - input[i];
    }
    for (size_t i = length - number; i < length; i++)
    {
        int32_t index = i + number - length;
        result[index] = input[i];
    }
}

void test_plaintext_expand()
{
    std::vector<uint64_t> input(N, 1);
    std::vector<std::vector<uint64_t> > output(N, std::vector<uint64_t>(N, 0));
    uint64_t modulus = bigMod;
    uint64_t length = N;

    sample_random(input, modulus);

    std::vector<uint64_t> temp1(N), temp2(N);
    
    int32_t l = log2(N);

    copy_value(output[0], input);
    for (size_t j = 0; j < l; j++)
    {
        int32_t index = (0x01 << (l - j)) + 1;
        for (size_t k = 0; k < (0x01 << j); k++)
        {
            mult_factor(temp1, output[k], 0x01 << j, length, modulus);

            automorphic(temp2, output[k], index, modulus);
            add(output[k], output[k], temp2, modulus);
            
            automorphic(temp2, temp1, index, modulus);
            add(output[k + (0x01 << j)], temp1, temp2, modulus);
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        output[i][0] = static_cast<uint128_t>(output[i][0]) * bNinv % bigMod;
    }

    std::cout << "Done." << std::endl;
}

/**
 * @brief test automorphic transform rns variant
 * 
 */
void test_expand()
{
    Secret answerKey(bigMod); // default rlwe secret, store in ntt form
    RlweCiphertext cipher1(N, bigMod);
    RlweCiphertext cipher2(N, bigMod2);
    std::vector<RlweCiphertext> cipher;
    cipher.push_back(cipher1);
    cipher.push_back(cipher2);

    std::vector<RlweCiphertext> expandCipher(2 * N, RlweCiphertext(N, bigMod));
    for (size_t i = 0; i < N; i++)
        expandCipher[i + N].setModulus(bigMod2);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, bigMod);
    for (size_t i = 0; i < N; i++)
        message[i] = 1000000;

    showLargeVector(message, "The message is ");

    // encrypt message
    encrypt_rns(cipher[0], cipher[1], answerKey, message);

    // evaluate automorphic transform
    AutoKeyRNS autokey;
    autokey.keyGen(answerKey, N); // N is for sample errors

    auto start = std::chrono::high_resolution_clock::now();
    evalExpandRNS(expandCipher, cipher, autokey);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "automorphic transform costs " << glapsed.count() << " ms." << std::endl;

    RlweCiphertext result;
    modSwitch(result, cipher[0], cipher[1]);
    // modSwitch(result, expandCipher[0], expandCipher[N]);

    // decrypt message
    decrypt(decryptd_message, result, answerKey);
    showLargeVector(decryptd_message, "The decrypted value is ");
}


int main(int argc, char** argv)
{
    srand(time(NULL));

    // test_crt();
    // test_plaintext_expand();
    test_expand();

    return 0;
}
