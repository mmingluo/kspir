#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <math.h>

#include "pir.h"
#include "utils.h"

void test_decompose()
{
    std::vector<uint64_t> a(N);
    sample_random(a, bigMod);

    std::vector<std::vector<uint64_t> > dec_a(ell, std::vector<uint64_t>(N, 0));
    decompose(dec_a, a, ell, Base, Bg);

    check_recontruct(dec_a, a, ell, Base, Bg);
}


/**
 * @brief test lwe encryption
 * 
 */
void test_lwe_encrypt()
{
    Secret answerKey(LWE, bigMod);
    RlweCiphertext lwes(N , bigMod);

    uint64_t message = 123456789;

    // encrypt message
    auto start = std::chrono::high_resolution_clock::now();
    encrypt(lwes.b.data(), lwes.a.data(), answerKey, message);
    auto stop = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "encrypt costs " << glapsed.count() << " us." << std::endl;
    
    lweToRlwe(lwes.a);

    // decrypt message
    start = std::chrono::high_resolution_clock::now();
    // uint64_t decryptd_message = decrypt(lwes.b.data(), lwes.a.data(), answerKey);
    std::vector<uint64_t> decryptd_message(N);
    decrypt(decryptd_message, lwes, answerKey);

    stop = std::chrono::high_resolution_clock::now();
    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "decrypt costs " << glapsed.count() << " us." << std::endl;

    std::cout << "The decrypted value is " << decryptd_message[0] << std::endl;
}

/**
 * @brief test rlwe encryption
 * 
 */
void test_rlwe_encrypt()
{
    Secret answerKey(bigMod); // default rlwe secret, store in ntt form
    RlweCiphertext rlwes(N , bigMod);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, bigMod);
    showLargeVector(message, "The message is ");

    // encrypt message
    auto start = std::chrono::high_resolution_clock::now();
    encrypt(rlwes, answerKey, message);
    auto stop = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "encrypt costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    start = std::chrono::high_resolution_clock::now();
    decrypt(decryptd_message, rlwes, answerKey);
    stop = std::chrono::high_resolution_clock::now();
    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "decrypt costs " << glapsed.count() << " us." << std::endl;

    showLargeVector(decryptd_message, "The decrypted value is ");
}

/**
 * @brief test automorphic transform
 * 
 */
void test_automorphic_transform()
{
    Secret answerKey(bigMod); // default rlwe secret, store in ntt form
    RlweCiphertext cipher(N , bigMod);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, bigMod);

    showLargeVector(message, "The message is ");

    // encrypt message
    encrypt(cipher, answerKey, message);

    // evaluate automorphic transform
    AutoKey autokey;
    autokey.keyGen(answerKey, 0x01 << 5);

    auto start = std::chrono::high_resolution_clock::now();
    evalAuto(cipher, 3, autokey);
    // sevalAuto(cipher, 5, autokey);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "automorphic transform costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message, cipher, answerKey);
    showLargeVector(decryptd_message, "The decrypted value is ");
}

/**
 * @brief test automorphic transform rns variant
 * 
 */
void test_automorphic_transform_rns()
{
    Secret answerKey(bigMod); // default rlwe secret, store in ntt form
    RlweCiphertext cipher1(N, bigMod);
    RlweCiphertext cipher2(N, bigMod2);
    std::vector<RlweCiphertext> cipher;
    cipher.push_back(cipher1);
    cipher.push_back(cipher2);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, bigMod);
    // for (size_t i = 0; i < N; i++)
    //    message[i] = 1000000;

    showLargeVector(message, "The message is ");

    // encrypt message
    encrypt_rns(cipher[0], cipher[1], answerKey, message);

    // evaluate automorphic transform
    AutoKeyRNS autokey;
    autokey.keyGen(answerKey, 0x01 << 3);

    auto start = std::chrono::high_resolution_clock::now();
    evalAutoRNS(cipher[0], cipher[1], 3, autokey);
    // evalTrNnRNS(cipher, 5, autokey);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "automorphic transform costs " << glapsed.count() << " us." << std::endl;

    RlweCiphertext result;
    modSwitch(result, cipher[0], cipher[1]);

    // decrypt message
    decrypt(decryptd_message, result, answerKey);
    showLargeVector(decryptd_message, "The decrypted value is ");
}


/**
 * @brief test external product
 * 
 */
void test_external_product()
{
    Secret answerKey(bigMod); // default rlwe secret, store in ntt form
    RlweCiphertext cipher(N , bigMod);
    RlweCiphertext result(N , bigMod);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, bigMod);
    showLargeVector(message, "The message is ");

    // encrypt message
    encrypt(cipher, answerKey, message);

    // evaluate external product
    RGSWCiphertext query(N, bigMod);
    query.keyGen(answerKey, 3, true);

    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < 1; i++)
{
    externalProduct(result, cipher, query);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "each external product costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message, result, answerKey);
    showLargeVector(decryptd_message, "The decrypted value is ");
}

void test_recover()
{
    Secret queryKey(bigMod);
    
    uint64_t col = 1;

    RlweCiphertext query1(N, bigMod);
    query(query1, queryKey, col);

    std::vector<uint64_t> message(N);
    message[col] = 999;

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();

    ntts.ComputeForward(query1.a.data(), query1.a.data(), 1, 1);
    ntts.ComputeForward(message.data(), message.data(), 1, 1);

    intel::hexl::EltwiseMultMod(query1.b.data(), query1.b.data(), message.data(), N,
                bigMod, 1);
    intel::hexl::EltwiseMultMod(query1.a.data(), query1.a.data(), message.data(), N,
                bigMod, 1);
    query1.setIsNtt(true);

#endif

    // recover
    std::vector<uint64_t> decrypted(N);
    recover(decrypted, query1, queryKey);

    std::cout << "the recovered value is " << decrypted[0] << std::endl;

}

void test_dummy_ksKeyGen()
{   
    Secret queryKey(bigMod), answerKey(bigMod);

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N, 0));
    sample_database(data);

    // std::vector<uint64_t> a(N), b(N);
    
    uint64_t row = 123, col = 254;

    RlweCiphertext query1(N, bigMod);
    query(query1, queryKey, col);
    RGSWCiphertext query2(N, bigMod);
    query2.keyGen(answerKey, row, true); // note there is reverse

    std::cout << "the wanted message is " << data[row][col] << std::endl;

    database_tontt(data);

    // build keyswitch key
    std::vector<std::vector<uint64_t> > ks(2 * ell * N, std::vector<uint64_t>(N, 0));
    dummy_ksKeyGen(ks, queryKey, answerKey, data);

    RlweCiphertext kskHint(N, bigMod);
    RlweCiphertext result(N, bigMod);

    ksKey_hint(kskHint.b, kskHint.a, query1.a, ks);

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();
    ntts.ComputeInverse(query1.b.data(), query1.b.data(), 1, 1);

    std::vector<uint64_t> decryptd_message(N);
    decrypt(decryptd_message, query1, queryKey);
    showLargeVector(decryptd_message, "The decrypted value is ");

    ntts.ComputeForward(ks[0].data(), ks[0].data(), 1, 1);
    ntts.ComputeForward(ks[1].data(), ks[1].data(), 1, 1);
    decrypt(decryptd_message.data(), ks[0].data(), ks[1].data(), answerKey);
    showLargeVector(decryptd_message, "The decrypted value is ");
#endif
}

int main(int argc, char** argv)
{
    // srand(time(NULL));
    
    // test_decompose();
    // test_lwe_encrypt();
    // test_rlwe_encrypt();
    // test_automorphic_transform();
    test_automorphic_transform_rns();

    // test_external_product();

    // test_dummy_ksKeyGen();
    // test_recover();

    return 0;
}
