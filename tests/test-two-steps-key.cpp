#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"
#include "crt.h"

void test_two_steps_bsgs_rns()
{
    Secret queryKey(crtMod, false); // stored in coefficient form

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    sample_database_bsgs(data);

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    int32_t N1 = 64; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "N1: " << N1 << ", N2: " << N2 << std::endl;

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, crtMod, N1);

    size_t num_words = N * N / 2;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    database_tocrt(datacrt, data_ntt, N1);

    // auto start_qu = std::chrono::high_resolution_clock::now();
    std::vector<RlweCiphertext> query1(1, RlweCiphertext(N, crtMod));
    query1.push_back(RlweCiphertext(N, bsMod));
    query_bsgs_rns(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGSRNS autoKey(N, crtMod, bsMod);
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i < N1; i++)
    {
        indexLists.push_back(pow_mod(5, i, 2 * N));
    }
    autoKey.keyGen(queryKey, indexLists, BabyStep);

    indexLists.clear();
    for (size_t i = 1; i < N2; i++)
    {
        indexLists.push_back(pow_mod(5, N1 * i, 2 * N));
    }
    autoKey.keyGen(queryKey, indexLists, GaintStep);

    // compute permutation matrix
    std::vector<std::vector<int32_t> > permutations(N1, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, N1, length);

    // the results
    RlweCiphertext result(N, crtMod);
    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs_rns_crt(result, query1, datacrt, autoKey, N1, permutations);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " matrix vector multiplication bsgs costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt_bsgs(decryptd_message, result, queryKey);

    std::cout << std::endl;

    std::cout << "the recovered value is " << decryptd_message[row] << std::endl;
    showLargeVector(decryptd_message, "result = ");
}

void generate_offline_key(AutoKeyBSGSRNS& offlineKey1, AutoKeyBSGSRNS& offlineKey2,
                          AutoKeyBSGSRNS& autoKey1, AutoKeyBSGSRNS& autoKey2,
                          Secret& queryKey, int32_t N1)
{
    std::vector<int32_t> indexLists;
    indexLists.push_back(pow_mod(5, 1, 2 * N));
    offlineKey1.keyGen(queryKey, indexLists, BabyStep);
    indexLists.clear();

    indexLists.push_back(pow_mod(5, N1 * 1, 2 * N));
    offlineKey2.keyGen(queryKey, indexLists, BabyStep);

    // the results
    // AutoKeyBSGSRNS autoKey(N, crtMod, bsMod);
    // the auxMod here is for modulus switching
    indexLists.clear();
    indexLists.push_back(pow_mod(5, 1, 2 * N));
    autoKey1.keyGen(queryKey, indexLists, BabyStep);
    indexLists.clear();

    indexLists.push_back(pow_mod(5, N1 * 1, 2 * N));
    autoKey2.keyGen(queryKey, indexLists, BabyStep);
}

/**
 * @brief generate autoKey from the offline keys
 * 
*/
void test_generate_autokey()
{
    Secret queryKey(crtMod, false); // stored in coefficient form

    int32_t N1 = 64; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "N1: " << N1 << ", N2: " << N2 << std::endl;

    AutoKeyBSGSRNS offlineKey1(N, crtMod, crtBaMod, 5, 0x01 << 18, 0x01 << 18);
    AutoKeyBSGSRNS offlineKey2(N, crtMod, auxMod, 5, 0x01 << 10, 0x01 << 14);
    AutoKeyBSGSRNS autoKey1(N, crtMod, crtBaMod, 3, (0x01 << 20) * auxMod, 0x01 << 20);
    AutoKeyBSGSRNS autoKey2(N, crtMod, auxMod, 3, (0x01 << 11) * auxMod, 0x01 << 15);
    generate_offline_key(offlineKey1, offlineKey2, autoKey1, autoKey2, queryKey, N1);

    AutoKeyBSGSRNS result(N, crtMod, bsMod);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    genAutoKeyFromOffline(result, autoKey1, autoKey2, offlineKey1, offlineKey2, N1);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " generate auto keys costs " << glapsed.count() << " us." << std::endl;
}

int main(int argc, char** argv)
{
    // srand(time(NULL));

    // test_two_steps_bsgs_rns();
    test_generate_autokey();

    return 0;
}
