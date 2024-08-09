#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"
#include "crt.h"

// 32 MB
void test_two_steps()
{
    Secret queryKey(bigMod);

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    sample_database_bsgs(data);

    // we let ntt(data) to be small
    /**
    std::vector<uint64_t> temp(N, 0);
    intel::hexl::NTT ntts(N, bigMod);
    for (size_t i = 0; i < N / 2; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            temp[j] = rand() & 0x01;
        }
        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
        for (size_t j = 0; j < N / 2; j++)
        {
            // diagonal line
            int row = (j - i + N/2) & (N/2 - 1); // the row is always smaller 
            result[i][j] = data[row][j];
            result[i][N - 1 - j] = data[N - 1 - row][j];
        }
    }
    **/

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, bigMod);

    // auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, bigMod);
    query_bsgs(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGS autoKey;
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i < N/2; i++)
    // for (size_t i = 1; i < 50; i++)
    {
        indexLists.push_back(pow_mod(5, i, 2 * N));
    }
    autoKey.keyGen(queryKey, indexLists);

    // the results
    RlweCiphertext result;
    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul(result, query1, data_ntt, autoKey);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " matrix vector multiplication costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt_bsgs(decryptd_message, result, queryKey);

    std::cout << std::endl;

    std::cout << "the recovered value is " << decryptd_message[row] << std::endl;
    showLargeVector(decryptd_message, "result = ");
}

// 32 MB
void test_two_steps_bsgs()
{
    Secret queryKey(bigMod);

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    sample_database_bsgs(data);

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    int32_t N1 = 64; // 64;
    int32_t N2 = N / 2 / N1;

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, bigMod, N1);

    // auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, bigMod);
    query_bsgs(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGS autoKey;
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

    // the results
    RlweCiphertext result;
    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs(result, query1, data_ntt, autoKey, N1);
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

// same as above, but using the crtmod
void test_two_steps_bsgs_crtmod()
{
    Secret queryKey(crtMod);

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    sample_database_bsgs(data);

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    int32_t N1 = 64; // 64;
    int32_t N2 = N / 2 / N1;

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, crtMod, N1);

    // auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, crtMod);
    query_bsgs(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGS autoKey(N, crtMod);
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

    // the results
    RlweCiphertext result(N, crtMod);
    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs(result, query1, data_ntt, autoKey, N1);
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


// 32 MB
void test_two_steps_bsgs_crt()
{
    Secret queryKey(crtMod);

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    sample_database_bsgs(data);

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    int32_t N1 = 32; // 64; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "N1: " << N1 << ", N2: " << N2 << std::endl;

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, crtMod, N1);

    size_t num_words = N * N / 2;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    database_tocrt(datacrt, data_ntt, N1);

    // auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, crtMod);
    query_bsgs(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGS autoKey(N, crtMod);
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

    // the results
    RlweCiphertext result(N, crtMod);
    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs_crt(result, query1, datacrt, autoKey, N1);
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

    int32_t N1 = 128; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "N1: " << N1 << ", N2: " << N2 << std::endl;

auto start_prep = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));
    database_tobsgsntt(data_ntt, data, crtMod, N1);

    size_t num_words = N * N / 2;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    database_tocrt(datacrt, data_ntt, N1);

auto stop_prep = std::chrono::high_resolution_clock::now();
    auto glapsed_prep = std::chrono::duration_cast<std::chrono::milliseconds>(stop_prep - start_prep);
    std::cout << " server preprocessing costs " << glapsed_prep.count() << " ms." << std::endl;

    // auto start_qu = std::chrono::high_resolution_clock::now();
    std::vector<RlweCiphertext> query1(1, RlweCiphertext(N, crtMod));
    query1.push_back(RlweCiphertext(N, bsMod));
    query_bsgs_rns(query1, queryKey, col);

    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGSRNS autoKey(N, crtMod, bsMod);
    // AutoKeyBSGSRNS autoKey(N, crtMod, bsMod, 7, 0x01 << 10, 0x01 << 10);
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i <= N1 / 2; i++)
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
    std::vector<std::vector<int32_t> > permutations(N1 + 1, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, N1 + 1, length);

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

void test_two_steps_bsgs_rns_large()
{
    // r = 16     : 256 MB
    // r = 64     :   1 GB
    // r = 64 * 8 :   8 GB 

    int32_t r = 64 * 8; // 256 MB
    
    Secret queryKey(crtMod, false); // stored in coefficient form

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    int32_t N1 = 128; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "r: " << r << ", N1: " << N1 << ", N2: " << N2 << std::endl;

    size_t num_words = N * N / 2;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words * r);

    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));


std::chrono::milliseconds glapsed_prep = std::chrono::milliseconds(0);

for (size_t k = 0; k < r; k++)
{
    // std::vector<std::vector<uint64_t> > data_temp(data.begin(), data.begin() + N);
    sample_database_bsgs(data);

    for (size_t i = 0; i < N; i++)
        data[i][col] = i + 1;
    // std::cout << "the wanted message is " << data[row][col] << std::endl;
auto start_prep = std::chrono::high_resolution_clock::now();
    database_tobsgsntt(data_ntt, data, crtMod, N1);
    database_tocrt(datacrt + num_words * k, data_ntt, N1);
auto stop_prep = std::chrono::high_resolution_clock::now();
glapsed_prep += std::chrono::duration_cast<std::chrono::milliseconds>(stop_prep - start_prep);
}

    std::cout << " server preprocessing costs " << glapsed_prep.count() << " ms." << std::endl;


    std::vector<RlweCiphertext> query1(1, RlweCiphertext(N, crtMod));
    query1.push_back(RlweCiphertext(N, bsMod));
    auto start_qu = std::chrono::high_resolution_clock::now();
    query_bsgs_rns(query1, queryKey, col);
    auto stop_qu = std::chrono::high_resolution_clock::now();

    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query costs " << glapsed_qu.count() << " us." << std::endl;


    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGSRNS autoKey(N, crtMod, bsMod);
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i <= N1 / 2; i++)
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
    // RlweCiphertext result(N, crtMod);
    std::vector<RlweCiphertext> result(r, RlweCiphertext(N, crtMod));

    std::vector<uint64_t> decryptd_message(length);

int ntimes = 1;

    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs_rns_crt_large(result, query1, datacrt, autoKey, N1, permutations, r);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " matrix vector multiplication bsgs costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    auto start_de = std::chrono::high_resolution_clock::now();
    decrypt_bsgs(decryptd_message, result[0], queryKey);
    auto stop_de = std::chrono::high_resolution_clock::now();

    auto glapsed_de = std::chrono::duration_cast<std::chrono::microseconds>(stop_de - start_de);
    std::cout << " decrypt costs " << glapsed_de.count() << " us." << std::endl;


    std::cout << std::endl;

    std::cout << "the recovered value is " << decryptd_message[row] << std::endl;
    showLargeVector(decryptd_message, "result = ");
}

int main(int argc, char** argv)
{
    // srand(time(NULL));
    
    // test_two_steps();
    // test_two_steps_bsgs();
    // test_two_steps_bsgs_crtmod(); // basic test but for crtMod
    // test_two_steps_bsgs_crt();
    // test_two_steps_bsgs_rns();
    test_two_steps_bsgs_rns_large(); // packing r basic database

    return 0;
}
