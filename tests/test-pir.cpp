
#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"
#include "crt.h"

void test_pir()
{
    // r = 16     : 256 MB
    // r = 64     :   1 GB
    // r = 64 * 8 :   8 GB 

    int32_t r = 16; // 256 MB
    
    Secret queryKey(crtMod, false); // stored in coefficient form

    // uint32_t row = rand() & (N - 1);
    uint32_t target_col = rand() % (N/2); // 123;
    uint32_t target_packing = rand() % r;


    uint64_t modulus = crtMod;
    int32_t N1 = 128; // 128; // 64;
    int32_t N2 = N / 2 / N1;
    std::cout << "Packing number: " << r << std::endl;
    std::cout << "Database configuration: " << r * N/2 << " * 8 KB, ";
    std::cout << "total database size " << 16 * r << " MB" << std::endl;
    std::cout << "BSGS parameters: (N1: " << N1 << ", N2: " << N2 << ")" << std::endl << std::endl;
    std::cout << "target_col: " << target_col << ", target_packing: " << target_packing << std::endl << std::endl;

    size_t num_words = N * N / 2;
    uint64_t* datacrt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words * r);

    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N/2, 0));
    std::vector<std::vector<uint64_t> > data_ntt(N/2, std::vector<uint64_t>(N, 0));


    /******** sample database ********/
    std::chrono::milliseconds glapsed_prep = std::chrono::milliseconds(0);
    std::vector<uint64_t> input_record(N); // store just for encoding
for (size_t k = 0; k < r; k++)
{
    sample_database_bsgs(data);

    auto start_prep = std::chrono::high_resolution_clock::now();
    // insert [1, 2, 3, 4,..., N] in the target position
    for (size_t i = 0; i < N; i += r) data[i + target_packing][target_col] = i + 1 + k;
    
    for (size_t i = 0; i < N; i++) input_record[i] = data[i][target_col];
    inverse_encode(input_record);
    for (size_t i = 0; i < N; i++) data[i][target_col] = input_record[i];

    // preprocessing database
    database_tobsgsntt(data_ntt, data, crtMod, N1);
    database_tocrt(datacrt + num_words * k, data_ntt, N1);
    auto stop_prep = std::chrono::high_resolution_clock::now();
    glapsed_prep += std::chrono::duration_cast<std::chrono::milliseconds>(stop_prep - start_prep);
}
    std::cout << " server preprocessing costs " << glapsed_prep.count() << " ms." << std::endl;


    /******** build query ********/
    auto start_qu = std::chrono::high_resolution_clock::now();
    std::vector<RlweCiphertext> query1(1, RlweCiphertext(N, crtMod));
    query1.push_back(RlweCiphertext(N, bsMod));
    RGSWCiphertext queryGsw(N, modulus, 2, 0x01 << 20, 0x01 << 18);
    // RGSWCiphertext queryGsw(N, modulus, 3, 0x01 << 11, 0x01 << 15);

    query_bsgs_rns(query1, queryKey, target_col);
    queryGsw.keyGen(queryKey, target_packing, true);
    auto stop_qu = std::chrono::high_resolution_clock::now();

    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query costs " << glapsed_qu.count() << " us." << std::endl;


    uint64_t length = queryKey.getLength();
    // uint64_t moudlus = answerKey.getModulus();

    AutoKeyBSGSRNS autoKey(N, crtMod, bsMod); // key-switching keys for BSGS algorithm
    autoKey.bsgsKeyGen(queryKey, N1);

    AutoKey packingKey(length, modulus, 4, 0, 0x01 << 14); // key-switching keys for packing algorithm
    packingKey.keyGen(queryKey, r, true);

    // compute permutation matrix
    std::vector<std::vector<int32_t> > permutations(N1, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, N1, length);

    // the results
    std::vector<RlweCiphertext> result(r, RlweCiphertext(N, crtMod));
    RlweCiphertext result_output(N, modulus);

    std::vector<uint64_t> decryptd_message(length);

    std::vector<RlweCiphertext> ext_rlwes;
    for (size_t i = 0; i < r; i++)
        ext_rlwes.push_back(RlweCiphertext(N, modulus));


    /******** server response ********/
int ntimes = 1;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    matrix_vector_mul_bsgs_rns_crt_large(result, query1, datacrt, autoKey, N1, permutations, r); // the first dimension folding
    for (size_t j = 0; j < r; j++)
        externalProduct(ext_rlwes[j], result[j], queryGsw); // the first dimension folding
    packingRLWEs(result_output, ext_rlwes, packingKey); // the packing algorithm
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    /******** decrypt message ********/
    auto start_de = std::chrono::high_resolution_clock::now();
    decrypt_bsgs_total(decryptd_message, result_output, queryKey, r);
    auto stop_de = std::chrono::high_resolution_clock::now();

    auto glapsed_de = std::chrono::duration_cast<std::chrono::microseconds>(stop_de - start_de);
    std::cout << " decrypt costs " << glapsed_de.count() << " us." << std::endl;

    // output result
    std::cout << std::endl;
    // std::cout << "the recovered value is " << decryptd_message[row] << std::endl;
    showLargeVector(decryptd_message, "the recovered result = ");
}


int main(int argc, char** argv)
{
    // srand(time(NULL));

    test_pir();

    return 0;
}
