#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "pir.h"
#include "utils.h"

void test_decompose()
{
    uint64_t modulus = crtMod;

    std::vector<uint64_t> a(N);
    sample_random(a, modulus);

    int32_t ellnum = 8; // 24;
    uint64_t base = 0; // 0x01 << 10; // 0x01 << 2; // 0x01 << 0;
    uint64_t BBg = 0x01 << 7; // 0x01 << 2;
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(N, 0));
    decompose(dec_a, a, ellnum, base, BBg, modulus);

    check_recontruct(dec_a, a, ellnum, base, BBg, modulus);
}

void test_decompose_bsgs()
{
    uint64_t modulus = crtMod;
    uint64_t bsModulus = bsMod;

    std::vector<uint64_t> a1(N), a2(N);
    sample_random(a1, modulus);
    sample_random(a2, bsModulus);

    int32_t ellnum = 4; // 8;
    uint64_t base = 1024; // 0;
    uint64_t BBg = 0x01 << 20; // 0x01 << 10;
    std::vector<std::vector<uint64_t> > dec_a1(ellnum, std::vector<uint64_t>(N, 0));
    std::vector<std::vector<uint64_t> > dec_a2(ellnum, std::vector<uint64_t>(N, 0));
    decompose_bsgs(dec_a1, dec_a2, a1, a2, ellnum, base, BBg);

    check_recontruct_bsgs(dec_a1, dec_a2, a1, a2, ellnum, base, BBg, modulus * (uint128_t)bsModulus);
}

void test_bgntt()
{
    uint64_t modulus = bsMod;

    std::vector<uint64_t> a(N);
    sample_random8_vector(a.data(), N);
    showLargeVector(a, "input  = ");

#ifdef INTEL_HEXL
    intel::hexl::NTT bsNtts(N, modulus);
    bsNtts.ComputeForward(a.data(), a.data(), 1, 1);
    bsNtts.ComputeInverse(a.data(), a.data(), 1, 1);

    showLargeVector(a, "result = ");
#endif
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
    // uint64_t modulus = bigMod;
    uint64_t modulus = crtMod;
    
    Secret answerKey(modulus); // default rlwe secret, store in ntt form
    RlweCiphertext cipher(N , modulus);
    RlweCiphertext result(N , modulus);

    std::vector<uint64_t> message(N), decryptd_message(N);
    sample_random(message, modulus);
    showLargeVector(message, "The message is ");

    // encrypt message
    encrypt(cipher, answerKey, message);

    // evaluate external product
    // RGSWCiphertext query(N, modulus, 3, 0x01 << 11, 0x01 << 15);
    RGSWCiphertext query(N, modulus, 2, 0x01 << 18, 0x01 << 16);
    auto start_qu = std::chrono::high_resolution_clock::now(); 
    query.keyGen(answerKey, 3, true);
    auto stop_qu = std::chrono::high_resolution_clock::now();

    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query costs " << glapsed_qu.count() << " us." << std::endl;

int32_t ntimes = 1; // 16 * 32;

    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    externalProduct(result, cipher, query);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "each external product costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    auto start_de = std::chrono::high_resolution_clock::now();
    decrypt(decryptd_message, result, answerKey);
    auto stop_de = std::chrono::high_resolution_clock::now();

    auto glapsed_de = std::chrono::duration_cast<std::chrono::microseconds>(stop_de - start_de);
    std::cout << " decrypt costs " << glapsed_de.count() << " us." << std::endl;

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

void test_intel_hexl()
{
    int32_t length = N;
    // uint64_t modulus = bigMod;
    uint64_t modulus = crtMod;

    uint64_t root_of_unity1 = intel::hexl::MinimalPrimitiveRoot(2 * N, crtq1);
    uint64_t root_of_unity2 = intel::hexl::MinimalPrimitiveRoot(2 * N, crtq2);
    uint64_t root_of_unity = crt_compose(root_of_unity1, root_of_unity2);
    std::cout << "root_of_unity = " << root_of_unity << std::endl;
    /***********************     test permutation      *************************/

    int32_t rotate_num = 2;

    std::vector<uint64_t> input(length, 0);
    std::vector<uint64_t> input_ntt(length, 0);
    std::vector<uint64_t> automorphiced(length, 0);
    std::vector<uint64_t> automorphiced_ntt(length, 0);


#ifdef INTEL_HEXL
    // intel::hexl::NTT ntts(length, modulus);
    intel::hexl::NTT ntts(length, modulus, root_of_unity_crt);
    sample_random8_vector(input.data(), length);

    // input_ntt
    ntts.ComputeForward(input_ntt.data(), input.data(), 1, 1);

    automorphic(automorphiced, input, pow_mod(5, rotate_num, 2 * length), modulus); // rotate 2
    ntts.ComputeForward(automorphiced_ntt.data(), automorphiced.data(), 1, 1);

    // automorphic elments has same ntt numbers

/***********
    std::vector<uint64_t> input2(length, 0);
    std::vector<uint64_t> input2_ntt(length, 0);
    input2[1] = 1;
    ntts.ComputeForward(input2_ntt.data(), input2.data(), 1, 1);

    std::cout << "RootOfUnity: " << ntts.GetRootOfUnityPower(0) << std::endl;
    std::cout << "MiniRootOfUnity: " << ntts.GetMinimalRootOfUnity() << std::endl;
    // for (size_t i = 0; i < length; i++)
    for (size_t i = length - 20; i < length; i ++)
    {
        for (size_t j = 1; j < 2 * N; j += 2)
        {
            if (pow_mod(ntts.GetMinimalRootOfUnity(), j, modulus) == input2_ntt[i])
            {
                std::cout << j << ", " << std::endl;
                break;
            }
        }
    }
*************/
    std::vector<int32_t> permutation(length, 0);
    compute_permutation(permutation, rotate_num, length);

    std::vector<uint64_t> result1(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        result1[i] = input_ntt[permutation[i]];
    }

    showLargeVector(automorphiced_ntt, "result  = ");
    showLargeVector(result1, "result1 = ");
    std::cout << " permutation finished." << std::endl;


    /***********************  test permutation matrix  *************************/
    int32_t max_index = 64;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int32_t> > permutations(max_index, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, max_index, length);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "generating permutation matrix costs " << glapsed.count() << " us." << std::endl;

    std::vector<uint64_t> result2(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        result2[i] = input_ntt[permutations[rotate_num][i]];
    }

    showLargeVector(result2, "result2 = ");
    std::cout << " permutation matrix finished." << std::endl;



    /***********************  test encode  *************************/
    /**    
    std::vector<uint64_t> input_single(length, 0);
    std::vector<uint64_t> automorphiced_single(length, 0);
    std::vector<uint64_t> result3(length, 0);

    int32_t encode_index = 3; // rand() & (N / 2 - 1);
    input_single[encode_index] = 111111;
    input_single[N - 1 - encode_index] = 111111;

    ntts.ComputeInverse(input_single.data(), input_single.data(), 1, 1);
    for (size_t i = 1; i < N / 2; i++)
    {
        automorphic(automorphiced_single, input_single, pow_mod(5, i, 2 * length), modulus); // rotate 2
        ntts.ComputeForward(automorphiced_single.data(), automorphiced_single.data(), 1, 1);
        intel::hexl::EltwiseAddMod(result3.data(), result3.data(), automorphiced_single.data(),
                                length, modulus);
    }

    std::cout << " test encode finished." << std::endl;
    **/
#endif
}

void test_bsgs()
{
    int32_t length = N;
    uint64_t modulus = bigMod;
    int32_t query_num = 3;

    std::vector<std::vector<uint64_t> > data(length, std::vector<uint64_t>(length/2, 0));
    sample_database_bsgs(data);
    // we set the colunm of query_num to be particular elements
    for (size_t i = 0; i < length; i++)
    {
        data[i][query_num] = i + 1000000000; // 111111;
    }
    std::vector<uint64_t> result(length, 0);

    // database precomputation
    std::vector<int32_t> query_encode(length, 0); // find encode positon
    compute_query_encode(query_encode, length);
    std::vector<int32_t> query_decode(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        // find_index[hexl_ntt_index[i] >> 0x01] = i;
        query_decode[query_encode[i]] = i;
    }

    std::vector<uint64_t> query(length, 0);
    query[query_encode[query_num]] = 1;
    query[query_encode[N - 1 - query_num]] = 1;

    std::vector<std::vector<uint64_t> > data_ntt(length/2, std::vector<uint64_t>(length, 0));
    for (size_t i = 0; i < N / 2; i++)
    {
        for (size_t j = 0; j < N / 2; j++)
        {
            // diagonal line
            int row = (j - i + N/2) & (N/2 - 1); // the row is always smaller 
            data_ntt[i][j] = data[row][j];

            data_ntt[i][N - 1 - j] = data[N - 1 - row][j];
            // data_ntt[i][j] = data[j][col];
            // data_ntt[i][N - 1 - j] = data[N - 1 - j][col];
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int32_t> > permutations(length/2, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, length/2, length);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "generating permutation matrix costs " << glapsed.count() << " us." << std::endl;

    // rotate 1 ~ N/2
    std::vector<uint64_t> temp1(length, 0);
    std::vector<uint64_t> temp2(length, 0);
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    auto start_ntt = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N / 2; i++)
    {
        // query permutation for data_ntt
        // copy data_ntt[i] to temp
        copy(data_ntt[i].begin(), data_ntt[i].end(), temp1.begin());
        // encode data_ntt[i]
        for (size_t j = 0; j < N; j++)
        {
            temp2[query_encode[j]] = temp1[j];
        }

        // (left) rotation for data_ntt[i]
        // ntts.ComputeInverse(temp2.data(), temp2.data(), 1, 1);
        // automorphic(temp1, temp2, pow_mod(5, i, 2 * length), modulus);
        // ntts.ComputeForward(data_ntt[i].data(), temp1.data(), 1, 1);
        /**/
        for (size_t j = 0; j < N; j++)
        {
            data_ntt[i][j] =  temp2[permutations[i][j]]; // i is rotation number
        }
        /**/
    }
    auto stop_ntt = std::chrono::high_resolution_clock::now();

    auto glapsed_ntt = std::chrono::duration_cast<std::chrono::microseconds>(stop_ntt - start_ntt);
    std::cout << "ntt costs " << glapsed_ntt.count() << " us." << std::endl;


    // matrix vector multiplication ?
    std::vector<uint64_t> query_temp(length, 0);
    for (size_t i = 0; i < N / 2; i++)
    {
        // ntts.ComputeInverse(query_temp.data(), query.data(), 1, 1);
        // automorphic(temp1, query_temp, pow_mod(5, i, 2 * length), modulus);
        // ntts.ComputeForward(query_temp.data(), temp1.data(), 1, 1);
        for (size_t j = 0; j < N; j++)
        {
            query_temp[j] = query[permutations[i][j]]; // i is rotation number
        }

        intel::hexl::EltwiseMultMod(temp1.data(), data_ntt[i].data(), query_temp.data(),
                                    length, modulus, 1);
        intel::hexl::EltwiseAddMod(result.data(), result.data(), temp1.data(),
                                    length, modulus);
    }
#endif 

        // query permutation for data_ntt
        copy(result.begin(), result.end(), temp1.begin());
        for (size_t j = 0; j < N; j++)
        {
            result[query_decode[j]] = temp1[j];
        }

    std::cout << " test bsgs encode finished." << std::endl;
}

void test_mul_error()
{
    int32_t length = N;
    uint64_t modulus = bigMod;

    std::vector<uint64_t> m1(length, 0);
    sample_random(m1, bsgsp);
    for (size_t j = 0; j < N; j++)
    {
        m1[j] = j + 1;
    }

    std::vector<uint64_t> input(length, 0);
    input[0] = 1;
    input[N - 1] = 1;

    // precompute and store in ntt form
    intel::hexl::NTT ntts(N, bigMod);
    intel::hexl::NTT nttp(N, bsgsp);

    nttp.ComputeInverse(m1.data(), m1.data(), 1, 1);
    uint64_t sub = modulus - bsgsp;
    for (size_t j = 0; j < N; j++)
    {
        if (m1[j] > bsgsp / 2)
            m1[j] += sub;
    }
    ntts.ComputeForward(m1.data(), m1.data(), 1, 1);

    nttp.ComputeInverse(input.data(), input.data(), 1, 1);
    for (auto & me : input)
    {
        /**
        if (me > bsgsp / 2)
        {
            me *= bsgsDelta;
            me += (bigMod - bsgsp * bsgsDelta);
        } else {
            me *= bsgsDelta;
        }
        **/
        // Revisiting rounding
        /**/
        if (me > bsgsp / 2)
        {
            me = ((long double)me - bsgsp) * bigMod / bsgsp + bigMod;
        } else {
            me = round((long double)me * bigMod / bsgsp);
        }
        /**/
    }
    ntts.ComputeForward(input.data(), input.data(), 1, 1);

    std::vector<uint64_t> result(length, 0);
    intel::hexl::EltwiseMultMod(result.data(), m1.data(), input.data(),
                                length, modulus, 1);
    ntts.ComputeInverse(result.data(), result.data(), 1, 1);

    for (auto & me : result)
    {
        if (me > bigMod / 2)
        {
            // ((double)me - bigMod) / bsgsDelta is always less than 0, so -.5 for rounding
            long double r = ((long double)me - bigMod) / bsgsDelta - .5;

            // TEST:
            uint64_t test = round(((long double)me - bigMod) / bsgsDelta) + bsgsp;
            if (test != (int64_t)r + bsgsp) std::cout << "error: " << me << std::endl; 

            me = (int64_t)r + bsgsp;
        } else {
            long double r2 = (long double)me / bsgsDelta + 0.5;

            // TEST:
            uint64_t test = round((long double)me / bsgsDelta);
            if (test != (uint64_t)r2) std::cout << "error: " << me << std::endl; 

            me = r2;
        }
    }

    nttp.ComputeForward(result.data(), result.data(), 1, 1);

    std::cout << " test multiplication error finished." << std::endl;
}

void test_nttp()
{   
    std::vector<uint64_t> input(N, bsgsp);
    std::vector<uint64_t> output(N, bsgsp);
    std::vector<uint64_t> input2(N, bsgsp);

    for (size_t i = 0; i < N; i++)
    {
        input[i] = i + 1;
    }

    intel::hexl::NTT nttp(N, bsgsp);
    nttp.ComputeForward(output.data(), input.data(), 1, 1);
    nttp.ComputeInverse(input2.data(), output.data(), 1, 1);

    std::cout << " test nttp finished." << std::endl;
}

void test_two_ntts()
{
    size_t length = 2048;

    // scratch = (uint64_t *)malloc(2 * poly_len * sizeof(uint64_t));

    std::vector<uint64_t> input1(length, 0);
    std::vector<uint64_t> input2(2 * length, 0);
    std::vector<uint64_t> input3(length, 0);

    // sample_random8_vector(input1.data(), length);
    sample_random(input1, bigMod);

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, bigMod);

    auto start = std::chrono::high_resolution_clock::now();
    ntts.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts.ComputeInverse(input1.data(), input1.data(), 1, 1);
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt costs " << glapsed.count() << " us." << std::endl;

    uint64_t modulus2 = 72057594038149121; // pow(2, 56)
    intel::hexl::NTT ntts2(length, modulus2);

    start = std::chrono::high_resolution_clock::now();
    ntts2.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts2.ComputeInverse(input1.data(), input1.data(), 1, 1);
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (modulus 56 bits) costs " << glapsed.count() << " us." << std::endl;
#endif
    auto start_crt = std::chrono::high_resolution_clock::now();
    computeForward(input2.data(), input1.data());
    computeInverse(input3.data(), input2.data());
    auto stop_crt  = std::chrono::high_resolution_clock::now();

    auto glapsed_crt = std::chrono::duration_cast<std::chrono::microseconds>(stop_crt - start_crt);
    std::cout << "crt ntt costs " << glapsed_crt.count() << " us." << std::endl;

    showLargeVector(input1, "input1 = ");
    showLargeVector(input3, "input3 = ");
}

void test_two_ntts4096()
{
    size_t length = 4096;

    // scratch = (uint64_t *)malloc(2 * poly_len * sizeof(uint64_t));

    std::vector<uint64_t> input1(length, 0);
    std::vector<uint64_t> input2(2 * length, 0);
    std::vector<uint64_t> input3(length, 0);

    // sample_random8_vector(input1.data(), length);
    sample_random(input1, bigMod);

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, bigMod);

    auto start = std::chrono::high_resolution_clock::now();
    ntts.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts.ComputeInverse(input1.data(), input1.data(), 1, 1);
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt costs " << glapsed.count() << " us." << std::endl;

    uint64_t modulus2 = 72057594038149121; // pow(2, 56)
    intel::hexl::NTT ntts2(length, modulus2);

    start = std::chrono::high_resolution_clock::now();
    ntts2.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts2.ComputeInverse(input1.data(), input1.data(), 1, 1);
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (modulus 56 bits) costs " << glapsed.count() << " us." << std::endl;
#endif
}

void computeMiniRoot()
{
    uint64_t length = 4096;
    uint64_t modulus1 = crtq1;
    uint64_t modulus2 = crtq2;
    uint64_t modulus = modulus1 * modulus2;
    uint64_t root_of_unity1 = intel::hexl::MinimalPrimitiveRoot(2 * length, modulus1);
    uint64_t root_of_unity2 = intel::hexl::MinimalPrimitiveRoot(2 * length, modulus2);

    uint64_t mini_root_of_unity = modulus;
    uint64_t current_root_of_unity1 = root_of_unity1;
    uint64_t current_root_of_unity2 = root_of_unity2;
    uint64_t root_of_unity1_powof2 = root_of_unity1 * root_of_unity1 % modulus1;
    uint64_t root_of_unity2_powof2 = root_of_unity2 * root_of_unity2 % modulus2;
    for (size_t i = 0; i < length; i++)
    { 
        for (size_t j = 0; j < length; j++)
        {
        uint64_t temp = crt_compose(current_root_of_unity1, current_root_of_unity2);
        if (temp < mini_root_of_unity) mini_root_of_unity = temp;
        current_root_of_unity2 = current_root_of_unity2 * root_of_unity2_powof2 % modulus2;
        }
        current_root_of_unity1 = current_root_of_unity1 * root_of_unity1_powof2 % modulus1;
  }

  std::cout << "mini_root_of_unity: " << mini_root_of_unity << std::endl;
  // mini_root_of_unity = 10297991595
}

// the minimum root of unity is better than default inverse crt?
void test_ntt_crt()
{
    size_t length = 4096;
    uint64_t modulus = crtMod;

    std::vector<uint64_t> input1(length, 0);

    sample_random8_vector(input1.data(), length);
    // sample_random(input1, modulus);

#ifdef INTEL_HEXL
    uint64_t root_of_unity = 10297991595; //minimum root of unity
    intel::hexl::NTT ntts(length, modulus, root_of_unity);

int ntimes = 10;
    auto start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < ntimes; i++)
{
    ntts.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts.ComputeInverse(input1.data(), input1.data(), 1, 1);
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt (mini root) costs " << glapsed.count() << " us." << std::endl;

    // default root of unity
    uint64_t root_of_unity1 = intel::hexl::MinimalPrimitiveRoot(2 * N, crtq1);
    uint64_t root_of_unity2 = intel::hexl::MinimalPrimitiveRoot(2 * N, crtq2);
    root_of_unity = crt_compose(root_of_unity1, root_of_unity2);

    intel::hexl::NTT ntts2(length, modulus, root_of_unity);

    start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    ntts2.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts2.ComputeInverse(input1.data(), input1.data(), 1, 1);
}
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "hexl ntt costs " << glapsed.count() << " us." << std::endl;
#endif
}

// compute root of unity of crtBaMod
void test_compute_root_of_unity()
{
    size_t length = 4096;
    uint64_t modulus = crtBaMod;

    std::vector<uint64_t> input1(length, 0);

    sample_random8_vector(input1.data(), length);
    // sample_random(input1, modulus);

#ifdef INTEL_HEXL

    // uint64_t root_of_unity1 = intel::hexl::MinimalPrimitiveRoot(2 * N, bsMod);
    // uint64_t root_of_unity2 = intel::hexl::MinimalPrimitiveRoot(2 * N, auxMod);
    uint64_t root_of_unity = 1486445687605966UL;

    intel::hexl::NTT ntts2(length, modulus, root_of_unity);

    ntts2.ComputeForward(input1.data(), input1.data(), 1, 1);
    ntts2.ComputeInverse(input1.data(), input1.data(), 1, 1);

    std::cout << "finished." << std::endl;
#endif
}

void test_inversedb()
{
    uint64_t length = N;
    std::vector<uint64_t> message(length, 0);

    for (size_t i = 0; i < length; i++)
        message[i] = i + 1;
    
    // encode
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length);

    std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_encode[j]] = temp[j];
    }

    intel::hexl::NTT nttp(length, bsgsp);
    nttp.ComputeInverse(message.data(), message.data(), 1, 1);

    // decode
    nttp.ComputeForward(message.data(), message.data(), 1, 1);

    std::vector<int32_t> query_decode(length, 0);
    compute_query_decode(query_decode, length);

    // std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_decode[j]] = temp[j];
    }

    showLargeVector(message, "Decrypted message:");
}

void test_forwarddb()
{
    uint64_t length = N;
    std::vector<uint64_t> message(length, 0);

    for (size_t i = 0; i < length; i++)
        message[i] = i + 1;
    
    // encode
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length);

    std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_encode[j]] = temp[j];
    }

    intel::hexl::NTT nttp(length, bsgsp);
    nttp.ComputeInverse(message.data(), message.data(), 1, 1);

    // decode
    nttp.ComputeForward(message.data(), message.data(), 1, 1);

    std::vector<int32_t> query_decode(length, 0);
    compute_query_decode(query_decode, length);

    // std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_decode[j]] = temp[j];
    }

    showLargeVector(message, "Decrypted message:");
}

int main(int argc, char** argv)
{
    // srand(time(NULL));
    
    // test_decompose();
    // test_decompose_bsgs();
    // test_bgntt();
    // test_lwe_encrypt();
    // test_rlwe_encrypt();
    // test_automorphic_transform();
    // test_automorphic_transform_rns();

    test_external_product();

    // test_dummy_ksKeyGen();
    // test_recover();

    // test_intel_hexl();
    // test_bsgs();
    // test_mul_error();
    // test_nttp();
    // test_two_ntts(); // hexl or crt?
    // test_two_ntts4096(); // 50 or 56 bits?
    // computeMiniRoot();
    // test_ntt_crt(); // mini root of unity of default root of unity
    // test_compute_root_of_unity();
    // test_inversedb();


#if defined(__AVX2__)
    std::cout << "have avx2" << std::endl;
#else
    std::cout << "warning: no avx2" << std::endl;
#endif

    return 0;
}
