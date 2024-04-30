#include <stdint.h>
#include <iostream>
#include <chrono>
#include <time.h>
#include <math.h>
#include <thread>
#include <cassert>

#include "params.h"
#include "answer.h"
#include "ntt.h"
#include "lwe.h"
#include "samples.h"
#include "secret.h"
#include "utils.h"
#include "encrypt.h"


#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

// #define DATA_INDEPENDENT_HANDLE
// #define MISSING_ADD

void decompose_a(uint64_t** dec_a, uint64_t* a)
{
    /*
    dec_a = new uint64_t*[N];
    for (size_t i = 0; i < N; i++)
        dec_a[i] = new uint64_t[l];
    */

    // decompose PVW LWE a

    for (size_t i = 0; i < N; i++) {
        int64_t r = a[i];
        if (r > mod/2)
            r -= mod/2;
        for (size_t j = 0; j < ell; j++) {
            // TODO:
            dec_a[i][0] = r % 256;
            dec_a[i][1] = (int64_t)round((double)r / 256) % 256;
            dec_a[i][2] = round((double)r / 65536);
        }
    }
}

#ifndef INTEL_HEXL
/**
 * @brief build KSkey for each query 
 * 
 * @param b_vec 
 * @param rlwe_a 
 * @param a 
 * @param ks 
 */
void ksKey_hint2(uint64_t* b_rlwe, uint64_t* a_rlwe, uint64_t* a, uint64_t** ks)
{
    uint64_t** dec_a;
    dec_a = new uint64_t*[N];
    for (size_t i = 0; i < N; i++)
        dec_a[i] = new uint64_t[ell];

    decompose_a(dec_a, a);
    for (size_t i = 0; i < N; i++) // 1, X, X^2, X^3, ...
    {
        for (size_t j = 0; j < N; j++)
        {
            for (size_t k = 0; k < ell; k++)
            {
                int index = 2 * ell * j + 2 * k;
                // rlwe_a[i] -= dec_a[j][k] * ks[index][i];
                // b_vec[i] -= dec_a[j][k] * ks[index + 1][i];
                b_rlwe[i] -= barret_reduce((__int128_t)dec_a[j][k] * ks[index][i]);
                a_rlwe[i] -= barret_reduce((__int128_t)dec_a[j][k] * ks[index + 1][i]);
            }
        }
    }

    for (size_t i = 0; i < N; i++)
        delete[] dec_a[i];
    delete[] dec_a;    
}
#endif

#ifdef INTEL_HEXL
int ksKey_hint(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& ks)
{
    lweToRlwe(a); // TODO: check if it's equal to rlweTolwe
    std::vector<std::vector<uint64_t> > dec_a(ell, std::vector<uint64_t>(N, 0));
    decompose(dec_a, a, ell, Base, Bg);

    // check_recontruct(dec_a, a, ell, Base, Bg);

    for (size_t j = 0; j < N; j++)
    {
        for (size_t k = 0; k < ell; k++)
        {
            int index = 2 * ell * j + 2 * k;

            // b_rlwe and a_rlwe are in coefficient form
            intel::hexl::EltwiseFMAMod(b_rlwe.data(), ks[index].data(), dec_a[k][j], b_rlwe.data(), N,
                bigMod, 1);
            intel::hexl::EltwiseFMAMod(a_rlwe.data(), ks[index + 1].data(), dec_a[k][j], a_rlwe.data(), N,
                bigMod, 1);
        }
    }

    return 0;
}

// use mult mod instead of fma mod
int ksKey_hint_variant(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& ks)
{
    lweToRlwe(a);
    std::vector<std::vector<uint64_t> > dec_a(ell, std::vector<uint64_t>(N, 0));
    decompose(dec_a, a, ell, 0x01 << 17, Bg);

    // std::vector<std::vector<uint64_t> > dec_a(ell * N, std::vector<uint64_t>(N));
    // decompose_variant(dec_a, a, ell, 0x01 << 17, Bg);

    std::vector<uint64_t> temp_a(N);
    std::vector<uint64_t> temp_b(N);
    std::vector<uint64_t> temp_dec(N);

    for (size_t j = 0; j < N; j++)
    {
        for (size_t k = 0; k < ell; k++)
        {
            int index = 2 * ell * j + 2 * k;
            
            element_to_vector(temp_dec, dec_a[k][j]);
            intel::hexl::EltwiseMultMod(temp_b.data(), ks[index].data(), temp_dec.data(), N, bigMod, 1);
            intel::hexl::EltwiseMultMod(temp_a.data(), ks[index + 1].data(), temp_dec.data(), N, bigMod, 1);
            for (size_t p = 0; p < N; p++)
            {
                b_rlwe[p] += temp_b[p];
                a_rlwe[p] += temp_a[p];
            }
            // intel::hexl::EltwiseAddMod(b_rlwe.data(), b_rlwe.data(), temp_b.data(), N, bigMod);
            // intel::hexl::EltwiseAddMod(a_rlwe.data(), a_rlwe.data(), temp_a.data(), N, bigMod);
        }
    }

    intel::hexl::EltwiseReduceMod(b_rlwe.data(), b_rlwe.data(), N, bigMod, bigMod, 1);
    intel::hexl::EltwiseReduceMod(a_rlwe.data(), a_rlwe.data(), N, bigMod, bigMod, 1);

    return 0;
}

// the first construction

int ksKey_hint_first_construction(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& pk,
                std::vector<std::vector<uint64_t> >& data,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    intel::hexl::NTT ntts(N, bigMod);
    ntts.ComputeForward(a.data(), a.data(), 1, 1);

    std::vector<std::vector<uint64_t> > temp(N, std::vector<uint64_t>(N, 0));    

    for (size_t i = 0; i < N; i++)
    {
        intel::hexl::EltwiseMultMod(temp[i].data(), a.data(), data[i].data(), N, bigMod, 1);
        
        ntts.ComputeInverse(temp[i].data(), temp[i].data(), 1, 1);
    }
    // TODO: check it
    for (size_t i = 0; i < N; i++)
    {
        lweToRlwe(temp[i]);
    }

    transpose(temp);

    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(N, 0));

    // check_recontruct(dec_a, a, ell, Base, Bg);
    std::vector<uint64_t> t(N);
    for (size_t j = 0; j < N; j++)
    {
        decompose(dec_a, temp[j], ellnum, base, BBg);
        for (size_t k = 0; k < ellnum; k++)
        {
            int index = 2 * ellnum * j + 2 * k;

            // b_rlwe and a_rlwe are in ntt form
            ntts.ComputeForward(dec_a[k].data(), dec_a[k].data(), 1, 1);

            intel::hexl::EltwiseMultMod(t.data(), pk[index].data(), dec_a[k].data(), N,
                bigMod, 1);
            intel::hexl::EltwiseAddMod(b_rlwe.data(), b_rlwe.data(), t.data(), N, bigMod);

            intel::hexl::EltwiseMultMod(t.data(), pk[index + 1].data(), dec_a[k].data(), N,
                bigMod, 1);
            intel::hexl::EltwiseAddMod(a_rlwe.data(), a_rlwe.data(), t.data(), N, bigMod);
        }
    }
    ntts.ComputeInverse(b_rlwe.data(), b_rlwe.data(), 1, 1);
    ntts.ComputeInverse(a_rlwe.data(), a_rlwe.data(), 1, 1);

    return 0;
}

/*
int ksKey_hint(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
                std::vector<uint64_t>& a, std::vector<std::vector<uint64_t> >& ks)
{
    std::vector<std::vector<uint64_t> > dec_a(ell, std::vector<uint64_t>(N, 0));

    decompose(dec_a, a, ell, 0x01 << 17, Bg);

    std::vector<uint64_t> temp_a(N);
    std::vector<uint64_t> temp_b(N);
    std::cout << "break 3" << std::endl;

    for (size_t j = 0; j < N; j++)
    {
        for (size_t k = 0; k < ell; k++)
        {
            int index = 2 * ell * j + 2 * k;
            intel::hexl::EltwiseFMAMod(temp_b.data(), ks[index].data(), dec_a[k][j], nullptr, N,
                bigMod, 1);
            intel::hexl::EltwiseFMAMod(temp_a.data(), ks[index + 1].data(), dec_a[k][j], nullptr, N,
                bigMod, 1);

            intel::hexl::EltwiseAddMod(b_rlwe.data(), b_rlwe.data(), temp_b.data(), N, bigMod);
            intel::hexl::EltwiseAddMod(b_rlwe.data(), b_rlwe.data(), temp_a.data(), N, bigMod);
        }
    }
    std::cout << "break 4" << std::endl;
    intel::hexl::EltwiseReduceMod(b_rlwe.data(), b_rlwe.data(), N, bigMod, 1, 1);
    intel::hexl::EltwiseReduceMod(a_rlwe.data(), a_rlwe.data(), N, bigMod, 1, 1);

    std::cout << "break 5" << std::endl;

    return 0;
}
*/
#endif



#ifndef INTEL_HEXL
/**s
 * @brief 
 * 
 * @param b_rlwe the result ciphertext, part b
 * @param a_rlwe the result ciphertext, part a
 * @param a the query(first), part a
 * @param b the query(first), part b
 * @param ks keySwitch key,
 * @param data the database, store in ntt form
 * 
 */
int answer(uint64_t* b_rlwe, uint64_t* a_rlwe, uint64_t* b, uint64_t* a,
           std::vector<std::vector<uint64_t> >& ks, std::vector<std::vector<uint64_t> >& data,
           bool isntt_b, bool isntt_data)
{
    if (!isntt_b | !isntt_data)
    {
        std::cout << "Error: b or data should be ntt form." << std::endl;
        return 0;
    }
    
    // ntt form
    uint64_t f = 11; // 1/n 
    uint64_t* b_vec = new uint64_t[N]; // store pvw lwe, part b
    for (size_t i = 0; i < N; i++)
    {
        uint64_t temp = 0;
        for (size_t j = 0; j < N; j++)
        {
            /**
             * @brief The each element of first line in the DFT^-1 matrix is 1.
             *        Thefore, we sum all elements and multi it with 1/n.
             *        note: the added temp should be smaller than 2^64.
             */
            temp += barret_reduce((__int128_t)b[i] * data[i][j]); // the i col of the database
        }
        b_vec[i]  = barret_reduce((__int128_t)temp * f);
    }

#ifndef DATA_INDEPENDENT_HANDLE
    ksKey_hint(b_rlwe, a_rlwe, a, ks);
#endif

    // 
    for (size_t i = 0; i < N; i++)
    {
        // TODO: mod subtraction
        b_rlwe[i] -= b_vec[i];
    }

    delete[] b_vec;
    // TODO: add another multiplication

    return 1;
}
#endif


/**
 * @brief in fact, these is nothing doing for a_rlwe and query_a
 * 
 */
#ifdef INTEL_HEXL
int answer_first_dimension(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
            std::vector<uint64_t>& query_b,// std::vector<uint64_t>& query_a,
            std::vector<std::vector<uint64_t> >& data,// std::vector<std::vector<uint64_t> >& ks,
            bool isntt_b, bool isntt_data)
{
    if (!isntt_b | !isntt_data)
    {
        std::cout << "Error: b or data should be ntt form." << std::endl;
        return 0;
    }

    // uint64_t* b_vec = new uint64_t[N]; // store pvw lwe, part b
    std::vector<uint64_t> b_vec(N);
    intel::hexl::EltwiseSubMod(a_rlwe.data(), b_vec.data(), a_rlwe.data(), N, bigMod); // -a_rlwe

    std::vector<uint64_t> t_i(N);

    for (size_t i = 0; i < N; i++)
    {

        uint64_t temp = 0;
        intel::hexl::EltwiseMultMod(t_i.data(), query_b.data(), data[i].data(), N,
                              bigMod, 1);
#ifndef MISSING_ADD
        for (size_t j = 0; j < N; j++)
        {
            //
            // * @brief The each element of first line in the DFT^-1 matrix is 1.
            // *        Therefore, we sum all elements and multiply it with 1/n.
            // *        note: the added temp should be smaller than 2^64.
            //
            temp += t_i[j];
        }
        b_vec[i] = temp;
#endif
    }

    // TODO: check if it's requried to a Reduce function for b_vec
    intel::hexl::EltwiseReduceMod(b_vec.data(), b_vec.data(), N, bigMod, bigMod, 1);
    intel::hexl::EltwiseFMAMod(b_vec.data(), b_vec.data(), bNinv, nullptr, N, bigMod, 1);

    intel::hexl::EltwiseSubMod(b_rlwe.data(), b_vec.data(), b_rlwe.data(), N, bigMod);

    // delete[] b_vec;
    // TODO: add another multiplication

    // auto stop = std::chrono::high_resolution_clock::now();

    // auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "hexl answer costs " << glapsed.count() << " us." << std::endl;

    return 1;
}
#endif

#ifdef INTEL_HEXL
int answer_two_dimensions(RlweCiphertext& result,
            std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe,
            std::vector<uint64_t>& b,
            std::vector<std::vector<uint64_t> >& data,
            RGSWCiphertext& query,
            bool isntt_b, bool isntt_data)
{
    if (!isntt_b | !isntt_data)
    {
        std::cout << "Error: b or data should be ntt form." << std::endl;
        return 0;
    }

    // store pvw lwe, part b
    std::vector<uint64_t> b_vec(N);
    intel::hexl::EltwiseSubMod(a_rlwe.data(), b_vec.data(), a_rlwe.data(), N, bigMod); // -a_rlwe

    std::vector<uint64_t> t_i(N);

    for (size_t i = 0; i < N; i++)
    {
        uint64_t temp = 0;
        intel::hexl::EltwiseMultMod(t_i.data(), b.data(), data[i].data(), b.size(),
                              bigMod, 1);
        for (size_t j = 0; j < N; j++)
        {
            temp += t_i[j];
        }
        b_vec[i] = temp;
    }

    intel::hexl::EltwiseReduceMod(b_vec.data(), b_vec.data(), N, bigMod, bigMod, 1);
    intel::hexl::EltwiseFMAMod(b_vec.data(), b_vec.data(), bNinv, nullptr, N, bigMod, 1);
    intel::hexl::EltwiseSubMod(b_rlwe.data(), b_vec.data(), b_rlwe.data(), N, bigMod);

    // the second dimension
    RlweCiphertext firstDim(N, bigMod);
    firstDim.b = b_rlwe;
    firstDim.a = a_rlwe;
    firstDim.setIsNtt(false);
    externalProduct(result, firstDim, query);

    return 1;
}
#endif

#ifdef INTEL_HEXL
int answer_packing(std::vector<uint64_t>& b_rlwe, std::vector<uint64_t>& a_rlwe, std::vector<uint64_t>& b, std::vector<uint64_t>& a,
                std::vector<std::vector<uint64_t> >& data, std::vector<std::vector<uint64_t> >& ks,
                bool isntt_b, bool isntt_data)
{
    if (!isntt_b | !isntt_data)
    {
        std::cout << "Error: b or data should be ntt form." << std::endl;
        return 0;
    }

    uint64_t f = 114521; // 1/n 
    uint64_t* b_vec = new uint64_t[N]; // store pvw lwe, part b

    std::vector<uint64_t> t_i(N);

    for (size_t i = 0; i < N; i++)
    {

    // auto start = std::chrono::high_resolution_clock::now();
        uint64_t temp = 0;
        intel::hexl::EltwiseMultMod(t_i.data(), b.data(), data[i].data(), b.size(),
                              bigMod, 1);
#ifndef MISSING_ADD
        for (size_t j = 0; j < N; j++)
        {
            temp += t_i[j];
        }
        b_vec[i]  = barret_reduce((__int128_t)temp * f);
#endif
    }

    // 
    for (size_t i = 0; i < N; i++)
    {
        // TODO: mod subtraction
        b_rlwe[i] -= b_vec[i];
    }

    delete[] b_vec;

    

    return 1;
}
#endif

/**
 * @brief 
 * 
 * @param pk
 * @param queryKey 
 * @param answerKey
 * @param data data[i] is always ntt form.
 */
void pkKeyGen(std::vector<std::vector<uint64_t> >& pk,
                 Secret& queryKey, Secret& answerKey,
                 int32_t ellnum, uint64_t base, uint64_t BBg)
{
        // queryKey.toCoeffForm();
#ifdef INTEL_HEXL
        std::vector<std::vector<uint64_t> > mess(N, std::vector<uint64_t>(N, 0));

        uint64_t length = queryKey.getLength();
        uint64_t modulus = queryKey.getModulus();

        intel::hexl::NTT ntts = queryKey.getNTT();

        if (queryKey.isNttForm())
        {
            queryKey.toCoeffForm();
        }

        std::vector<uint64_t> mult_factor = powerOfBg(base, BBg, ellnum);

        // encrypt power of message
        std::vector<uint64_t> temp(N);
        for (size_t i = 0; i < N; i++)
        {        
            std::vector<uint64_t> mess(N);
            mess[0] = queryKey.getData(i);
            for (size_t j = 0; j < ellnum; j++)
            {
                intel::hexl::EltwiseFMAMod(temp.data(), mess.data(),
                                    mult_factor[j],
                                    nullptr, length, modulus, 1);
                int index = 2 * ellnum * i + 2 * j;

                //TODO: add another encrypt to encrypt cofficient, and this would be faster
                ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
                encrypt(pk[index].data(), pk[index + 1].data(), answerKey, temp.data());
                
                // return to coeffcient form
                // ntts.ComputeInverse(pk[index].data(), pk[index].data(), 1, 1);
                // ntts.ComputeInverse(pk[index + 1].data(), pk[index + 1].data(), 1, 1);
            }
        }
#endif
}

/**
 * @brief 
 * 
 * @param ks 
 * @param queryKey 
 * @param answerKey 
 * @param data data[i] is always ntt form.
 */
void dummy_ksKeyGen(std::vector<std::vector<uint64_t> >& ks,
                 Secret& queryKey, Secret& answerKey,
                 std::vector<std::vector<uint64_t> >& data)
{ 
        // queryKey.toCoeffForm();
#ifdef INTEL_HEXL
        std::vector<std::vector<uint64_t> > mess(N, std::vector<uint64_t>(N, 0));

        uint64_t length = queryKey.getLength();
        uint64_t modulus = queryKey.getModulus();

        intel::hexl::NTT ntts = queryKey.getNTT();
        for (size_t i = 0; i < N; i++)
        {
            intel::hexl::EltwiseMultMod(mess[i].data(), data[i].data(), queryKey.getData().data(),
                              queryKey.getLength(), queryKey.getModulus(), 1);

            ntts.ComputeInverse(mess[i].data(), mess[i].data(), 1, 1);
        }
        
        transpose(mess);

        std::vector<uint64_t> mult_factor = powerOfBg(Base, Bg, ell);

        // encrypt power of message
        std::vector<uint64_t> temp(N);
        for (size_t i = 0; i < N; i++)
        {        
            for (size_t j = 0; j < ell; j++)
            {
                intel::hexl::EltwiseFMAMod(temp.data(), mess[i].data(),
                                    mult_factor[j],
                                    nullptr, length, modulus, 1);
                int index = 2 * ell * i + 2 * j;

                //TODO: add another encrypt to encrypt cofficient, and this would be faster
                ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
                encrypt(ks[index].data(), ks[index + 1].data(), answerKey, temp.data());
                
                // return to coeffcient form
                ntts.ComputeInverse(ks[index].data(), ks[index].data(), 1, 1);
                ntts.ComputeInverse(ks[index + 1].data(), ks[index + 1].data(), 1, 1);
            }
        }
#endif
}

/**
 * @brief inline. evaluate automorphic transform
 * 
 * @param result always store in coeffiecints form
 * @param index 
 * @param autokey store in ntt form
 */
void evalAuto(RlweCiphertext& result, const int32_t index, const AutoKey& autokey)
{
    uint64_t length = result.getLength();
    uint64_t modulus = result.getModulus();
    int32_t ellnum = autokey.getEllnum();

    if (result.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = autokey.getNTT();
        ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
        ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
        result.setIsNtt(false);
#endif
    }

    std::vector<uint64_t> temp_a(length), temp_b(length);

    // (temp_a, temp_b) = automorphic^()(a, b)
    for (size_t i = 0; i < length; i++)
    {
#if N == 2048
        uint64_t destination = (i * index) & 0x0fff; // mod 2N
#elif N == 256
        uint64_t destination = (i * index) & 0x01ff; // mod 2N
#elif N == 4096
        uint64_t destination = (i * index) & 0x1fff; // mod 2N
#else
        uint64_t destination = (i * index) % (2 * length);
#endif
        if (destination >= length)
        {
            temp_a[destination - length] = modulus - result.a[i];
            temp_b[destination - length] = modulus - result.b[i];
        } else {
            temp_a[destination] = result.a[i];
            temp_b[destination] = result.b[i];
        }
        result.a[i] = 0;
    }

    // keyswitch
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));
    decompose(dec_a, temp_a, ellnum, autokey.getBase(), autokey.getBg(), modulus);

    std::vector<RlweCiphertext> autokey_index = autokey.keyMap.at(index);

    std::vector<uint64_t> temp(length, 0);
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = autokey.getNTT();
    ntts.ComputeForward(result.b.data(), temp_b.data(), 1, 1);

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(dec_a[i].data(), dec_a[i].data(), 1, 1);

        // (0, b) - g^-1(a) * key
        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i].a.data(), length, modulus, 1);

        intel::hexl::EltwiseSubMod(result.a.data(), result.a.data(), temp.data(), length, modulus);


        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i].b.data(), length, modulus, 1);

        intel::hexl::EltwiseSubMod(result.b.data(), result.b.data(), temp.data(), length, modulus);
    }

    ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
    ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
#endif
}

/**
 * @brief evaluate trace function
 * 
 * @param result 
 * @param logNn 
 * @param autokey 
 */
void evalTrNn(RlweCiphertext& result, int32_t log2Nn, const AutoKey& autokey)
{
    uint64_t length = result.getLength();
    uint64_t modulus = result.getModulus();

    std::vector<uint64_t> temp_a(length);
    std::vector<uint64_t> temp_b(length);
    
    for (size_t k = 1; k <= log2Nn; k++)
    {
    // auto start = std::chrono::high_resolution_clock::now();
        int32_t index = (result.getLength() >> k << 1) + 1;
        copy(result.a.begin(), result.a.end(), temp_a.begin());
        copy(result.b.begin(), result.b.end(), temp_b.begin());

        evalAuto(result, index, autokey);

#ifdef INTEL_HEXL
        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp_a.data(),
                            length, modulus);
        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp_b.data(),
                            length, modulus);
#endif

    // auto stop = std::chrono::high_resolution_clock::now();

    // auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "evalAuto costs " << glapsed.count() << " us." << std::endl;
    }
}

/**
 * @brief the core part of packingLWEs algorithm.
 *  Refer to Algorithm 2 in [https://eprint.iacr.org/2020/015].
 * @param result
 * @param cipher1  the input RLWE ciphertext 2
 * @param curl 
 */
void convertTwoLWEs(RlweCiphertext& result, RlweCiphertext& cipher,
                    const uint64_t twol, const AutoKey& autokey)
{   
    int32_t XNdiv2l = N / twol;
    
    cipher.multConst(XNdiv2l);

    RlweCiphertext temp = result;

    result.subAndEqual(cipher);
    cipher.addAndEqual(temp);

    int32_t index = twol + 1;
    evalAuto(result, index, autokey);

    result.addAndEqual(cipher);
}

/**
 * @brief 
 * 
 * @param result response of answer
 * @param lwes the requried packing lwe ciphertexts (note: store them in RlweCiphertext).
 *      The size should be power of two. lwes should be store in coefficient form.
 */
void packingLWEs(RlweCiphertext& result, std::vector<RlweCiphertext>& lwes,
                 const AutoKey& autokey)
{
    // auto start = std::chrono::high_resolution_clock::now();
    
    uint64_t lwesNum = lwes.size(); // power of two

    // the input lwes should be in coeffcient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = autokey.getNTT();
    for (size_t i = 0; i < lwesNum; i++)
    {
    if (lwes[i].getIsNtt())
        {
            ntts.ComputeInverse(lwes[i].a.data(), lwes[i].a.data(), 1, 1);
            ntts.ComputeInverse(lwes[i].b.data(), lwes[i].b.data(), 1, 1);
            lwes[i].setIsNtt(false);
#endif
        }
    }

    for (size_t i = lwesNum/2; i > 0; i >>= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            convertTwoLWEs(lwes[j], lwes[i + j], lwesNum / i, autokey);
        }
    }
    /*
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "evalTrNn costs " << glapsed.count() << " us." << std::endl;
    */
    result = lwes[0];

    /**
     *  lazy variant in [https://eprint.iacr.org/2020/015].
     *  Because we don't have to let unwanted coefficients to be zeros. Therefore, we let
     *  evalTrNn function invisible.
     **/

    // evalTrNn(result, log2(result.getLength() / lwesNum), autokey);
}

/*******************************************************************/
/*********************** ONE-TIME SETUP ****************************/

/**
 * @brief note: in the one-time setup phase, the server should homomorphically evaluate to
 * obtain ks. In other words, is to replace dummy_ksKeyGen function.
 *
 */

void evalAutoRNS(RlweCiphertext& result1, RlweCiphertext& result2,
                    const int32_t index, const AutoKeyRNS& autokey)
{
    uint64_t length = result1.getLength();
    assert(length == result2.getLength());

    uint64_t modulus1 = result1.getModulus();
    uint64_t modulus2 = result2.getModulus();
    int32_t ellnum = autokey.getEllnum();

    if (result1.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts1 = autokey.getNTT1();
        ntts1.ComputeInverse(result1.a.data(), result1.a.data(), 1, 1);
        ntts1.ComputeInverse(result1.b.data(), result1.b.data(), 1, 1);
        result1.setIsNtt(false);
#endif
    }

    if (result2.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts2 = autokey.getNTT2();
        ntts2.ComputeInverse(result2.a.data(), result2.a.data(), 1, 1);
        ntts2.ComputeInverse(result2.b.data(), result2.b.data(), 1, 1);
        result2.setIsNtt(false);
#endif
    }

    std::vector<uint64_t> temp_a1(length), temp_b1(length), temp_a2(length), temp_b2(length);

    // (temp_a, temp_b) = automorphic^()(a, b)
    for (size_t i = 0; i < length; i++)
    {
#if N == 2048
        uint64_t destination = (i * index) & 0x0fff; // mod 2N
#elif N == 256
        uint64_t destination = (i * index) & 0x01ff; // mod 2N
#elif N == 4096
        uint64_t destination = (i * index) & 0x1fff; // mod 2N
#else
        uint64_t destination = (i * index) % (2 * length);
#endif
        if (destination >= length)
        {
            temp_a1[destination - length] = modulus1 - result1.a[i];
            temp_b1[destination - length] = modulus1 - result1.b[i];
            temp_a2[destination - length] = modulus2 - result2.a[i];
            temp_b2[destination - length] = modulus2 - result2.b[i];
        } else {
            temp_a1[destination] = result1.a[i];
            temp_b1[destination] = result1.b[i];
            temp_a2[destination] = result2.a[i];
            temp_b2[destination] = result2.b[i];
        }
        result1.a[i] = 0;
        result2.a[i] = 0;
    }

    // keyswitch
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));
    std::vector<std::vector<uint64_t> > dec_a2(ellnum, std::vector<uint64_t>(length, 0));
    decompose_crt(dec_a, dec_a2, temp_a1, temp_a2, ellnum, autokey.getBase(), autokey.getBg());

    std::vector<RlweCiphertext> autokey_index = autokey.keyMap.at(index);

    std::vector<uint64_t> temp(length, 0);
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = autokey.getNTT1();
    ntts1.ComputeForward(result1.b.data(), temp_b1.data(), 1, 1);

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts1.ComputeForward(dec_a[i].data(), dec_a[i].data(), 1, 1);

        // (0, b) - g^-1(a) * key
        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i].a.data(), length, modulus1, 1);

        intel::hexl::EltwiseSubMod(result1.a.data(), result1.a.data(), temp.data(), length, modulus1);


        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i].b.data(), length, modulus1, 1);

        intel::hexl::EltwiseSubMod(result1.b.data(), result1.b.data(), temp.data(), length, modulus1);
    }

    ntts1.ComputeInverse(result1.a.data(), result1.a.data(), 1, 1);
    ntts1.ComputeInverse(result1.b.data(), result1.b.data(), 1, 1);
#endif

    // TOOD: check another crt part
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts2 = autokey.getNTT2();
    ntts2.ComputeForward(result2.b.data(), temp_b2.data(), 1, 1);

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts2.ComputeForward(dec_a[i].data(), dec_a2[i].data(), 1, 1);

        // (0, b) - g^-1(a) * key
        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i + ellnum].a.data(), length, modulus2, 1);

        intel::hexl::EltwiseSubMod(result2.a.data(), result2.a.data(), temp.data(), length, modulus2);


        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    autokey_index[i + ellnum].b.data(), length, modulus2, 1);

        intel::hexl::EltwiseSubMod(result2.b.data(), result2.b.data(), temp.data(), length, modulus2);
    }

    ntts2.ComputeInverse(result2.a.data(), result2.a.data(), 1, 1);
    ntts2.ComputeInverse(result2.b.data(), result2.b.data(), 1, 1);
#endif
}

/**
 * @brief evaluate expand function, rns varian
 * 
 * @param result size is 2 * N 
 * @param input size is 2
 * @param autokey 
 */
void evalExpandRNS(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
                    const AutoKeyRNS& autokey)
{
    // the input ciphertexts should be coefficient form
    if (input[0].getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts1 = autokey.getNTT1();
        ntts1.ComputeInverse(input[0].a.data(), input[0].a.data(), 1, 1);
        ntts1.ComputeInverse(input[0].b.data(), input[0].b.data(), 1, 1);
        input[0].setIsNtt(false);
#endif
    }

    if (input[1].getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts2 = autokey.getNTT2();
        ntts2.ComputeInverse(input[1].a.data(), input[1].a.data(), 1, 1);
        ntts2.ComputeInverse(input[1].b.data(), input[1].b.data(), 1, 1);
        input[1].setIsNtt(false);
#endif
    }

    int32_t l = log2(N);

    // copy_value(output[0], input);
    result[0] = input[0];
    result[N] = input[1];
    RlweCiphertext temp1, temp2;

    for (size_t j = 0; j < l; j++)
    {
        int32_t index = (0x01 << (l - j)) + 1;
        for (size_t k = 0; k < (0x01 << j); k++)
        {
            // mult_factor(temp1, output[k], 0x01 << j, length, modulus);

            // automorphic(temp2, output[k], index, modulus);
            // add(output[k], output[k], temp2, modulus);
            
            // automorphic(temp2, temp1, index, modulus);
            // add(output[k + (0x01 << j)], temp1, temp2, modulus);

            temp1 = result[k];
            temp2 = result[k + N];

            evalAutoRNS(result[k], result[k + N], index, autokey);
            addRnsCiphertext(result[k], result[k + N], temp1, temp2);

            temp1.multNegConst(0x01 << j);
            temp2.multNegConst(0x01 << j);

            int32_t second_index = k + (0x01 << j);
            evalAutoRNS(temp1, temp2, index, autokey);
            addRnsCiphertext(result[second_index], result[second_index + N],
                                temp1, temp2);
        }
    }

    // N^-1 * result
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            result[i].a[j] = static_cast<uint128_t>(result[i].a[j]) * bNinv % bigMod;
            result[i].b[j] = static_cast<uint128_t>(result[i].b[j]) * bNinv % bigMod;

            result[i + N].a[j] = static_cast<uint128_t>(result[i + N].a[j]) * bNinv2 % bigMod2;
            result[i + N].b[j] = static_cast<uint128_t>(result[i + N].b[j]) * bNinv2 % bigMod2;
        }
    }
}

// TODO: check or remove it
/**
 * @brief evaluate trace function, rns variant
 * 
 * @param result 
 * @param logNn 
 * @param autokey 
 */
void evalTrNnRNS(std::vector<RlweCiphertext>& result, int32_t log2Nn, const AutoKeyRNS& autokey)
{
    uint64_t length = result[0].getLength();
    uint64_t modulus1 = result[0].getModulus();
    uint64_t modulus2 = result[1].getModulus();

    std::vector<uint64_t> temp_a1(length);
    std::vector<uint64_t> temp_b1(length);

    std::vector<uint64_t> temp_a2(length);
    std::vector<uint64_t> temp_b2(length);
    
    for (size_t k = 1; k <= log2Nn; k++)
    {
        int32_t index = (length >> k << 1) + 1;
        copy(result[0].a.begin(), result[0].a.end(), temp_a1.begin());
        copy(result[0].b.begin(), result[0].b.end(), temp_b1.begin());
        copy(result[1].a.begin(), result[1].a.end(), temp_a2.begin());
        copy(result[1].b.begin(), result[1].b.end(), temp_b2.begin());

        evalAutoRNS(result[0],result[1], index, autokey);

#ifdef INTEL_HEXL
        intel::hexl::EltwiseAddMod(result[0].a.data(), result[0].a.data(), temp_a1.data(),
                            length, modulus1);
        intel::hexl::EltwiseAddMod(result[0].b.data(), result[0].b.data(), temp_b1.data(),
                            length, modulus1);

        intel::hexl::EltwiseAddMod(result[1].a.data(), result[1].a.data(), temp_a2.data(),
                            length, modulus2);
        intel::hexl::EltwiseAddMod(result[1].b.data(), result[1].b.data(), temp_b2.data(),
                            length, modulus2);

#endif

    }
}

void publicKeyGen(std::vector<std::vector<uint64_t> >& publicKey,
                Secret& queryKey, Secret& answerKey)
{
    uint64_t length = answerKey.getLength();
    uint64_t modulus1 = answerKey.getModulus();
    uint64_t modulus2 = bigMod2;

    if (queryKey.isNttForm())
    {
        queryKey.toCoeffForm();
    }

    if (answerKey.isNttForm())
    {
        answerKey.toCoeffForm();
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = answerKey.getNTT();
    intel::hexl::NTT ntts2(length, modulus2);

    std::vector<uint64_t> secret1(length);
    std::vector<uint64_t> secret2(length);

    for (size_t i = 0; i < answerKey.getLength(); i++)
    {
        secret1[i] = answerKey.getData(i);
        secret2[i] = answerKey.getData(i);
    }
    guass_to_modulus(secret2.data(), modulus1, modulus2);
    ntts1.ComputeForward(secret1.data(), secret1.data(), 1, 1);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    std::vector<uint64_t> message(length);
    std::vector<uint64_t> temp(length);
    std::vector<uint64_t> mult_factor = powerOfBg(Base, Bg, ell);

    uint64_t crt_interval = 2 * ell * N;

    for (size_t i = 0; i < N; i++)
    {
        message[0] = queryKey.getData(i);
        for (size_t j = 0; j < ell; j++)
        {
            intel::hexl::EltwiseFMAMod(temp.data(), message.data(),
                                    mult_factor[j],
                                    nullptr, length, modulus1, 1);
            // encode for q2 * temp
            encode_crt(temp);

            int32_t index = 2 * ell * i + 2 * j;

            // CRT encryption
            encrypt(publicKey[index].data(), publicKey[index + 1].data(),
                    publicKey[index + crt_interval].data(), publicKey[index + 1 + crt_interval].data(),
                    answerKey, temp.data(),
                    secret1.data(), secret2.data(), // should be ntt form
                    ntts1, ntts2, modulus1, modulus2);
        }
    }

#endif
}

/**
 * @brief 
 * 
 * @param pk the client send pk to the server. dim1: 2*ell*N * 2, dim2: N
 * @param input1 the input rlwe ciphertexts, store in coefficient form. dim: 2*ell
 */
void expand(std::vector<std::vector<uint64_t> >& pk, std::vector<RlweCiphertext>& input,
            AutoKey& autokey)
{
    uint64_t length = input[0].getLength();
    // uint64_t modulus = input[0].getModulus();
    
    for (size_t i = 0; i < 2 * ell; i++)
    {
        lweToRlwe(input[i].a);

        RlweCiphertext temp(input[i].a, input[i].b);

        evalTrNn(temp, length, autokey);
    }

}

/**
 * @brief 
 * 
 * @param ks dim1: 2*ell*N * 2(CRT), dim2: N
 * @param pk the client send pk to the server. dim1: 2*ell*N * 2, dim2: N
 * @param data 
 * @param s_index RLWE(s[s_index])
 */
void evalSingleKsKeyGen(std::vector<std::vector<uint64_t> >& ks,
                 std::vector<std::vector<uint64_t> >& pk,
                 std::vector<std::vector<uint64_t> >& data,
                 int32_t s_index)
{
#ifdef INTEL_HEXL
    uint64_t modulus1 = bigMod;
    uint64_t modulus2 = bigMod2;

    /**
     * homomorphically evaluate native multiplication.
     *            s0       s1        s2         s3          s(n-1)
     * st0[0]    t0[0]   -t0[n-1]  -t0[n-2]  -t0[n-3]  ...  -t0[1]
     * st0[1]    t0[1]    t0[0]    -t0[n-1]  -t0[n-2]  .... 
     * st0[2]    t0[2]    t0[1]     t0[0]    -t0[n-1]  ....
     * st0[3]    t0[3]    t0[2]     t0[1]     t0[0]    ....
     * ...
     * st0[n-1]  t0[n-1]  t0[n-2]   t0[n-3]   t0[n-4]  ...   t0[0]
    */
    std::vector<uint64_t> temp(N);

    bool reverse = false;
    int32_t data_index = 0;

    uint64_t crt_interval = 2 * ell * N;
    for (size_t j = 0; j < N; j++)
    {
        compute_indicator(data_index, reverse, j, s_index);
        for (size_t k = 0; k < ell; k++)
        {
            int32_t pk_index = 2 * ell * j + 2 * k;
            int32_t output_index = 2 * ell * s_index + 2 * k;

            intel::hexl::EltwiseMultMod(temp.data(), pk[pk_index].data(), data[data_index].data(),
                            N, modulus1, 1);
            if (!reverse)
            intel::hexl::EltwiseAddMod(ks[output_index].data(), ks[output_index].data(), temp.data(),
                            N, modulus1);
            else
            intel::hexl::EltwiseSubMod(ks[output_index].data(), ks[output_index].data(), temp.data(),
                            N, modulus1);

            intel::hexl::EltwiseMultMod(temp.data(), pk[pk_index + 1].data(), data[data_index].data(),
                            N, modulus1, 1);
            if (!reverse)
            intel::hexl::EltwiseAddMod(ks[output_index + 1].data(), ks[output_index + 1].data(), temp.data(),
                            N, modulus1);
            else
            intel::hexl::EltwiseSubMod(ks[output_index + 1].data(), ks[output_index + 1].data(), temp.data(),
                            N, modulus1);

            // handle CRT 2
            intel::hexl::EltwiseMultMod(temp.data(), pk[pk_index + crt_interval].data(), data[data_index + N].data(),
                            N, modulus2, 1);
            if (!reverse)
            intel::hexl::EltwiseAddMod(ks[output_index + crt_interval].data(), ks[output_index + crt_interval].data(),
                            temp.data(), N, modulus2);
            else
            intel::hexl::EltwiseSubMod(ks[output_index + crt_interval].data(), ks[output_index + crt_interval].data(),
                            temp.data(), N, modulus2);

            intel::hexl::EltwiseMultMod(temp.data(), pk[pk_index + 1 + crt_interval].data(), data[data_index + N].data(),
                            N, modulus2, 1);
            if (!reverse)
            intel::hexl::EltwiseAddMod(ks[output_index + 1 + crt_interval].data(), ks[output_index + 1 + crt_interval].data(),
                            temp.data(), N, modulus2);
            else
            intel::hexl::EltwiseSubMod(ks[output_index + 1 + crt_interval].data(), ks[output_index + 1 + crt_interval].data(),
                            temp.data(), N, modulus2);
        }
    }
#endif
}

void evalKsKeyGenForThreads(std::vector<std::vector<uint64_t> >& ks,
                 std::vector<std::vector<uint64_t> >& pk,
                 std::vector<std::vector<uint64_t> >& data,
                 size_t first, size_t last)
{
    for (size_t i = first; i < last; i++)
    {
        evalSingleKsKeyGen(ks, pk, data, i);
    }
}


void evalKsKeyGen(std::vector<std::vector<uint64_t> >& ks,
                 std::vector<std::vector<uint64_t> >& pk,
                 std::vector<std::vector<uint64_t> >& data)
{
#ifndef THREADS
    for (size_t i = 0; i < N; i++)
    {
        evalSingleKsKeyGen(ks, pk, data, i);
    }
#endif

#ifdef THREADS

    size_t threads_num = THREADS_NUM; //8, 16, 32, power of two
    std::thread threads[threads_num];
    size_t threadPerCount = N / threads_num;

    for (size_t i = 0; i < threads_num; i++)
    {
        size_t first = threadPerCount * i;
        size_t last = threadPerCount * (i + 1);

        threads[i] = std::thread(evalKsKeyGenForThreads, ref(ks), ref(pk),
                            ref(data), first, last);
    }


    for (size_t i = 0; i < threads_num; i++)
    {
        threads[i].join();
    }

#endif

}


/**
 * @brief rns modswitch
 * 
 * @param result dim1: 2 * ell * N, dim2: N. store in ntt form
 * @param ks dim1: 2 * ell * N * 2(CRT), dim2: N. store in ntt form
 * @param 
 */
void modSwitch(std::vector<std::vector<uint64_t> >& result,
                std::vector<std::vector<uint64_t> >& ks)
{
    // rns modswitch:
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    uint64_t length = N;
    uint64_t modulus1 = bigMod;

    uint64_t crt_interval = 2 * ell * N;

    intel::hexl::NTT ntts1(N, bigMod);
    intel::hexl::NTT ntts2(N, bigMod2);

    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    for (size_t i = 0; i < result.size(); i++)
    {

#ifdef INTEL_HEXL
        ntts1.ComputeInverse(temp1.data(), ks[i].data(), 1, 1);
        ntts2.ComputeInverse(temp2.data(), ks[i + crt_interval].data(), 1, 1);

        intel::hexl::EltwiseSubMod(result[i].data(), temp1.data(), temp2.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result[i].data(), result[i].data(), mod2inv,
                        nullptr, length, modulus1, 1);

        ntts1.ComputeForward(result[i].data(), result[i].data(), 1, 1);
#endif
    }
}

/**
 * @brief rns modswitch
 * 
 * @param result store in ntt form
 * @param input1 rlwe ciphertext in modulus1
 * @param input2 rlwe ciphertext in modulus2
 */
void modSwitch(RlweCiphertext& result, RlweCiphertext& input1, RlweCiphertext& input2)
{
    // rns modswitch:
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    uint64_t length = N;
    uint64_t modulus1 = bigMod;

    intel::hexl::NTT ntts1(N, bigMod);
    intel::hexl::NTT ntts2(N, bigMod2);

    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    if (input1.getIsNtt())
    {
        ntts1.ComputeInverse(input1.a.data(), input1.a.data(), 1, 1);
        ntts1.ComputeInverse(input1.b.data(), input1.b.data(), 1, 1);
    }
    if (input2.getIsNtt())
    {
        ntts2.ComputeInverse(input2.a.data(), input2.a.data(), 1, 1);
        ntts2.ComputeInverse(input2.b.data(), input2.b.data(), 1, 1);
    }

#ifdef INTEL_HEXL
        intel::hexl::EltwiseSubMod(result.a.data(), input1.a.data(), input2.a.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result.a.data(), result.a.data(), mod2inv,
                        nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result.a.data(), result.a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result.b.data(), input1.b.data(), input2.b.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result.b.data(), result.b.data(), mod2inv,
                        nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result.b.data(), result.b.data(), 1, 1);

        result.setIsNtt(true);
#endif
}
