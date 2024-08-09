#include <stddef.h>
#include <vector>
#include <chrono>
#include <assert.h>
#include <string.h>

#include "lwe.h"
#include "utils.h"
#include "twosteps.h"
#include "encrypt.h"
#include "crt.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

#include <iostream>

#define FBSGS
// #define TIME_COLLECT

// note database should be its signed variant
void sample_database_bsgs(std::vector<std::vector<uint64_t> >& data)
{
    // the database has N * N entries, and each entry has 16 bits    
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N / 2; j++)
        {
            sample_rej:
                uint64_t r = rand() & 0x0ffff;
            if (r > bsgsp)
            {
                goto sample_rej;
            }
            data[i][j] = r;
        }
    }
}


void database_tobsgsntt(std::vector<std::vector<uint64_t> >& result,
                        std::vector<std::vector<uint64_t> >& data, uint64_t modulus, int32_t N1)
{
    uint64_t length = N;

    // in fact, the number of hexl_ntt_index are conjugated for complementary positions,
    // 0 ~ 4095, 1 ~ 4094, 2 ~ 4093, \cdots
    // conjugated = [4095, 4094, 4093, 4092, \cdots, 2048]

    // database precomputation
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length); 

    // for i = 1, result = [d, a, b, c, f, g, h, e]
    // data = [ 0, a, 0, 0]
    //        [ 0, 0, b, 0]
    //        [ 0, 0, 0, c]
    //        [ d, 0, 0, 0]
    //        -------------
    //        [ e, 0, 0, 0]
    //        [ 0, 0, 0, f]
    //        [ 0, 0, g, 0]
    //        [ 0, h, 0, 0]
    for (size_t i = 0; i < N / 2; i++)
    {
        for (size_t j = 0; j < N / 2; j++)
        {
            // diagonal line
            int row = (j - i + N/2) & (N/2 - 1); // the row is always smaller 
            result[i][j] = data[row][j];
            result[i][N - 1 - j] = data[N - 1 - row][j];
        }
    }

    std::vector<std::vector<int32_t> > permutations(length/2, std::vector<int32_t>(length, 0));
    // compute_permutation_matrix(permutations, length/2, length);
    compute_permutation_matrix(permutations, N1, length);

    // rotate 1 ~ N/2
    std::vector<uint64_t> temp1(length, 0);
    std::vector<uint64_t> temp2(length, 0);
    /**
    for (size_t i = 0; i < N / 2; i++)
    {
        // query permutation for data_ntt
        // copy data_ntt[i] to temp
        copy(result[i].begin(), result[i].end(), temp1.begin());
        // encode data_ntt[i]
        for (size_t j = 0; j < N; j++)
        {
            temp2[query_encode[j]] = temp1[j];
        }

        for (size_t j = 0; j < N; j++)
        {
            result[i][j] =  temp2[permutations[i][j]]; // i is rotation number
        }
    }
    **/
    int32_t N2 = N / 2 / N1;
    for (size_t k = 0; k < N2; k++)
    {
        for (size_t i = 0; i < N1; i++)
        {
            // query permutation for data_ntt
            // copy data_ntt[i] to temp
            int32_t kN1i = k * N1 + i; 
            copy(result[kN1i].begin(), result[kN1i].end(), temp1.begin());
            // encode data_ntt[i]
            for (size_t j = 0; j < N; j++)
            {
                temp2[query_encode[j]] = temp1[j];
            }

            for (size_t j = 0; j < N; j++)
            {
                result[kN1i][j] =  temp2[permutations[i][j]]; // i is rotation number
            }
        }
    }

/** if hexl_ntt_index = [5, 25, 125, \cdots]
    for (size_t i = 0; i < N / 2; i++)
    {
        for (size_t j = 0; j < N / 2; j++)
        {
            // cycle diagonal line
            // result[i][j] = data[j][(i + j) % (N/2)];
            int32_t col = (i + j) & (N/2 - 1);
            result[i][j] = data[j][col];

            result[i][N - 1 - j] = data[N - 1 - j][col];
        }
        // for i = 1, result = [a, b, c, d, e, f, g, h]
        // data = [ 0, a, 0, 0]
        //        [ 0, 0, b, 0]
        //        [ 0, 0, 0, c]
        //        [ d, 0, 0, 0]
        //        -------------
        //        [ e, 0, 0, 0]
        //        [ 0, 0, 0, f]
        //        [ 0, 0, g, 0]
        //        [ 0, h, 0, 0]
    }
**/
#ifdef INTEL_HEXL
    // encode
    intel::hexl::NTT nttp(N, bsgsp);
    for (size_t i = 0; i < N / 2; i++)
    {
        nttp.ComputeInverse(result[i].data(), result[i].data(), 1, 1);
    }

    // to signed
    // the database has N/2 * N entries, and each entry has X bits
    uint64_t sub = modulus - bsgsp;

    for (size_t i = 0; i < N / 2; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (result[i][j] > bsgsp / 2)
                result[i][j] += sub;
        }
    }

    // precompute and store in ntt form
#ifdef USE_CRTMOD
    intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
#else
    intel::hexl::NTT ntts(N, bigMod);
#endif
    for (size_t i = 0; i < N / 2; i++)
    {
        ntts.ComputeForward(result[i].data(), result[i].data(), 1, 1);
    }
#endif
}


void inverse_encode(std::vector<uint64_t>& message)
{
    uint64_t length = N;

    intel::hexl::NTT nttp(length, bsgsp);
    nttp.ComputeForward(message.data(), message.data(), 1, 1);

    std::vector<int32_t> query_decode(length, 0);
    compute_query_decode(query_decode, length);

    std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_decode[j]] = temp[j];
    }
}

void query_bsgs(RlweCiphertext& cipher, Secret& queryKey, uint64_t row)
{
    // uint64_t length = queryKey.getLength();
    uint64_t length = N;
    uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length);

    // two rotations are conjugated
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length);

    // message[row] = 1;
    // message[N - 1 - row] = 1;
    message[query_encode[row]] = 1;
    message[query_encode[N - 1 - row]] = 1;

    // intel::hexl::NTT nttp(length, bsgsp);
    intel::hexl::NTT nttp(N, bsgsp);
    nttp.ComputeInverse(message.data(), message.data(), 1, 1);

    // me \in (-p/2. p/2] * Delta
    // if me \in (p/2, p/2], (me - p) * Delta + bigMod
    /**/
    for (auto & me : message)
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
        // refer to `Revisiting BGV and BFV'
        // VARIANT A
        /**
        if (me > bsgsp / 2)
        {
            me = ((long double)me - bsgsp) * bigMod / bsgsp + bigMod;
        } else {
            me = round((long double)me * bigMod / bsgsp);
        }
        **/
        // VARIANT A
        me = round((long double)me * modulus / bsgsp); // TODO: it's same?
    }
    /**/
    /**
    for (auto & me : message)
    {
        me *= bsgsDelta;
    }
    **/
    /**
    for (size_t i = 0; i < N; i++)
    {
        if (message[i] > bsgsp / 2)
        {
            message[i] *= bsgsDelta;
            message[i] += (bigMod - bsgsp * bsgsDelta);
        } else {
            message[i] *= bsgsDelta;
        }
    }
    **/

    // store a in coefficient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();

    ntts.ComputeForward(message.data(), message.data(), 1, 1);
    encrypt(cipher.b.data(), cipher.a.data(), queryKey, message.data());

    // ntts.ComputeInverse(cipher.a.data(), cipher.a.data(), 1, 1);
    cipher.setIsNtt(true);
#endif
}

void query_bsgs_rns(std::vector<RlweCiphertext>& cipher, Secret& queryKey, uint64_t row)
{
    // uint64_t length = queryKey.getLength();
    uint64_t length = N;
    uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length);

    // two rotations are conjugated
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length);

    // message[row] = 1;
    // message[N - 1 - row] = 1;
    message[query_encode[row]] = 1;
    message[query_encode[N - 1 - row]] = 1;

#ifdef INTEL_HEXL
    intel::hexl::NTT nttp(N, bsgsp);
    nttp.ComputeInverse(message.data(), message.data(), 1, 1);

    // me \in (-p/2. p/2] * Delta
    // if me \in (p/2, p/2], (me - p) * Delta + bigMod
    /**/
    for (auto & me : message)
    {
        // VARIANT A
        me = round((long double)me * modulus / bsgsp); // TODO: it's same?
    }

    // store a in coefficient form

    // intel::hexl::NTT ntts = queryKey.getNTT();

    // ntts.ComputeForward(message.data(), message.data(), 1, 1);
    // encrypt(cipher.b.data(), cipher.a.data(), queryKey, message.data());
    encrypt_rns_bsgs(cipher, queryKey, message);

    // cipher.setIsNtt(true);
#endif
}

void decrypt_bsgs(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret)
{
    uint64_t length = secret.getLength();
    uint64_t modulus = secret.getModulus();

    if (!cipher.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = secret.getNTT();
        ntts.ComputeForward(cipher.a.data(), cipher.a.data(), 1, 1);
        ntts.ComputeForward(cipher.b.data(), cipher.b.data(), 1, 1);
        cipher.setIsNtt(true);
#endif
    }
    
    decrypt(message.data(), cipher.b.data(), cipher.a.data(), secret);

    // me shoule be in (-p/2, p/2)
    /**/
    for (auto & me : message)
    {
        if (me > modulus / 2)
        {
            // ((double)me - bigMod) / bsgsDelta is always less than 0, so -.5 for rounding
            // long double r = ((long double)me - bigMod) / bsgsDelta - .5;
            // me = (int64_t)r + bsgsp;
            // TEST:
            // uint64_t test = round(((long double)me - bigMod) / bsgsDelta) + bsgsp;
            // if (test != (int64_t)r + bsgsp) std::cout << "error: " << me << std::endl; 

            me = round(((long double)me - modulus) / bsgsDelta) + bsgsp;
        } else {
            // long double r2 = (long double)me / bsgsDelta + 0.5;
            // TEST:
            // uint64_t test =  round((long double)me / bsgsDelta);
            // if (test != (uint64_t)r2) std::cout << "error: " << me << std::endl; 

            // me = r2;
            me = round((long double)me / bsgsDelta);
        }
    }
    /**/
    /**
    for (auto & me : message)
    {
        long double r2 = (long double)me / bsgsDelta;
        me = round(r2);
    }
    **/

    intel::hexl::NTT nttp(length, bsgsp);
    nttp.ComputeForward(message.data(), message.data(), 1, 1);

    // query decode
    std::vector<int32_t> query_decode(length, 0);
    compute_query_decode(query_decode, length);

    std::vector<uint64_t> temp(length, 0);
    copy(message.begin(), message.end(), temp.begin());
    for (size_t j = 0; j < N; j++)
    {
        message[query_decode[j]] = temp[j];
    }
}

/**
 * decrypt for total protocol
 */
void decrypt_bsgs_total(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret,
                        const int32_t rlwesNum)
{
    uint64_t modulus = secret.getModulus();

    if (!cipher.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = secret.getNTT();
        ntts.ComputeForward(cipher.a.data(), cipher.a.data(), 1, 1);
        ntts.ComputeForward(cipher.b.data(), cipher.b.data(), 1, 1);
        cipher.setIsNtt(true);
#endif
    }
    
    decrypt(message.data(), cipher.b.data(), cipher.a.data(), secret);

    multConst(message, QInv(rlwesNum, modulus), modulus);

    // me shoule be in (-p/2, p/2)
    /**/
    for (auto & me : message)
    {
        if (me > modulus / 2)
        {
            me = round(((long double)me - modulus) / bsgsDelta) + bsgsp;
        } else {
            me = round((long double)me / bsgsDelta);
        }
    }
}

inline void autoInNtt(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,
                      const std::vector<int32_t>& permutation)
{
    int32_t length = input.size();
    for (size_t i = 0; i < length; i++)
    {
        result[i] = input[permutation[i]];
    }
}

// voud rotate(RlweCiphertext& rotated_cipher, const RlweCiphertext& input, size_t i)
/**
 * @brief inline. evaluate automorphic transform
 * 
 * @param result always store in coeffiecints form
 * @param index 
 * @param autokey store in ntt form
 */
void evalAuto(RlweCiphertext& result, RlweCiphertext& input,
              const int32_t index, const AutoKeyBSGS& autokey,
              const int32_t ellnum, const uint64_t PP, const uint64_t BBg)
{
    // use the second parameters
    uint64_t length = input.getLength();
    uint64_t modulus = input.getModulus();
    // int32_t ellnum = autokey.getEllnumGS();

    if (input.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = autokey.getNTT();
        ntts.ComputeInverse(input.a.data(), input.a.data(), 1, 1);
        ntts.ComputeInverse(input.b.data(), input.b.data(), 1, 1);
        input.setIsNtt(false);
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
            temp_a[destination - length] = modulus - input.a[i];
            temp_b[destination - length] = modulus - input.b[i];
        } else {
            temp_a[destination] = input.a[i];
            temp_b[destination] = input.b[i];
        }
        result.a[i] = 0;
    }

    // keyswitch
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));
    decompose(dec_a, temp_a, ellnum, PP, BBg, modulus);

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

    // return in ntt form
    result.setIsNtt(true);
    // ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
    // ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
#endif
}

void evalAutoCRT(RlweCiphertext& result, std::vector<std::vector<uint64_t> >& input,
              const int32_t N2, const AutoKeyBSGS& autokey,
              const int32_t ellnum, const uint64_t PP, const uint64_t BBg)
{
    // use the second parameters
    size_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus();
    // int32_t ellnum = autokey.getEllnumGS();

    int32_t N1 = N / 2 / N2;

#ifdef INTEL_HEXL

    intel::hexl::NTT ntts = autokey.getNTT();
    std::vector<uint64_t> temp_a(length), temp_b(length);

for (size_t cipher_num = 1; cipher_num < N2; cipher_num++)
{
    ntts.ComputeInverse(input[cipher_num].data(), input[cipher_num].data(), 1, 1);
    ntts.ComputeInverse(input[cipher_num].data() + length, input[cipher_num].data() + length, 1, 1);

    size_t index = pow_mod(5, N1 * cipher_num, 2 * N);
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));

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
            temp_a[destination - length] = modulus - input[cipher_num][i];
            temp_b[destination - length] = modulus - input[cipher_num][i + length];
        } else {
            temp_a[destination] = input[cipher_num][i];
            temp_b[destination] = input[cipher_num][i + length];
        }
        // result.a[i] = 0;
    }

    // keyswitch

    decompose(dec_a, temp_a, ellnum, PP, BBg, modulus);

    std::vector<RlweCiphertext> autokey_index = autokey.keyMap.at(index);

    std::vector<uint64_t> temp(length, 0);

    ntts.ComputeForward(temp_b.data(), temp_b.data(), 1, 1);    
    intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp_b.data(), length, modulus);

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
}
    // return in ntt form
    // result.setIsNtt(true);
#endif
}

//  evalAutoRNSCRT is as same as evalAutoCRT other that the input is AutoKeyBSGSRNS here
//  TODO: is slower than evalAutoCRT?
void evalAutoRNSCRT(RlweCiphertext& result, std::vector<std::vector<uint64_t> >& input,
              const int32_t N2, const AutoKeyBSGSRNS& autokey,
              const int32_t ellnum, const uint64_t PP, const uint64_t BBg)
{
    // use the second parameters
    size_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus();
    // int32_t ellnum = autokey.getEllnumGS();

    int32_t N1 = N / 2 / N2;

#ifdef INTEL_HEXL

    intel::hexl::NTT ntts = autokey.getNTT();
    std::vector<uint64_t> temp_a(length), temp_b(length);

for (size_t cipher_num = 1; cipher_num < N2; cipher_num++)
{
    ntts.ComputeInverse(input[cipher_num].data(), input[cipher_num].data(), 1, 1);
    ntts.ComputeInverse(input[cipher_num].data() + length, input[cipher_num].data() + length, 1, 1);

    size_t index = pow_mod(5, N1 * cipher_num, 2 * N);
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));

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
            temp_a[destination - length] = modulus - input[cipher_num][i];
            temp_b[destination - length] = modulus - input[cipher_num][i + length];
        } else {
            temp_a[destination] = input[cipher_num][i];
            temp_b[destination] = input[cipher_num][i + length];
        }
        // result.a[i] = 0;
    }

    // keyswitch

    decompose(dec_a, temp_a, ellnum, PP, BBg, modulus);

    std::vector<RlweCiphertext> autokey_index = autokey.keyMap.at(index);

    std::vector<uint64_t> temp(length, 0);

    ntts.ComputeForward(temp_b.data(), temp_b.data(), 1, 1);    
    intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp_b.data(), length, modulus);

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
}
    // return in ntt form
    // result.setIsNtt(true);
#endif
}

void evalAutoS(std::vector<RlweCiphertext>& result, RlweCiphertext& input,
              const int32_t N1, const AutoKeyBSGS& autoKey)
{
    // use the first parameters
    int32_t ellnum = autoKey.getEllnumBS();
    uint64_t base = autoKey.getBaseBS();
    uint64_t BBg = autoKey.getBgBS();

    uint64_t length = input.getLength();
    uint64_t modulus = input.getModulus();

    // compute permutation matrix
    std::vector<std::vector<int32_t> > permutations(N1, std::vector<int32_t>(length, 0));
    compute_permutation_matrix(permutations, N1, length);

    std::vector<std::vector<uint64_t> > decomposed(2 * ellnum, std::vector<uint64_t>(N, 0));
    intel::hexl::NTT ntts = autoKey.getNTT();

    //  the input a should be store in coefficient form
    std::vector<uint64_t> temp_input(N);
    if (input.getIsNtt())
    {
        ntts.ComputeInverse(temp_input.data(), input.a.data(), 1, 1);
        decompose(decomposed, temp_input, ellnum, base, BBg, modulus);
        for (size_t i = 1; i < N1; i++)
        {
            // std::copy(input.b.begin(), input.b.end(), result[i - 1].b.begin());
            for (size_t k = 0; k < length; k++)
            {
                result[i - 1].b[k] = input.b[permutations[i][k]];
            }
        }
    } else {
        decompose(decomposed, input.a, ellnum, base, BBg, modulus);
        ntts.ComputeForward(temp_input.data(), input.b.data(), 1, 1);
        for (size_t i = 1; i < N1; i++)
        {
            // std::copy(temp_input.begin(), temp_input.end(), result[i - 1].b.begin());
            for (size_t k = 0; k < length; k++)
            {
                result[i - 1].b[k] = temp_input[permutations[i][k]];
            }
        }
    }

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(decomposed[i].data(), decomposed[i].data(), 1, 1);
    }

    std::vector<uint64_t> autoed(N);
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);
    // RlweCiphertext tempCpher;
    for (size_t i = 1; i < N1; i++)
    {
        std::vector<RlweCiphertext> autoKey_index = autoKey.keyMap.at(pow_mod(5, i, 2 * N));

        for (size_t j = 0; j < ellnum; j++)
        {
            // autoInNtt(autoed, decomposed[j], pow_mod(5, i, 2 * N));
            // TODO: avx2, avx512f
            for (size_t k = 0; k < length; k++)
            {
                autoed[k] = decomposed[j][permutations[i][k]];
            }

            intel::hexl::EltwiseMultMod(temp1.data(), autoed.data(),
                                    autoKey_index[j].a.data(), length, modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed.data(),
                                    autoKey_index[j].b.data(), length, modulus, 1);
            intel::hexl::EltwiseSubMod(result[i - 1].a.data(), result[i - 1].a.data(), temp1.data(), length, modulus);
            intel::hexl::EltwiseSubMod(result[i - 1].b.data(), result[i - 1].b.data(), temp2.data(), length, modulus);
        }
        result[i - 1].setIsNtt(true);
    }
}

inline void automorphic_rns(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,
                const int32_t index, const uint64_t modulus)
{
    // std::vector<uint64_t> temp(N);
    // std::copy(data.begin(), data.end(), temp.begin());
    
    for (size_t i = 0; i < N; i++)
    {
#if N == 2048
        uint64_t destination = (i * index) & 0x0fff; // mod 2N
#elif N == 256
        uint64_t destination = (i * index) & 0x01ff; // mod 2N
#elif N == 4096
        uint64_t destination = (i * index) & 0x1fff; // mod 2N
#else
        uint64_t destination = (i * index) & (2 * length - 1); // we suppose length is a power of 2
#endif
        if (destination >= N)
        {
            // result[destination - N] = ( modulus - input[i]) % modulus;
            uint64_t temp = modulus - input[i];
            result[destination - N] = temp != modulus ? temp: 0; // = temp is also OK?
        } else {
            result[destination] = input[i];
        }
    }
}

void evalAutoSRNS(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
              const int32_t N1, const AutoKeyBSGSRNS& autoKey,
              std::vector<std::vector<int32_t> >& permutations)
{
    // use the first parameters
    int32_t ellnum = autoKey.getEllnumBS();
    uint64_t base = autoKey.getBaseBS();
    uint64_t BBg = autoKey.getBgBS();

    uint64_t length = input[0].getLength();
    uint64_t modulus = input[0].getModulus();
    uint64_t bsModulus = autoKey.getBSModulus();

    std::vector<RlweCiphertext> result_aux(N1 - 1, RlweCiphertext(length, bsModulus));

    // compute permutation matrix

    std::vector<std::vector<uint64_t> > decomposed(ellnum, std::vector<uint64_t>(N, 0));
    std::vector<std::vector<uint64_t> > decomposed_aux(ellnum, std::vector<uint64_t>(N, 0));
    intel::hexl::NTT ntts = autoKey.getNTT();
    intel::hexl::NTT bsNtts = autoKey.getBSNTT();
    //  the input a should be store in coefficient form
    if (input[0].getIsNtt())
    {
        ntts.ComputeInverse(input[0].a.data(), input[0].a.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].a.data(), input[1].a.data(), 1, 1);
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
        ntts.ComputeInverse(input[0].b.data(), input[0].b.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].b.data(), input[1].b.data(), 1, 1);
        input[0].setIsNtt(false);
        input[1].setIsNtt(false);

        // automomorph and permutation, who is faster?
        /**
        for (size_t i = 1; i < N1; i++)
        {
            // std::copy(input.b.begin(), input.b.end(), result[i - 1].b.begin());
            for (size_t k = 0; k < length; k++)
            {
                result[i - 1].b[k] = input.b[permutations[i][k]];
            }
        }
        **/
    } else {
        // TODO: add more
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
        /**
        ntts.ComputeForward(temp_input.data(), input.b.data(), 1, 1);
        for (size_t i = 1; i < N1; i++)
        {
            // std::copy(temp_input.begin(), temp_input.end(), result[i - 1].b.begin());
            for (size_t k = 0; k < length; k++)
            {
                result[i - 1].b[k] = temp_input[permutations[i][k]];
            }
        }
        **/
    }

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(decomposed[i].data(), decomposed[i].data(), 1, 1);
        bsNtts.ComputeForward(decomposed_aux[i].data(), decomposed_aux[i].data(), 1, 1); /// do first
    }

    std::vector<uint64_t> autoed1(N);
    std::vector<uint64_t> autoed2(N);
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);
    // RlweCiphertext tempCpher;
    for (size_t i = 1; i < N1; i++)
    {        
        std::vector<RlweCiphertext> autoKey_index = autoKey.keyMap.at(pow_mod(5, i, 2 * N));

        for (size_t j = 0; j < ellnum; j++)
        {
            // autoInNtt(autoed, decomposed[j], pow_mod(5, i, 2 * N));
            // TODO: avx2, avx512f
            for (size_t k = 0; k < length; k++)
            {
                autoed1[k] = decomposed[j][permutations[i][k]];
                autoed2[k] = decomposed_aux[j][permutations[i][k]];
            }

            intel::hexl::EltwiseMultMod(temp1.data(), autoed1.data(),
                                    autoKey_index[j].a.data(), length, modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed1.data(),
                                    autoKey_index[j].b.data(), length, modulus, 1);
            intel::hexl::EltwiseSubMod(result[i - 1].a.data(), result[i - 1].a.data(), temp1.data(), length, modulus);
            intel::hexl::EltwiseSubMod(result[i - 1].b.data(), result[i - 1].b.data(), temp2.data(), length, modulus);

            // handle bs modulus
            intel::hexl::EltwiseMultMod(temp1.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].a.data(), length, bsModulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].b.data(), length, bsModulus, 1);
            intel::hexl::EltwiseSubMod(result_aux[i - 1].a.data(), result_aux[i - 1].a.data(), temp1.data(), length, bsModulus);
            intel::hexl::EltwiseSubMod(result_aux[i - 1].b.data(), result_aux[i - 1].b.data(), temp2.data(), length, bsModulus);
        }
    }

    // rns modswitch
    for (size_t i = 0; i < N1 - 1; i++)
    {
        ntts.ComputeInverse(result[i].a.data(), result[i].a.data(), 1, 1);
        ntts.ComputeInverse(result[i].b.data(), result[i].b.data(), 1, 1);

        bsNtts.ComputeInverse(result_aux[i].a.data(), result_aux[i].a.data(), 1, 1);
        bsNtts.ComputeInverse(result_aux[i].b.data(), result_aux[i].b.data(), 1, 1);

        // we add (0, b) before the rns modswitch
        // (0, b) - g^-1(a) * evk
        automorphic_rns(temp1, input[0].b, pow_mod(5, i + 1, 2 * N), modulus);
        automorphic_rns(temp2, input[1].b, pow_mod(5, i + 1, 2 * N), bsModulus);
        intel::hexl::EltwiseAddMod(result[i].b.data(), result[i].b.data(), temp1.data(), length, modulus);
        intel::hexl::EltwiseAddMod(result_aux[i].b.data(), result_aux[i].b.data(), temp2.data(), length, bsModulus);

        // rns modswitch:
        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result[i].a.data(), result[i].a.data(), result_aux[i].a.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[i].a.data(), result[i].a.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[i].a.data(), result[i].a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result[i].b.data(), result[i].b.data(), result_aux[i].b.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[i].b.data(), result[i].b.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[i].b.data(), result[i].b.data(), 1, 1);

        result[i].setIsNtt(true);
    }

    // rns modswitch for input
    intel::hexl::EltwiseSubMod(input[0].a.data(), input[0].a.data(), input[1].a.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(input[0].a.data(), input[0].a.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(input[0].a.data(), input[0].a.data(), 1, 1);

    // handle b
    intel::hexl::EltwiseSubMod(input[0].b.data(), input[0].b.data(), input[1].b.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(input[0].b.data(), input[0].b.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(input[0].b.data(), input[0].b.data(), 1, 1);
    input[0].setIsNtt(true);
}

// the variant of evalAutoSRNS. the key-switching keys would be much smaller 
void evalAutoSRNSFirst(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& intermediate,
              std::vector<RlweCiphertext>& input,
              const int32_t N1, const AutoKeyBSGSRNS& autoKey,
              std::vector<std::vector<int32_t> >& permutations)
{
    // use the first parameters
    int32_t ellnum = autoKey.getEllnumBS();
    uint64_t base = autoKey.getBaseBS();
    uint64_t BBg = autoKey.getBgBS();

    uint64_t length = input[0].getLength();
    uint64_t modulus = input[0].getModulus();
    uint64_t bsModulus = autoKey.getBSModulus();

    std::vector<RlweCiphertext> result_aux(N1 - 1, RlweCiphertext(length, bsModulus));

    // compute permutation matrix

    std::vector<std::vector<uint64_t> > decomposed(ellnum, std::vector<uint64_t>(N, 0));
    std::vector<std::vector<uint64_t> > decomposed_aux(ellnum, std::vector<uint64_t>(N, 0));
    intel::hexl::NTT ntts = autoKey.getNTT();
    intel::hexl::NTT bsNtts = autoKey.getBSNTT();
    //  the input a should be store in coefficient form
    if (input[0].getIsNtt())
    {
        ntts.ComputeInverse(input[0].a.data(), input[0].a.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].a.data(), input[1].a.data(), 1, 1);
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
        ntts.ComputeInverse(input[0].b.data(), input[0].b.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].b.data(), input[1].b.data(), 1, 1);
        input[0].setIsNtt(false);
        input[1].setIsNtt(false);

    } else {
        // TODO: add more
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
    }

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(decomposed[i].data(), decomposed[i].data(), 1, 1);
        bsNtts.ComputeForward(decomposed_aux[i].data(), decomposed_aux[i].data(), 1, 1); /// do first
    }

    std::vector<uint64_t> autoed1(N);
    std::vector<uint64_t> autoed2(N);
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);
    // RlweCiphertext tempCpher;
    for (size_t i = 1; i < N1; i++)
    {        
        std::vector<RlweCiphertext> autoKey_index = autoKey.keyMap.at(pow_mod(5, i, 2 * N));

        for (size_t j = 0; j < ellnum; j++)
        {
            // autoInNtt(autoed, decomposed[j], pow_mod(5, i, 2 * N));
            // TODO: avx2, avx512f
            for (size_t k = 0; k < length; k++)
            {
                autoed1[k] = decomposed[j][permutations[i][k]];
                autoed2[k] = decomposed_aux[j][permutations[i][k]];
            }

            intel::hexl::EltwiseMultMod(temp1.data(), autoed1.data(),
                                    autoKey_index[j].a.data(), length, modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed1.data(),
                                    autoKey_index[j].b.data(), length, modulus, 1);
            intel::hexl::EltwiseSubMod(result[i - 1].a.data(), result[i - 1].a.data(), temp1.data(), length, modulus);
            intel::hexl::EltwiseSubMod(result[i - 1].b.data(), result[i - 1].b.data(), temp2.data(), length, modulus);

            // handle bs modulus
            intel::hexl::EltwiseMultMod(temp1.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].a.data(), length, bsModulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].b.data(), length, bsModulus, 1);
            intel::hexl::EltwiseSubMod(result_aux[i - 1].a.data(), result_aux[i - 1].a.data(), temp1.data(), length, bsModulus);
            intel::hexl::EltwiseSubMod(result_aux[i - 1].b.data(), result_aux[i - 1].b.data(), temp2.data(), length, bsModulus);
        }
    }

    // output intermediate result
    /**
    intermediate[0] = result[N1 - 2];
    intermediate[1] = result_aux[N1 - 2];
    ntts.ComputeInverse(intermediate[0].a.data(), intermediate[0].a.data(), 1, 1);
    ntts.ComputeInverse(intermediate[0].b.data(), intermediate[0].b.data(), 1, 1);
    bsNtts.ComputeInverse(intermediate[1].a.data(), intermediate[1].a.data(), 1, 1);
    bsNtts.ComputeInverse(intermediate[1].b.data(), intermediate[1].b.data(), 1, 1);
    // we add (0, b) before the rns modswitch
    // (0, b) - g^-1(a) * evk
    automorphic_rns(temp1, input[0].b, pow_mod(5, N1 - 1, 2 * N), modulus);
    automorphic_rns(temp2, input[1].b, pow_mod(5, N1 - 1, 2 * N), bsModulus);
    intel::hexl::EltwiseAddMod(intermediate[0].b.data(), intermediate[0].b.data(), temp1.data(), length, modulus);
    intel::hexl::EltwiseAddMod(intermediate[1].b.data(), intermediate[1].b.data(), temp2.data(), length, bsModulus);
    **/

    // rns modswitch
    for (size_t i = 0; i < N1 - 1; i++)
    {
        ntts.ComputeInverse(result[i].a.data(), result[i].a.data(), 1, 1);
        ntts.ComputeInverse(result[i].b.data(), result[i].b.data(), 1, 1);

        bsNtts.ComputeInverse(result_aux[i].a.data(), result_aux[i].a.data(), 1, 1);
        bsNtts.ComputeInverse(result_aux[i].b.data(), result_aux[i].b.data(), 1, 1);

        // we add (0, b) before the rns modswitch
        // (0, b) - g^-1(a) * evk
        automorphic_rns(temp1, input[0].b, pow_mod(5, i + 1, 2 * N), modulus);
        automorphic_rns(temp2, input[1].b, pow_mod(5, i + 1, 2 * N), bsModulus);
        intel::hexl::EltwiseAddMod(result[i].b.data(), result[i].b.data(), temp1.data(), length, modulus);
        intel::hexl::EltwiseAddMod(result_aux[i].b.data(), result_aux[i].b.data(), temp2.data(), length, bsModulus);

        if (i == N1 - 2)
        {
            intermediate[0] = result[N1 - 2];
            intermediate[1] = result_aux[N1 - 2];
        }

        // rns modswitch:
        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result[i].a.data(), result[i].a.data(), result_aux[i].a.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[i].a.data(), result[i].a.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[i].a.data(), result[i].a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result[i].b.data(), result[i].b.data(), result_aux[i].b.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[i].b.data(), result[i].b.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[i].b.data(), result[i].b.data(), 1, 1);

        result[i].setIsNtt(true);
    }

    // rns modswitch for input
    intel::hexl::EltwiseSubMod(input[0].a.data(), input[0].a.data(), input[1].a.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(input[0].a.data(), input[0].a.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(input[0].a.data(), input[0].a.data(), 1, 1);

    // handle b
    intel::hexl::EltwiseSubMod(input[0].b.data(), input[0].b.data(), input[1].b.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(input[0].b.data(), input[0].b.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(input[0].b.data(), input[0].b.data(), 1, 1);
    input[0].setIsNtt(true);
}

// the variant of evalAutoSRNS. the key-switching keys would be much smaller 
void evalAutoSRNSSecond(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
              const int32_t N1, const AutoKeyBSGSRNS& autoKey,
              std::vector<std::vector<int32_t> >& permutations)
{
    // use the first parameters
    int32_t ellnum = autoKey.getEllnumBS();
    uint64_t base = autoKey.getBaseBS();
    uint64_t BBg = autoKey.getBgBS();

    uint64_t length = input[0].getLength();
    uint64_t modulus = input[0].getModulus();
    uint64_t bsModulus = autoKey.getBSModulus();

    std::vector<RlweCiphertext> result_aux(N1, RlweCiphertext(length, bsModulus));

    // compute permutation matrix

    std::vector<std::vector<uint64_t> > decomposed(ellnum, std::vector<uint64_t>(N, 0));
    std::vector<std::vector<uint64_t> > decomposed_aux(ellnum, std::vector<uint64_t>(N, 0));
    intel::hexl::NTT ntts = autoKey.getNTT();
    intel::hexl::NTT bsNtts = autoKey.getBSNTT();
    //  the input a should be store in coefficient form
    if (input[0].getIsNtt())
    {
        ntts.ComputeInverse(input[0].a.data(), input[0].a.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].a.data(), input[1].a.data(), 1, 1);
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
        ntts.ComputeInverse(input[0].b.data(), input[0].b.data(), 1, 1);
        bsNtts.ComputeInverse(input[1].b.data(), input[1].b.data(), 1, 1);
        input[0].setIsNtt(false);
        input[1].setIsNtt(false);

    } else {
        // TODO: add more
        decompose_bsgs(decomposed, decomposed_aux, input[0].a, input[1].a, ellnum, base, BBg);
    }

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(decomposed[i].data(), decomposed[i].data(), 1, 1);
        bsNtts.ComputeForward(decomposed_aux[i].data(), decomposed_aux[i].data(), 1, 1); /// do first
    }

    std::vector<uint64_t> autoed1(N);
    std::vector<uint64_t> autoed2(N);
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);
    // RlweCiphertext tempCpher;
    for (size_t i = 0; i < N1; i++)
    {        
        std::vector<RlweCiphertext> autoKey_index = autoKey.keyMap.at(pow_mod(5, i + 1, 2 * N));

        for (size_t j = 0; j < ellnum; j++)
        {
            // autoInNtt(autoed, decomposed[j], pow_mod(5, i, 2 * N));
            // TODO: avx2, avx512f
            for (size_t k = 0; k < length; k++)
            {
                autoed1[k] = decomposed[j][permutations[i + 1][k]];
                autoed2[k] = decomposed_aux[j][permutations[i + 1][k]];
            }

            intel::hexl::EltwiseMultMod(temp1.data(), autoed1.data(),
                                    autoKey_index[j].a.data(), length, modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed1.data(),
                                    autoKey_index[j].b.data(), length, modulus, 1);
            intel::hexl::EltwiseSubMod(result[N1 + i - 1].a.data(), result[N1 + i - 1].a.data(), temp1.data(), length, modulus);
            intel::hexl::EltwiseSubMod(result[N1 + i - 1].b.data(), result[N1 + i - 1].b.data(), temp2.data(), length, modulus);

            // handle bs modulus
            intel::hexl::EltwiseMultMod(temp1.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].a.data(), length, bsModulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), autoed2.data(),
                                    autoKey_index[j + ellnum].b.data(), length, bsModulus, 1);
            intel::hexl::EltwiseSubMod(result_aux[i].a.data(), result_aux[i].a.data(), temp1.data(), length, bsModulus);
            intel::hexl::EltwiseSubMod(result_aux[i].b.data(), result_aux[i].b.data(), temp2.data(), length, bsModulus);
        }
    }

    // rns modswitch
    for (size_t i = 0; i < N1; i++)
    {
        size_t rindex = N1 + i - 1;
        
        ntts.ComputeInverse(result[rindex].a.data(), result[rindex].a.data(), 1, 1);
        ntts.ComputeInverse(result[rindex].b.data(), result[rindex].b.data(), 1, 1);

        bsNtts.ComputeInverse(result_aux[i].a.data(), result_aux[i].a.data(), 1, 1);
        bsNtts.ComputeInverse(result_aux[i].b.data(), result_aux[i].b.data(), 1, 1);

        // we add (0, b) before the rns modswitch
        // (0, b) - g^-1(a) * evk
        automorphic_rns(temp1, input[0].b, pow_mod(5, i + 1, 2 * N), modulus);
        automorphic_rns(temp2, input[1].b, pow_mod(5, i + 1, 2 * N), bsModulus);
        intel::hexl::EltwiseAddMod(result[rindex].b.data(), result[rindex].b.data(), temp1.data(), length, modulus);
        intel::hexl::EltwiseAddMod(result_aux[i].b.data(), result_aux[i].b.data(), temp2.data(), length, bsModulus);

        // rns modswitch:
        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result[rindex].a.data(), result[rindex].a.data(), result_aux[i].a.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[rindex].a.data(), result[rindex].a.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[rindex].a.data(), result[rindex].a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result[rindex].b.data(), result[rindex].b.data(), result_aux[i].b.data(),
                            length, modulus);
        intel::hexl::EltwiseFMAMod(result[rindex].b.data(), result[rindex].b.data(), bsModInv,
                            nullptr, length, modulus, 1);
        ntts.ComputeForward(result[rindex].b.data(), result[rindex].b.data(), 1, 1);

        result[rindex].setIsNtt(true);
    }
}

void evalAutoSRNSSmallKeys(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
              const int32_t N1, const AutoKeyBSGSRNS& autoKey,
              std::vector<std::vector<int32_t> >& permutations)
{
    std::vector<RlweCiphertext> intermediate(1, RlweCiphertext(N, crtMod));
    intermediate.push_back(RlweCiphertext(N, bsMod));
    evalAutoSRNSFirst(result, intermediate, input, N1 / 2, autoKey, permutations);
    evalAutoSRNSSecond(result, intermediate, N1 / 2, autoKey, permutations);
}

/**
 * @brief the input ciphertext rotate N/2 times
 * @param data have store in diagonal and ntt form
 * 
*/
void matrix_vector_mul(RlweCiphertext& result, RlweCiphertext& input,
                        std::vector<std::vector<uint64_t> >& data,
                        const AutoKeyBSGS& autoKey)
{
    RlweCiphertext input_temp = input;
    
    if (!input_temp.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = autoKey.getNTT();
        ntts.ComputeForward(input_temp.a.data(), input_temp.a.data(), 1, 1);
        ntts.ComputeForward(input_temp.b.data(), input_temp.b.data(), 1, 1);
        input_temp.setIsNtt(true);
#endif
    }

    RlweCiphertext rotated_cipher;
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    // the first line
    intel::hexl::EltwiseMultMod(result.a.data(), input_temp.a.data(), data[0].data(), N,
                                bigMod, 1);
    intel::hexl::EltwiseMultMod(result.b.data(), input_temp.b.data(), data[0].data(), N,
                                bigMod, 1);
    /**/
    for (size_t i = 1; i < N/2; i++)
    {
        // rotate(rotated_cipher, input, i)
        evalAuto(rotated_cipher, input, pow_mod(5, i, 2 * N), autoKey,
                 autoKey.getEllnumBS(), autoKey.getBaseBS(), autoKey.getBgBS());

        intel::hexl::EltwiseMultMod(temp1.data(), rotated_cipher.a.data(), data[i].data(), N,
                                    bigMod, 1);
        intel::hexl::EltwiseMultMod(temp2.data(), rotated_cipher.b.data(), data[i].data(), N,
                                    bigMod, 1);
        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp1.data(), N, bigMod);
        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp2.data(), N, bigMod);
    }
    /**/
    result.setIsNtt(true);
}


/**
 * @brief the input ciphertext rotate N/2 times
 * @param data have store in diagonal and ntt form
 * 
*/
void matrix_vector_mul_bsgs(RlweCiphertext& result, RlweCiphertext& input,
                        std::vector<std::vector<uint64_t> >& data,
                        const AutoKeyBSGS& autoKey, const int32_t N1)
{
    RlweCiphertext input_temp = input;
    uint64_t modulus = autoKey.getModulus();
    
    if (!input_temp.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = autoKey.getNTT();
        ntts.ComputeForward(input_temp.a.data(), input_temp.a.data(), 1, 1);
        ntts.ComputeForward(input_temp.b.data(), input_temp.b.data(), 1, 1);
        input_temp.setIsNtt(true);
#endif
    }

    int32_t N2 = N / 2 / N1;

    std::vector<RlweCiphertext> rotated_cipher(N1 - 1, RlweCiphertext(N, autoKey.getModulus()));
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    // the first line
    // intel::hexl::EltwiseMultMod(result.a.data(), input.a.data(), data[0].data(), N,
    //                             bigMod, 1);
    // intel::hexl::EltwiseMultMod(result.b.data(), input.b.data(), data[0].data(), N,
    //                             bigMod, 1);
#ifdef TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
#ifndef FBSGS
    for (size_t i = 1; i < N1; i++)
    {
        // rotate(rotated_cipher, input, i)
        evalAuto(rotated_cipher[i-1], input, pow_mod(5, i, 2 * N), autoKey,
                autoKey.getEllnumBS(), autoKey.getBaseBS(), autoKey.getBgBS());
    }    
#else
    evalAutoS(rotated_cipher, input, N1, autoKey);
    /**
    intel::hexl::NTT ntts = autoKey.getNTT();
    for (size_t i = 0; i < N1 - 1; i++)
    {
        ntts.ComputeInverse(rotated_cipher[i].a.data(), rotated_cipher[i].a.data(), 1, 1);
        ntts.ComputeInverse(rotated_cipher[i].b.data(), rotated_cipher[i].b.data(), 1, 1);
    }
    for (size_t i = 0; i < N1 - 1; i++)
    {
        ntts.ComputeForward(rotated_cipher[i].a.data(), rotated_cipher[i].a.data(), 1, 1);
        ntts.ComputeForward(rotated_cipher[i].b.data(), rotated_cipher[i].b.data(), 1, 1);
    }
    **/
#endif
#ifdef TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " evalAutoS costs " << glapsed.count() << " us." << std::endl;
#endif
#ifdef TIME_COLLECT
    start = std::chrono::high_resolution_clock::now();
#endif
    RlweCiphertext temp_cipher(N, modulus);
    for (size_t i = 0; i < N2; i++)
    {   
        intel::hexl::EltwiseMultMod(temp_cipher.a.data(), input_temp.a.data(), data[i * N1].data(), N,
                                    modulus, 1);
        intel::hexl::EltwiseMultMod(temp_cipher.b.data(), input_temp.b.data(), data[i * N1].data(), N,
                                    modulus, 1);
        for (size_t j = 1; j < N1; j++)
        {
            intel::hexl::EltwiseMultMod(temp1.data(), rotated_cipher[j - 1].a.data(), data[i * N1 + j].data(), N,
                                        modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), rotated_cipher[j - 1].b.data(), data[i * N1 + j].data(), N,
                                        modulus, 1);
            intel::hexl::EltwiseAddMod(temp_cipher.a.data(), temp_cipher.a.data(), temp1.data(), N, modulus);
            intel::hexl::EltwiseAddMod(temp_cipher.b.data(), temp_cipher.b.data(), temp2.data(), N, modulus);
        }
        temp_cipher.setIsNtt(true);
        /**/
        RlweCiphertext temp_cipher1 = temp_cipher;
        if (i != 0) evalAuto(temp_cipher, temp_cipher1, pow_mod(5, N1 * i, 2 * N), autoKey,
                             autoKey.getEllnumGS(), autoKey.getBaseGS(), autoKey.getBgGS());
        /**/
        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp_cipher.a.data(), N, modulus);
        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp_cipher.b.data(), N, modulus);
    }
    result.setIsNtt(true);
#ifdef TIME_COLLECT
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " core costs " << glapsed.count() << " us." << std::endl;
#endif
}

inline uint64_t barrett_raw_u64(uint64_t input, uint64_t const_ratio_1, uint64_t modulus)
{
    unsigned long long tmp[2];
    tmp[1] = static_cast<unsigned long long>(((static_cast<__uint128_t>(input) * static_cast<__uint128_t>(const_ratio_1)) >> 64));

    // Barrett subtraction
    tmp[0] = input - tmp[1] * modulus;

    // One more subtraction is enough
    return (tmp[0] >= modulus ? (tmp[0] - modulus) : (tmp[0]));
}

inline uint64_t barrett_coeff(uint64_t val, size_t n)
{
    // return val % moduli[n];
    return (n == 0) ? barrett_raw_u64(val, cr1_p, crtq1) : barrett_raw_u64(val, cr1_b, crtq2);
    return val;
}

void reorientCipher(uint64_t* cipherbuf, RlweCiphertext& input,
                    std::vector<RlweCiphertext>& rotated_cipher, int32_t N1)
{
    assert(rotated_cipher.size() == N1 - 1);

    for (size_t i = 0; i < N; i++)
    {
        cipherbuf[i * 2 * N1] = barrett_coeff(input.a[i], 0) | barrett_coeff(input.a[i], 1) << 32;
        cipherbuf[i * 2 * N1 + 1] = barrett_coeff(input.b[i], 0) | barrett_coeff(input.b[i], 1) << 32;
        for (size_t j = 1; j < N1; j++)
        {
            uint64_t temp = rotated_cipher[j - 1].a[i];
            cipherbuf[i * 2 * N1 + 2 * j] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
            temp = rotated_cipher[j - 1].b[i];
            cipherbuf[i * 2 * N1 + 2 * j + 1] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
        }
    }
}

void matrix_vector_mul_bsgs_crt(RlweCiphertext& result, RlweCiphertext& input,
                        uint64_t* data, const AutoKeyBSGS& autoKey, const int32_t N1)
{
    RlweCiphertext input_temp = input;
    
    intel::hexl::NTT ntts = autoKey.getNTT(); // current
    if (!input_temp.getIsNtt())
    {
#ifdef INTEL_HEXL
        // intel::hexl::NTT ntts = autoKey.getNTT();
        ntts.ComputeForward(input_temp.a.data(), input_temp.a.data(), 1, 1);
        ntts.ComputeForward(input_temp.b.data(), input_temp.b.data(), 1, 1);
        input_temp.setIsNtt(true);
#endif
    }

    int32_t N2 = N / 2 / N1;

    std::vector<RlweCiphertext> rotated_cipher(N1 - 1, RlweCiphertext(N, autoKey.getModulus()));
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

#ifdef TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
#ifndef FBSGS
    for (size_t i = 1; i < N1; i++)
    {
        // rotate(rotated_cipher, input, i)
        evalAuto(rotated_cipher[i-1], input, pow_mod(5, i, 2 * N), autoKey,
                autoKey.getEllnumBS(), autoKey.getBaseBS(), autoKey.getBgBS());
    }    
#else
    evalAutoS(rotated_cipher, input, N1, autoKey);
#endif
#ifdef TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " evalAutoS costs " << glapsed.count() << " us." << std::endl;
#endif

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N * 2);
    reorientCipher(cipherbuf, input_temp, rotated_cipher, N1);

    // std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(4 * N, 0));
    std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(2 * N, 0));
#ifdef TIME_COLLECT
    start = std::chrono::high_resolution_clock::now();
#endif
    // fastMultiplyQueryByDatabaseDim1(output, data, cipherbuf, N1, N2);
    fastMultiplyQueryByDatabaseDim1InvCRT(output, data, cipherbuf, N1, N2);

    memcpy(result.a.data(), output[0].data(), N * sizeof(uint64_t));
    memcpy(result.b.data(), output[0].data() + N, N * sizeof(uint64_t));
    result.setIsNtt(true);
    // output is stored in ntt form

    evalAutoCRT(result, output, N2, autoKey,
                autoKey.getEllnumGS(), autoKey.getBaseGS(), autoKey.getBgGS());

#ifdef TIME_COLLECT
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " core costs " << glapsed.count() << " us." << std::endl;
#endif
}

void matrix_vector_mul_bsgs_rns_crt(RlweCiphertext& result, std::vector<RlweCiphertext>& input,
                        uint64_t* data, const AutoKeyBSGSRNS& autoKey, const int32_t N1,
                        std::vector<std::vector<int32_t> >& permutations)
{
    int32_t N2 = N / 2 / N1;

    std::vector<RlweCiphertext> rotated_cipher(N1 - 1, RlweCiphertext(N, autoKey.getModulus()));
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

#ifdef TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
    // evalAutoSRNS(rotated_cipher, input, N1, autoKey, permutations);
    evalAutoSRNSSmallKeys(rotated_cipher, input, N1, autoKey, permutations);
#ifdef TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " evalAutoS costs " << glapsed.count() << " us." << std::endl;
#endif

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N * 2);
    reorientCipher(cipherbuf, input[0], rotated_cipher, N1);

int32_t ntimes = 64 * 8;
    // std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(4 * N, 0));
    std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(2 * N, 0));
#ifdef TIME_COLLECT
    start = std::chrono::high_resolution_clock::now();
#endif
for (size_t i = 0; i < ntimes; ++i)
{
    // fastMultiplyQueryByDatabaseDim1(output, data, cipherbuf, N1, N2);
    fastMultiplyQueryByDatabaseDim1InvCRT(output, data, cipherbuf, N1, N2);

    memcpy(result.a.data(), output[0].data(), N * sizeof(uint64_t));
    memcpy(result.b.data(), output[0].data() + N, N * sizeof(uint64_t));
    result.setIsNtt(true);
    // output is stored in ntt form

    evalAutoRNSCRT(result, output, N2, autoKey,
                autoKey.getEllnumGS(), autoKey.getBaseGS(), autoKey.getBgGS());
}
#ifdef TIME_COLLECT
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " core costs " << glapsed.count() << " us." << std::endl;
#endif
}


void matrix_vector_mul_bsgs_rns_crt_large(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
                        uint64_t* data, const AutoKeyBSGSRNS& autoKey, const int32_t N1,
                        std::vector<std::vector<int32_t> >& permutations, int32_t r)
{
    int32_t N2 = N / 2 / N1;

    std::vector<RlweCiphertext> rotated_cipher(N1 - 1, RlweCiphertext(N, autoKey.getModulus()));
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

#ifdef TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
    // evalAutoSRNS(rotated_cipher, input, N1, autoKey, permutations);
    evalAutoSRNSSmallKeys(rotated_cipher, input, N1, autoKey, permutations);
#ifdef TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "   evalAutoS costs " << glapsed.count() << " us." << std::endl;
#endif

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N * 2);
    reorientCipher(cipherbuf, input[0], rotated_cipher, N1);

// int32_t ntimes = 16;
    // std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(4 * N, 0));
    std::vector<std::vector<uint64_t> > output(N2, std::vector<uint64_t>(2 * N, 0));
#ifdef TIME_COLLECT
    start = std::chrono::high_resolution_clock::now();
#endif

    size_t num_words = N * N / 2;
for (size_t i = 0; i < r; ++i)
{
    // fastMultiplyQueryByDatabaseDim1(output, data, cipherbuf, N1, N2);
    fastMultiplyQueryByDatabaseDim1InvCRT(output, data + i * num_words, cipherbuf, N1, N2);

    memcpy(result[i].a.data(), output[0].data(), N * sizeof(uint64_t));
    memcpy(result[i].b.data(), output[0].data() + N, N * sizeof(uint64_t));
    result[i].setIsNtt(true);
    // output is stored in ntt form

    evalAutoRNSCRT(result[i], output, N2, autoKey,
                autoKey.getEllnumGS(), autoKey.getBaseGS(), autoKey.getBgGS());
}
#ifdef TIME_COLLECT
    stop  = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "   core costs " << glapsed.count() << " us." << std::endl;
#endif
}

/**
 * @brief evaluate automorphic transform for 
 * 
 * @param result always store in coeffiecints form
 * @param index 
 * @param autokey store in ntt form
 */
void evalAutoRNSKey(RlweCiphertext& result1, RlweCiphertext& result2,
              RlweCiphertext& input1, RlweCiphertext& input2,
              const int32_t index, const AutoKeyBSGSRNS& autokey,
              const int32_t ellnum, const uint64_t PP, const uint64_t BBg)
{
    uint64_t length = input1.getLength();
    assert(length == input2.getLength());

    uint64_t modulus1 = input1.getModulus();
    uint64_t modulus2 = input2.getModulus();

    intel::hexl::NTT ntts1 = autokey.getNTT();
    intel::hexl::NTT ntts2 = autokey.getBSNTT();

    if (input1.getIsNtt())
    {
#ifdef INTEL_HEXL

        ntts1.ComputeInverse(input1.a.data(), input1.a.data(), 1, 1);
        ntts1.ComputeInverse(input1.b.data(), input1.b.data(), 1, 1);
        input1.setIsNtt(false);
#endif
    }

    if (input2.getIsNtt())
    {
#ifdef INTEL_HEXL
        ntts2.ComputeInverse(input2.a.data(), input2.a.data(), 1, 1);
        ntts2.ComputeInverse(input2.b.data(), input2.b.data(), 1, 1);
        input2.setIsNtt(false);
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
            temp_a1[destination - length] = modulus1 - input1.a[i];
            temp_b1[destination - length] = modulus1 - input1.b[i];
            temp_a2[destination - length] = modulus2 - input2.a[i];
            temp_b2[destination - length] = modulus2 - input2.b[i];
        } else {
            temp_a1[destination] = input1.a[i];
            temp_b1[destination] = input1.b[i];
            temp_a2[destination] = input2.a[i];
            temp_b2[destination] = input2.b[i];
        }
        result1.a[i] = 0;
        result2.a[i] = 0;
    }

    // keyswitch
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));
    std::vector<std::vector<uint64_t> > dec_a2(ellnum, std::vector<uint64_t>(length, 0));
    // assert(modulus2 == bsMod || modulus2 == crtBaMod);
    assert(modulus2 == auxMod || modulus2 == crtBaMod);
    // if (modulus2 == bsMod)
    if (modulus2 == auxMod)
    {
        // decompose_bsgs(dec_a, dec_a2, temp_a1, temp_a2, ellnum, PP, BBg);
        decompose_bsgs_aux(dec_a, dec_a2, temp_a1, temp_a2, ellnum, PP, BBg);
    } else {
        decompose_bsgs_ba(dec_a, dec_a2, temp_a1, temp_a2, ellnum, PP, BBg);
    }

    std::vector<RlweCiphertext> autokey_index = autokey.keyMap.at(index);

    std::vector<uint64_t> temp(length, 0);
#ifdef INTEL_HEXL
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
 * @brief genreating the autokeys from the offline key
 * @param result1 as input, for baby-step
 * @param result2 as input, for gaint-step
 * @param offlineKey  
 * 
*/
void genAutoKeyFromOffline(AutoKeyBSGSRNS& result, AutoKeyBSGSRNS& result1, AutoKeyBSGSRNS& result2,
              const AutoKeyBSGSRNS& offlineKey1, const AutoKeyBSGSRNS& offlineKey2,
              const int32_t N1)
{
    uint64_t length = result1.getLength();
    int32_t N2 = N / 2 / N1;

    int32_t ellnum = result1.getEllnumBS();
    std::vector<RlweCiphertext> key_for_index(ellnum, RlweCiphertext(N, result1.getModulus()));
    for (size_t i = 0; i < ellnum; i++)
    {
        key_for_index.push_back(RlweCiphertext(N, result1.getBSModulus()));
    }

    // handle baby-step autoKeys
    for (size_t i = 1; i < N1 - 1; i++)
    {
        int32_t index = pow_mod(5, i, 2 * N);
        std::vector<RlweCiphertext> input = result1.keyMap.at(index);

        for (size_t j = 0; j < ellnum; j++)
        {
            evalAutoRNSKey(key_for_index[j], key_for_index[j + ellnum],
                        input[j], input[j + ellnum],
                        5, offlineKey1,
                        offlineKey1.getEllnumBS(), offlineKey1.getBaseBS(), offlineKey1.getBgBS());
        }
        // store in result1
        int32_t next_index = pow_mod(5, i + 1, 2 * N);
        result1.keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(next_index, key_for_index));
    }

    // modulus switching
    std::vector<RlweCiphertext> temp1(ellnum, RlweCiphertext(N, result2.getModulus()));
    for (size_t i = 0; i < ellnum; i++)
    {
        temp1.push_back(RlweCiphertext(N, result2.getBSModulus()));
    }
    // modulus switching
    //  q_3^-1 * (a1 - a3) (mod q_1), q_3^-1 * (a2 - a3) (mod q_2) 
    intel::hexl::NTT ntts1 = result1.getNTT();
    intel::hexl::NTT ntts2 = result2.getBSNTT();
    for (size_t i = 0; i < N1 - 1; i++)
    {
        int32_t index = pow_mod(5, i + 1, 2 * N);
        std::vector<RlweCiphertext> input = result1.keyMap.at(index);
        for (size_t j = 0; j < ellnum; j++)
        {
            std::vector<uint64_t> mod2(N, 0);
            std::vector<uint64_t> mod3(N, 0);
            intel::hexl::EltwiseReduceMod(mod2.data(), input[j + ellnum].a.data(), length, bsMod, 1, 1);
            intel::hexl::EltwiseReduceMod(mod3.data(), input[j + ellnum].a.data(), length, auxMod, 1, 1);
            // q_3^-1 * (a1 - a3) (mod q_1)
            intel::hexl::EltwiseSubMod(temp1[j].a.data(), input[j].a.data(), mod3.data(),
                                length, crtMod);
            intel::hexl::EltwiseFMAMod(temp1[j].a.data(), temp1[j].a.data(), auxModInv,
                                nullptr, length, crtMod, 1);
            ntts1.ComputeForward(temp1[j].a.data(), temp1[j].a.data(), 1, 1);

            // q_3^-1 * (a2 - a3) (mod q_2)
            intel::hexl::EltwiseSubMod(temp1[j + ellnum].a.data(), mod2.data(), mod3.data(),
                                length, bsMod);
            intel::hexl::EltwiseFMAMod(temp1[j + ellnum].a.data(), temp1[j + ellnum].a.data(), auxModInvBs,
                                nullptr, length, bsMod, 1);
            ntts2.ComputeForward(temp1[j + ellnum].a.data(), temp1[j + ellnum].a.data(), 1, 1);

            // handle b
            intel::hexl::EltwiseReduceMod(mod2.data(), input[j + ellnum].b.data(), length, bsMod, 1, 1);
            intel::hexl::EltwiseReduceMod(mod3.data(), input[j + ellnum].b.data(), length, auxMod, 1, 1);
            // q_3^-1 * (b1 - b3) (mod q_1)
            intel::hexl::EltwiseSubMod(temp1[j].b.data(), input[j].b.data(), mod3.data(),
                                length, crtMod);
            intel::hexl::EltwiseFMAMod(temp1[j].b.data(), temp1[j].b.data(), auxModInv,
                                nullptr, length, crtMod, 1);
            ntts1.ComputeForward(temp1[j].b.data(), temp1[j].b.data(), 1, 1);
            // q_3^-1 * (b2 - b3) (mod q_2)
            intel::hexl::EltwiseSubMod(temp1[j + ellnum].b.data(), mod2.data(), mod3.data(),
                                length, bsMod);
            intel::hexl::EltwiseFMAMod(temp1[j + ellnum].b.data(), temp1[j + ellnum].b.data(), auxModInvBs,
                                nullptr, length, bsMod, 1);
            ntts2.ComputeForward(temp1[j + ellnum].b.data(), temp1[j + ellnum].b.data(), 1, 1);

            temp1[j].setIsNtt(true);
            temp1[j + ellnum].setIsNtt(true);
        }
        result.keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, temp1));
    }

    // handle gaint-step autoKeys
    ellnum = result2.getEllnumBS();
    // std::vector<RlweCiphertext> key_for_index2(2 * ellnum);
    std::vector<RlweCiphertext> key_for_index2(ellnum, RlweCiphertext(N, result2.getModulus()));
    for (size_t i = 0; i < ellnum; i++)
    {
        key_for_index2.push_back(RlweCiphertext(N, result2.getBSModulus()));
    }
    for (size_t i = 1; i < N2 - 1; i++)
    {
        int32_t index = pow_mod(5, i * N1, 2 * N);
        std::vector<RlweCiphertext> input = result2.keyMap.at(index);

        for (size_t j = 0; j < ellnum; j++)
        {
            evalAutoRNSKey(key_for_index2[j], key_for_index2[j + ellnum],
                        input[j], input[j + ellnum],
                        pow_mod(5, N1, 2 * N), offlineKey2,
                        offlineKey2.getEllnumBS(), offlineKey2.getBaseBS(), offlineKey2.getBgBS());
        }
        // store in result2
        int32_t next_index = pow_mod(5, (i + 1) * N1, 2 * N);
        result2.keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(next_index, key_for_index2));
    }

    // modulus switching
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    std::vector<RlweCiphertext> temp2(ellnum, RlweCiphertext(N, result2.getModulus()));
    for (size_t i = 0; i < N2 - 1; i++)
    {
        int32_t index = pow_mod(5, (i + 1) * N1, 2 * N);
        std::vector<RlweCiphertext> input = result2.keyMap.at(index);
        for (size_t j = 0; j < ellnum; j++)
        {
            intel::hexl::EltwiseSubMod(temp2[j].a.data(), input[j].a.data(), input[j + ellnum].a.data(),
                                length, crtMod);
            intel::hexl::EltwiseFMAMod(temp2[j].a.data(), temp2[j].a.data(), auxModInv,
                                nullptr, length, crtMod, 1);
            ntts1.ComputeForward(temp2[j].a.data(), temp2[j].a.data(), 1, 1);

            // handle b
            intel::hexl::EltwiseSubMod(temp2[j].b.data(), input[j].b.data(), input[j + ellnum].b.data(),
                                length, crtMod);
            intel::hexl::EltwiseFMAMod(temp2[j].b.data(), temp2[j].b.data(), auxModInv,
                                nullptr, length, crtMod, 1);
            ntts1.ComputeForward(temp2[j].b.data(), temp2[j].b.data(), 1, 1);
            temp2[j].setIsNtt(true);
        }
        result.keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, temp2));
    }
}
