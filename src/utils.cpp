#include <iostream>
#include <math.h>
#include <cassert>

#include "params.h"
#include "utils.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

void copy_vector(int64_t* result, uint64_t* a)
{
    for (size_t i = 0; i < N; i++)
    {
        result[i] = a[i];
    }
}

void copy_vector(uint64_t* result, int64_t* a, uint64_t modulus)
{
    for (size_t i = 0; i < N; i++)
    {
        // after montgomery reduction, -q < a[i] < q
        result[i] = a[i] < 0 ? a[i] + modulus: a[i];
    }
}

template<typename T1, typename T2>
void copy_vector(T1* result, T2* a)
{
    for (size_t i = 0; i < N; i++)
    {
        result[i] = a[i];
    }
}

template<typename T>
void showLargeVector(std::vector<T>& vals)
{
	std::cout << "[";
	std::cout << vals[0];
	if (vals.size() <= 16)
    {
		for (long i = 1; i < vals.size(); ++i) std::cout << ", " << vals[i];
	} else
    {
		for (long i = 1; i < 16; ++i) std::cout << ", " << vals[i];
		std::cout << ", ..., " << vals[vals.size()-1];
	}
	std::cout << "]" << std::endl;
}

void showLargeVector(std::vector<uint64_t>& vals, std::string ss)
{
	std::cout << ss;
    std::cout << "[";
	std::cout << vals[0];
	if (vals.size() <= 16)
    {
		for (long i = 1; i < vals.size(); ++i) std::cout << ", " << vals[i];
	} else
    {
		for (long i = 1; i < 16; ++i) std::cout << ", " << vals[i];
		std::cout << ", ..., " << vals[vals.size()-1];
	}
	std::cout << "]" << std::endl;
}

void showLargeIntervalVector(std::vector<uint64_t>& vals, int32_t interval, std::string ss)
{
	std::cout << ss;
    std::cout << "[";
	std::cout << vals[0];
	if (vals.size() <= 16 * interval)
    {
		for (long i = interval; i < vals.size(); i += interval) std::cout << ", " << vals[i];
	} else
    {
		for (long i = interval; i < 16 * interval; i += interval) std::cout << ", " << vals[i];
		std::cout << ", ..., " << vals[vals.size() - interval];
	}
	std::cout << "]" << std::endl;
}

template<typename T>
void showVector(std::vector<T>& vals, std::string ss)
{
	std::cout << ss;
    std::cout << "[";
	std::cout << vals[0];
	for (long i = 1; i < vals.size(); ++i)
    {
		std::cout << ", " << vals[i];
	}
	std::cout << "]" << std::endl;
}

void showVector(std::vector<uint64_t>& vals, std::string ss)
{
	std::cout << ss;
    std::cout << "[";
	std::cout << vals[0];
	for (long i = 1; i < vals.size(); ++i)
    {
		std::cout << ", " << vals[i];
	}
	std::cout << "]" << std::endl;
}

/**
 * @brief return inv(lwesnum, Q). Readers can add the following table.
 * 
 * @param lwenum 
 * @return int32_t 
 */
uint64_t QInv(int lwenum)
{
    uint64_t result = 0;

#if bigMod == 1125899906826241
    static_assert(bigMod == 1125899906826241, "the function is related to specific modulus.");
    switch (lwenum)
    {
    case 1: result = 1;
        break;
    case 2: result = 562949953413121;
        break;
    case 4: result = 844424930119681;
        break;
    case 8: result = 985162418472961;
        break;
    case 32: result = 1090715534737921;
        break;
    case 64: result = 1108307720782081;
        break;
    case 256: result = 1121501860315201;
        break;
    case 4096: result = 1125625028919301;
        break;
    default: result = 0;
        break;
    }
#endif

#if bigMod == 18014398509450241
    static_assert(bigMod == 18014398509450241, "the function is related to specific modulus.");
    switch (lwenum)
    {
    case 1: result = 1;
        break;
    case 2: result = 9007199254725121;
        break;
    case 4: result = 13510798882087681;
        break;
    case 8: result = 15762598695768961;
        break;
    case 32: result = 17451448556029921;
        break;
    case 64: result = 17732923532740081;
    default: result = 0;
        break;
    }
#endif

    return result;
}


void database_to_signed(std::vector<std::vector<uint64_t> >& data, uint64_t pbits, uint64_t modulus)
{
    // the database has N * N entries, and each entry has 16 bits
    uint64_t pmodulus = 0x01 << pbits;
    uint64_t halfp = 0x01 << (pbits - 1);

    uint64_t sub = modulus - pmodulus;
    
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (data[i][j] >= halfp)
                data[i][j] += sub;
        }
    }
}

void database_to_rnssigned(std::vector<std::vector<uint64_t> >& data, uint64_t pbits,
                            uint64_t modulus1, uint64_t modulus2)
{
    // the database has N * N entries, and each entry has 16 bits
    uint64_t pmodulus = 0x01 << pbits;
    uint64_t halfp = 0x01 << (pbits - 1);

    // handle crt
    uint64_t sub1 = modulus1 - pmodulus;
    uint64_t sub2 = modulus2 - pmodulus;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (data[i][j] >= halfp)
            {
                data[i + N][j] = data[i][j] + sub2;
                data[i][j] += sub1;
            }
            else
            {
                data[i + N][j] = data[i][j];
            }
        }
    }
}

void database_tontt(std::vector<std::vector<uint64_t> >& data)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(N, bigMod);

    for (size_t i = 0; i < N; i++)
    {
        ntts.ComputeForward(data[i].data(), data[i].data(), 1, 1);
    }
#endif
}

void database_to_rnsntt(std::vector<std::vector<uint64_t> >& data)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1(N, bigMod);
    intel::hexl::NTT ntts2(N, bigMod2);

    for (size_t i = 0; i < N; i++)
    {
        ntts1.ComputeForward(data[i].data(), data[i].data(), 1, 1);
        ntts2.ComputeForward(data[i + N].data(), data[i + N].data(), 1, 1);
    }
#endif
}

void data_to_setupdata(std::vector<std::vector<uint64_t> >& setup_data,
                        std::vector<std::vector<uint64_t> >& data,
                        uint64_t pbits, uint64_t modulus1, uint64_t modulus2)
{
    // copy and tranpose
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            setup_data[j][i] = data[i][j];
        }
    }
    database_to_rnssigned(setup_data, pbits, modulus1, modulus2);
    database_to_rnsntt(setup_data);
}

void negate(uint64_t* result, uint64_t length, uint64_t modulus)
{
    for (size_t i = 0; i < length; i++)
    {
        result[i] = modulus - result[i];
    }
}

void multConst(std::vector<uint64_t>& result, uint64_t cosntNum)
{
#ifdef INTEL_HEXL
    intel::hexl::EltwiseFMAMod(result.data(), result.data(), cosntNum,
                nullptr, N, bigMod, 1);
#endif
    // for (auto iter = result.begin(), iter != result.end()
}

void lweToRlwe(std::vector<uint64_t>& result)
{
    uint64_t length = result.size();
    uint64_t modulus = bigMod;

    std::vector<uint64_t> temp(length);

    copy(result.begin(), result.end(), temp.begin());

    for (size_t i = 1; i < length; i++)
    {   
        result[i] = modulus - temp[length - i];
    }
}

void transpose(std::vector<std::vector<uint64_t> >& a)
{
    if (a.size() != a[0].size())
    {
        std::cout << "Error: the transposed matrix should be square.";
    }

    for (size_t i = 0; i < a.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            uint64_t temp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = temp;
        }
    }

}

std::vector<uint64_t> powerOfBg(uint64_t base, uint64_t BBg, int32_t ellnum)
{
    std::vector<uint64_t> result;
    if (base == 0)
        base = 1;
    result.push_back(base);

    for (size_t i = 0; i < ellnum - 1; i++)
    {
        result.push_back(base *= BBg);
    }

    return result;
}

std::vector<uint128_t> powerOfBg128(uint64_t base, uint64_t BBg, int32_t ellnum)
{
    std::vector<uint128_t> result;
    if (base == 0)
        base = 1;
    result.push_back(base);

    uint128_t temp = base;
    for (size_t i = 0; i < ellnum - 1; i++)
    {
        result.push_back(temp *= BBg);
    }

    return result;
}

uint128_t compute_mult_qinvqi(uint64_t input, uint64_t qinv, uint64_t qi)
{
    
    uint128_t modulus = static_cast<uint128_t>(bigMod) * bigMod2;

    uint128_t two128_mod = (static_cast<uint128_t>(two128_mod_hignbits) << 63) + 
                static_cast<uint128_t>(two128_mod_lowbits);
    uint128_t lowbits_mask = (static_cast<uint128_t>(0x01) << 64) - 1;

    /**
     * (2^64 a + b) * c
     * = 2^64 ac + bc
     * = 2^128 ac_hign + 2^64 ac_low + bc
     * 
     * */ 
    uint128_t input_mult_factor = static_cast<uint128_t>(qi) * qinv;
    uint128_t input_mult_factor_hignbits = input_mult_factor >> 64; // a
    uint128_t input_mult_factor_lowbits = input_mult_factor & lowbits_mask; // b

    uint128_t temp = static_cast<uint128_t>(input) * input_mult_factor_hignbits; // ac

    uint128_t result = (temp >> 64) * static_cast<uint128_t>(two128_mod) % modulus;
    result += (temp << 64) % modulus;
    result += input * input_mult_factor_lowbits;
    result %= modulus;

    return result;
}

// q1 and q2 should not be too large
// TODO: check it
/**
 * input[i] = input1[i] * q2inv * q2 + input2[i] * q1inv * q1 (mod q1*q2)
*/
void crt_inv(std::vector<uint128_t>& result, const std::vector<uint64_t>& input1,
            const std::vector<uint64_t>& input2, uint64_t modulus1, uint64_t modulus2)
{
    assert(modulus1 == bigMod);
    assert(modulus2 == bigMod2);
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    uint64_t q2inv = mod2inv;
    uint64_t q1inv = mod1inv;

    for (size_t i = 0; i < N; i++)
    {
        result[i] = compute_mult_qinvqi(input1[i], q2inv, modulus2);
        result[i] += compute_mult_qinvqi(input2[i], q1inv, modulus1);
        result[i] %= modulus; 
    }
}

/*
// powofBg < modulus
uint64_t powerOfBg(uint64_t base, uint64_t BBg, int32_t num)
{
    for (size_t i = 0; i < num; i++)
    {
        base *= BBg;
    }

    return base;
}
*/

void decompose(uint64_t** result, const uint64_t* input,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus = bigMod;

    uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    if (base == 0)
        baseBits = 0;

    int64_t nativeSubgBits = 64 - gBits;

    int64_t d = 0;
    for (size_t i = 0; i < length; i++)
    {
        if (input[i] > half_modulus)
            d = input[i] - modulus;
        else 
            d = input[i];

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant      
            int64_t r = d << nativeSubgBits;
            r >>= nativeSubgBits;

            d -= r;
            d >>= gBits;

            if (r >= 0)
                result[j][i] = r;
            else
                result[j][i] = r + modulus;
        }
        d = 0;
    }

}


void decompose(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus = bigMod;

    uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    if (base == 0)
        baseBits = 0;

    int64_t nativeSubgBits = 64 - gBits;

    int64_t d = 0;
    for (size_t i = 0; i < length; i++)
    {
        if (input[i] > half_modulus)
            d = input[i] - modulus;
            // d = input[i];
        else 
            d = input[i];

        // TODO: check
        // the remove bits should be signed
        if (base > 0)
        {
            int64_t rr = d & (base - 1);
            if (rr >= base/2) rr -= base;
            d -= rr;
            d >>= baseBits;
        }

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            int64_t r = d << nativeSubgBits;
            r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;

            d -= r;
            d >>= gBits;

            if (r >= 0)
                result[j][i] = r;
            else
                result[j][i] = r + modulus;
        }
        d = 0;
    }

}

std::vector<uint64_t> recontruct(std::vector<std::vector<uint64_t> >& dec_a,
                            int32_t ellnum, uint64_t base, uint64_t BBg)
{
    uint64_t length = N;
    
    std::vector<uint64_t> result(length);

    std::vector<uint64_t> mult_factor = powerOfBg(base, BBg, ellnum);;
    for (size_t i = 0; i < length; i++)
    {
        for (size_t j = 0; j < ellnum; j++)
        {
            result[i] += ((uint128_t)mult_factor[j] * dec_a[j][i] % bigMod);
            result[i] %= bigMod;
        }
    }

    return result;
}

void check_recontruct(std::vector<std::vector<uint64_t> >& dec_a, const std::vector<uint64_t>& a,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    std::vector<uint64_t> temp = recontruct(dec_a, ellnum, base, BBg);

    std::cout << "err = [";
    for (size_t i = 0; i < N - 1; i++)
    {
        int64_t err = (int64_t)a[i] - temp[i];
        std::cout << err << ", ";
    }
    std::cout << (int64_t)a[N-1] - (int64_t)temp[N-1] << "]" << std::endl;
}

// too slow
void decompose_variant(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus = bigMod;

    uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    if (base == 0)
        baseBits = 0;

    int64_t nativeSubgBits = 64 - gBits;

    int64_t d = 0;
    for (size_t i = 0; i < length; i++)
    {
        if (input[i] > half_modulus)
            d = input[i] - modulus;
            // d = input[i];
        else 
            d = input[i];

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant      
            int64_t r = d << nativeSubgBits;
            r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;

            d -= r;
            d >>= gBits;

            if (r >= 0)
                // result[j][i] = r;
            {
                for (size_t k = 0; k < length; k++)
                    result[j * N + i][k] = r;
            }
            else {
                // result[j][i] = r + modulus;
                for (size_t k = 0; k < length; k++)
                    result[j * N + i][k] = r + modulus;
            }
        }
        d = 0;
    }

}

void decompose_rlwe(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input1,
                const std::vector<uint64_t>& input2, int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus = bigMod;

    uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    int64_t nativeSubgBits = 64 - gBits;

    int64_t d = 0;
    // handle input1
        for (size_t i = 0; i < length; i++)
        {
            if (input1[i] > half_modulus)
                d = input1[i] - modulus;
            else
                d = input1[i];

            d >>= baseBits;

            // note: the signed digit decompose is optimal
            for (size_t j = 0; j < ellnum; j++)
            {
                // Faster variant      
                int64_t r = d << nativeSubgBits;
                r >>= nativeSubgBits;

                d -= r;
                d >>= gBits;

                if (r >= 0)
                    result[j][i] = r;
                else
                    result[j][i] = r + modulus;
            }
            d = 0;
        }

    // handle input2
        for (size_t i = 0; i < length; i++)
        {
            if (input2[i] > half_modulus)
                d = input2[i] - modulus;
            else
                d = input2[i];

            d >>= baseBits;

            // note: the signed digit decompose is optimal
            for (size_t j = 0; j < ellnum; j++)
            {
                // Faster variant      
                int64_t r = d << nativeSubgBits;
                r >>= nativeSubgBits;

                d -= r;
                d >>= gBits;

                if (r >= 0)
                    result[ellnum + j][i] = r;
                else
                    result[ellnum + j][i] = r + modulus;
            }
            d = 0;
        }

}

/**
 * @brief decompose for rns RLWE ciphertexts(two part), not rns-decompose
*/
void decompose_crt(std::vector<std::vector<uint64_t> >& result1, std::vector<std::vector<uint64_t> >& result2,
                const std::vector<uint64_t>& input1, const std::vector<uint64_t>& input2,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus1 = bigMod;
    uint64_t modulus2 = bigMod2;
    uint128_t modulus = static_cast<uint128_t>(bigMod) * bigMod2;

    std::vector<uint128_t> input(length);
    crt_inv(input, input1, input2, modulus1, modulus2);

    uint128_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    if (base == 0)
        baseBits = 0;

    int128_t nativeSubgBits = 128 - gBits;

    int128_t d = 0;
    for (size_t i = 0; i < length; i++)
    {
        if (input[i] > half_modulus)
            d = input[i] - modulus;
        else 
            d = input[i];

        // TODO: check
        // the remove bits should be signed
        if (base > 0)
        {
            int128_t rr = d & (base - 1);
            if (rr >= base/2) rr -= base;
            d -= rr;
            d >>= baseBits;
        }

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            int128_t r = d << nativeSubgBits;
            r >>= nativeSubgBits;

            d -= r;
            d >>= gBits;

            if (r >= 0)
            {
                result1[j][i] = r;
                result2[j][i] = r;
            }
            else
            {
                // note: the returned value is in [0, modulus1 - 1] and [0, modulus2 - 1]
                result1[j][i] = r + modulus1;
                result2[j][i] = r + modulus2;
            }

        }
        d = 0;
    }

}

/*
void element_to_vector(std::vector<uint64_t>& result, uint64_t input)
{
    for (auto& re : result)
        re = input;
    // for_each(result.begin(),)    
}
*/
/**
 * @brief 
 * 
 * @param data default length of N
 * @param index 
 * @param modulus
 */
void automorphic(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,
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
        uint64_t destination = (i * index) % (2 * length);
#endif
        if (destination >= N)
        {
            result[destination - N] = (modulus - input[i]) % modulus;
        } else {
            result[destination] = input[i];
        }
    }
}

/*
uint64_t modPowerofN(uint64_t input)
{
    uint64_t result;
#if N == 2048
        uint64_t result = input & 0x0fff; // mod 2N
#elif N == 256
        uint64_t result = input & 0x01ff; // mod 2N
#elif N == 4096
        uint64_t result = input & 0x1fff; // mod 2N
#else
        uint64_t result = input % (2 * N);
#endif

    return result;
}
*/

void encode_crt(std::vector<uint64_t>& result)
{
    // uint128_t temp = static_cast<uint128_t>(bigMod2) * mod2inv % bigMod;
    
    for (auto iter = result.begin(); iter != result.end(); iter++)
    {
        *iter = static_cast<uint128_t>(*iter) * bigMod2 % bigMod;
    }
}

void compute_indicator(int32_t& data_index, bool& reverse, size_t j, int32_t s_index)
{
    if (j > s_index)
    {
        reverse = true;
        data_index = N - j + s_index;
    } else
    {
        reverse = false;
        data_index = s_index - j;
    }
}

/*
int64_t mod_inverse(int64_t a, int64_t b)
{
    // usint b0 = b;
    int64_t t, q;
    int64_t x0 = 0, x1 = 1;
    if (b == 1)
        return 1;
    while (a > 1) {
        q = a / b;
        t = b, b = a % b, a = t;
        t = x0, x0 = x1 - q * x0, x1 = t;
    }
    // if (x1 < 0) x1 += b0;
    // TODO: x1 is never < 0

    return x1;
}
*/
