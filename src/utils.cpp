#include <iostream>
#include <math.h>
#include <cassert>

#include "params.h"
#include "utils.h"
#include "crt.h"

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
uint64_t QInv(int lwenum, uint64_t modulus)
{
    uint64_t result = 0;

    if (modulus == bigMod)
    {
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
    case 16: result = 1055531162649601;
        break;
    case 32: result = 1090715534737921;
        break;
    case 64: result = 1108307720782081;
        break;
    case 256: result = 1121501860315201;
        break;
    case 512: result = 1123700883570721;
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
    } else if (modulus == crtMod) {
        switch (lwenum)
        {
        case 1: result = 1;
            break;
        case 2: result = 33487344869801985;
            break;
        case 4: result = 50231017304702977;
            break;
        case 8: result = 58602853522153473;
            break;
        case 16: result = 62788771630878721;
            break;
        case 32: result = 64881730685241345;
            break;
        case 64: result = 65928210212422657;
            break;
        case 256: result = 66713069857808641;
            break;
        case 512: result = 66843879798706305;
            break;
        case 4096: result = 66958338496991761;
            break;
        default: result = 0;
            break;
        }
    } else {
        result = 0;
    }

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

void multConst(std::vector<uint64_t>& result, uint64_t cosntNum, uint64_t modulus)
{
#ifdef INTEL_HEXL
    intel::hexl::EltwiseFMAMod(result.data(), result.data(), cosntNum,
                nullptr, N, modulus, 1);
#endif
    // for (auto iter = result.begin(), iter != result.end()
}

void lweToRlwe(std::vector<uint64_t>& result, uint64_t modulus)
{
    uint64_t length = result.size();

    // uint64_t modulus = bigMod;

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

uint128_t compute_mult_qinvqi_bsgs(uint64_t input, uint64_t qinv, uint64_t qi, uint64_t modulus2 = bsMod)
{
    uint128_t modulus = static_cast<uint128_t>(crtMod) * modulus2;

    uint128_t two128_mod = (static_cast<uint128_t>(two128_mod_hignbits_bs) << 63) + 
                static_cast<uint128_t>(two128_mod_lowbits_bs);
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

uint128_t compute_mult_qinvqi_bsgs_ba(uint64_t input, uint64_t qinv, uint64_t qi, uint64_t modulus2 = bsMod)
{
    uint128_t modulus = static_cast<uint128_t>(crtMod) * modulus2;

    uint128_t two128_mod = (static_cast<uint128_t>(two128_mod_hignbits_bsaux) << 63) + 
                static_cast<uint128_t>(two128_mod_lowbits_bsaux);
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

uint128_t compute_mult_qinvqi_bsgs_aux(uint64_t input, uint64_t qinv, uint64_t qi, uint64_t modulus2 = auxMod)
{
    uint128_t modulus = static_cast<uint128_t>(crtMod) * modulus2;

    uint128_t two128_mod = (static_cast<uint128_t>(two128_mod_hignbits_aux) << 63) + 
                static_cast<uint128_t>(two128_mod_lowbits_aux);
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

/**
 * input[i] = input1[i] * q2inv * q2 + input2[i] * q1inv * q1 (mod q1*q2)
*/
void crt_inv_bsgs(std::vector<uint128_t>& result, const std::vector<uint64_t>& input1,
            const std::vector<uint64_t>& input2, uint64_t modulus1, uint64_t modulus2)
{
    assert(modulus1 == crtMod);
    assert(modulus2 == bsMod);
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    uint64_t q2inv = bsModInv;
    uint64_t q1inv = crtModInv;

    for (size_t i = 0; i < N; i++)
    {
        result[i] = compute_mult_qinvqi_bsgs(input1[i], q2inv, modulus2);
        result[i] += compute_mult_qinvqi_bsgs(input2[i], q1inv, modulus1);
        result[i] %= modulus; 
    }
}

// aux modulus variant of crt_inv_bsgs
void crt_inv_bsgs_aux(std::vector<uint128_t>& result, const std::vector<uint64_t>& input1,
            const std::vector<uint64_t>& input2, uint64_t modulus1, uint64_t modulus2)
{
    assert(modulus1 == crtMod);
    assert(modulus2 == auxMod);
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    // uint64_t q2inv = bsModInv;
    // uint64_t q1inv = crtModInv;
    uint64_t q2inv = auxModInv;
    uint64_t q1inv = crtModInvAux;

    for (size_t i = 0; i < N; i++)
    {
        result[i] = compute_mult_qinvqi_bsgs_aux(input1[i], q2inv, modulus2);
        result[i] += compute_mult_qinvqi_bsgs_aux(input2[i], q1inv, modulus1);
        result[i] %= modulus; 
    }
}

void crt_inv_bsgs_ba(std::vector<uint128_t>& result, const std::vector<uint64_t>& input1,
            const std::vector<uint64_t>& input2, uint64_t modulus1, uint64_t modulus2)
{
    assert(modulus1 == crtMod);
    assert(modulus2 == crtBaMod);
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    uint64_t q2inv = crtBaModInv;
    uint64_t q1inv = crtModInvBa;

    for (size_t i = 0; i < N; i++)
    {
        result[i] = compute_mult_qinvqi_bsgs_ba(input1[i], q2inv, modulus2, modulus2);
        result[i] += compute_mult_qinvqi_bsgs_ba(input2[i], q1inv, modulus1, modulus2);
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


void decompose_openfhe(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
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

// #define DECOMPOSE_VARIANT_A 

#ifdef DECOMPOSE_VARIANT_A
// VARIANT A:
// refer to OpenFHE
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
#else
// VARIANT B:
// refer to tfhe \url{https://github.com/tfhe/tfhe}
void decompose(std::vector<std::vector<uint64_t> >& result, const std::vector<uint64_t>& input,
                int32_t ellnum, uint64_t base, uint64_t BBg, uint64_t modulus)
{
    // TODO: here use directly macro
    uint64_t length = N;
    // uint64_t modulus = bigMod;

    // uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = log2(base);
    if (base == 0)
        baseBits = 0;

    int64_t nativeSubgBits = 64 - gBits;

    // offset = BBg/2 * (multfactors) + base/2
    // We use unsignd digit decompose firstly and then we substract Bg/2 for each decomposed element.
    // Therefore, adding an offset ensures it always correct
    uint64_t offset = 0;
    for (size_t i = 0; i < ellnum; i++)
    {
        offset += (uint64_t)0x01 << (baseBits + i * gBits);
    }
    offset <<= (gBits - 1); // offset *= BBg/2
    offset += base / 2;

    uint64_t d;
    for (size_t i = 0; i < length; i++)
    {
        d = input[i] + offset;
        if (d > modulus)
            // d - modulus < modulus, so the following decomposition is valid
            d -= modulus;

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            // int64_t r = d << nativeSubgBits;
            // r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;
            int64_t r = d << nativeSubgBits >> nativeSubgBits; // least gBits bits
            r -= BBg / 2;

            d >>= gBits;

            if (r >= 0)
                result[j][i] = r;
            else
                result[j][i] = r + modulus;
        }
    }
}
#endif


std::vector<uint64_t> recontruct(std::vector<std::vector<uint64_t> >& dec_a,
                            int32_t ellnum, uint64_t base, uint64_t BBg, uint64_t modulus)
{
    uint64_t length = N;
    
    std::vector<uint64_t> result(length);

    std::vector<uint64_t> mult_factor = powerOfBg(base, BBg, ellnum);;
    for (size_t i = 0; i < length; i++)
    {
        for (size_t j = 0; j < ellnum; j++)
        {
            result[i] += ((uint128_t)mult_factor[j] * dec_a[j][i] % modulus);
            result[i] %= modulus;
        }
    }

    return result;
}

std::vector<uint128_t> recontruct_bsgs(std::vector<std::vector<uint64_t> >& dec_a1,
                            std::vector<std::vector<uint64_t> >& dec_a2,
                            int32_t ellnum, uint64_t base, uint64_t BBg, uint128_t modulus)
{
    uint64_t length = N;
    assert(modulus == crtMod * (uint128_t)bsMod);

    std::vector<uint128_t> result(length);

    std::vector<uint128_t> mult_factor = powerOfBg128(base, BBg, ellnum);;
    for (size_t i = 0; i < length; i++)
    {
        int128_t temp = 0;
        for (size_t j = 0; j < ellnum; j++)
        {
            // decomposed two values is equivalent in the sense of signed decompoition
            int64_t temp1 = dec_a1[j][i];
            int64_t temp2 = dec_a2[j][i];
            if (temp1 != temp2)
            {
                temp1 -= crtMod;
                temp2 -= bsMod;
                assert(temp1 == temp2);
            }

            temp += ((int128_t)mult_factor[j] * temp1);
        }
        // result[i] = (temp % modulus + modulus) % modulus;
        result[i] = (temp + modulus) % modulus;
    }

    return result;
}

void check_recontruct(std::vector<std::vector<uint64_t> >& dec_a, const std::vector<uint64_t>& a,
                int32_t ellnum, uint64_t base, uint64_t BBg, uint64_t modulus)
{
    std::vector<uint64_t> temp = recontruct(dec_a, ellnum, base, BBg, modulus);

    /**
    std::cout << "err = [";
    for (size_t i = 0; i < N - 1; i++)
    {
        int64_t err = (int64_t)a[i] - temp[i];
        std::cout << err << ", ";
    }
    std::cout << (int64_t)a[N-1] - (int64_t)temp[N-1] << "]" << std::endl;
    **/
    std::vector<int64_t> error(N);
    for (size_t i = 0; i < N; i++)
    {
        int64_t err = (int64_t)a[i] - temp[i];
        error.push_back(err);
    }
    std::cout << "err = ";
    showLargeVector(error);
}

void check_recontruct_bsgs(std::vector<std::vector<uint64_t> >& dec_a1, std::vector<std::vector<uint64_t> >& dec_a2,
                const std::vector<uint64_t>& a1, const std::vector<uint64_t>& a2,
                int32_t ellnum, uint64_t base, uint64_t BBg, uint128_t modulus)
{
    std::vector<uint128_t> temp = recontruct_bsgs(dec_a1, dec_a2, ellnum, base, BBg, modulus);

    std::vector<int64_t> error1;
    std::vector<int64_t> error2;
    for (size_t i = 0; i < N; i++)
    {
        int64_t err = (int64_t)a1[i] - temp[i] % crtMod;
        error1.push_back(err);

        err = (int64_t)a2[i] - temp[i] % bsMod;
        error2.push_back(err);
    }
    std::cout << "err1 = ";
    showLargeVector(error1);

    std::cout << "err2 = ";
    showLargeVector(error2);
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
                const std::vector<uint64_t>& input2, int32_t ellnum, uint64_t base, uint64_t BBg,
                uint64_t modulus)
{
    // TODO: here use directly macro
    uint64_t length = N;
    // uint64_t modulus = bigMod;

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

/**
 * @brief decompose for rns RLWE ciphertexts(two part), not rns-decompose
*/
// VARIANT 2
void decompose_bsgs(std::vector<std::vector<uint64_t> >& result1, std::vector<std::vector<uint64_t> >& result2,
                const std::vector<uint64_t>& input1, const std::vector<uint64_t>& input2,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus1 = crtMod;
    uint64_t modulus2 = bsMod;
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    std::vector<uint128_t> input(length);
    crt_inv_bsgs(input, input1, input2, modulus1, modulus2);

    // uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = 0;
    if (base != 0)
    {
        baseBits = log2(base);
    }

    int32_t nativeSubgBits = 128 - gBits;

    // offset = BBg/2 * (multfactors) + base/2
    // We use unsignd digit decompose firstly and then we substract Bg/2 for each decomposed element.
    // Therefore, adding an offset ensures it always correct
    uint128_t offset = 0;
    for (size_t i = 0; i < ellnum; i++)
    {
        offset += (uint128_t)0x01 << (baseBits + i * gBits);
    }
    offset <<= (gBits - 1); // offset *= BBg/2
    offset += base / 2;

    uint128_t d;
    for (size_t i = 0; i < length; i++)
    {
        d = input[i] + offset;
        if (d > modulus)
            // d - modulus < modulus, so the following decomposition is valid
            d -= modulus;

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            // int64_t r = d << nativeSubgBits;
            // r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;
            int64_t r = d << nativeSubgBits >> nativeSubgBits; // least gBits bits
            r -= BBg / 2;

            d >>= gBits;

            if (r >= 0)
            {
                result1[j][i] = r;
                result2[j][i] = r;
            } else {
                result1[j][i] = r + modulus1;
                result2[j][i] = r + modulus2;
            }
        }
    }
}

// aux modulus variant of decompose_bsgs
void decompose_bsgs_aux(std::vector<std::vector<uint64_t> >& result1, std::vector<std::vector<uint64_t> >& result2,
                const std::vector<uint64_t>& input1, const std::vector<uint64_t>& input2,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus1 = crtMod;
    uint64_t modulus2 = auxMod;
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    std::vector<uint128_t> input(length);
    crt_inv_bsgs_aux(input, input1, input2, modulus1, modulus2);

    // uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = 0;
    if (base != 0)
    {
        baseBits = log2(base);
    }

    int32_t nativeSubgBits = 128 - gBits;

    // offset = BBg/2 * (multfactors) + base/2
    // We use unsignd digit decompose firstly and then we substract Bg/2 for each decomposed element.
    // Therefore, adding an offset ensures it always correct
    uint128_t offset = 0;
    for (size_t i = 0; i < ellnum; i++)
    {
        offset += (uint128_t)0x01 << (baseBits + i * gBits);
    }
    offset <<= (gBits - 1); // offset *= BBg/2
    offset += base / 2;

    uint128_t d;
    for (size_t i = 0; i < length; i++)
    {
        d = input[i] + offset;
        if (d > modulus)
            // d - modulus < modulus, so the following decomposition is valid
            d -= modulus;

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            // int64_t r = d << nativeSubgBits;
            // r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;
            int64_t r = d << nativeSubgBits >> nativeSubgBits; // least gBits bits
            r -= BBg / 2;

            d >>= gBits;

            if (r >= 0)
            {
                result1[j][i] = r;
                result2[j][i] = r;
            } else {
                result1[j][i] = r + modulus1;
                result2[j][i] = r + modulus2;
            }
        }
    }
}

/**
 * @brief decompose for rns RLWE ciphertexts(two part), not rns-decompose
*/
// VARIANT 2
void decompose_bsgs_ba(std::vector<std::vector<uint64_t> >& result1, std::vector<std::vector<uint64_t> >& result2,
                const std::vector<uint64_t>& input1, const std::vector<uint64_t>& input2,
                int32_t ellnum, uint64_t base, uint64_t BBg)
{
    // TODO: here use directly macro
    uint64_t length = N;
    uint64_t modulus1 = crtMod;
    uint64_t modulus2 = crtBaMod;
    uint128_t modulus = static_cast<uint128_t>(modulus1) * modulus2;

    std::vector<uint128_t> input(length);
    crt_inv_bsgs_ba(input, input1, input2, modulus1, modulus2);

    // uint64_t half_modulus = modulus >> 1;
    int64_t gBits = log2(BBg);
    int64_t baseBits = 0;
    if (base != 0)
    {
        baseBits = log2(base);
    }

    int32_t nativeSubgBits = 128 - gBits;

    // offset = BBg/2 * (multfactors) + base/2
    // We use unsignd digit decompose firstly and then we substract Bg/2 for each decomposed element.
    // Therefore, adding an offset ensures it always correct
    uint128_t offset = 0;
    for (size_t i = 0; i < ellnum; i++)
    {
        offset += (uint128_t)0x01 << (baseBits + i * gBits);
    }
    offset <<= (gBits - 1); // offset *= BBg/2
    offset += base / 2;

    uint128_t d;
    for (size_t i = 0; i < length; i++)
    {
        d = input[i] + offset;
        if (d > modulus)
            // d - modulus < modulus, so the following decomposition is valid
            d -= modulus;

        d >>= baseBits;

        // note: the signed digit decompose is optimal
        for (size_t j = 0; j < ellnum; j++)
        {
            // Faster variant(refer to OpenFHE) 
            // int64_t r = d << nativeSubgBits;
            // r >>= nativeSubgBits;

            // in theory
            // int64_t r = d & ((0x01 << gBits) - 1);
            // if (r >= BBg/2) r -= BBg;
            int64_t r = d << nativeSubgBits >> nativeSubgBits; // least gBits bits
            r -= BBg / 2;

            d >>= gBits;

            if (r >= 0)
            {
                result1[j][i] = r;
                result2[j][i] = r;
            } else {
                result1[j][i] = r + modulus1;
                result2[j][i] = r + modulus2;
            }
        }
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
 * 
 * @note the input and result must be two vectors
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

/**
 * @brief a ^ b % mod_number
 * we suppose that mod_number is a power of 2
*/
int32_t pow_mod(int32_t a, int32_t b, int32_t mod_number)
{
    int32_t result = 1;
    for (size_t i = 0; i < b; i++)
    {
        // result = result * a % mod_number;
        // result = (result * a) & (mod_number - 1);
        // result = (result * a) & (mod_number - 1);
        result *= a;
    }
    
    // return static_cast<int32_t>(pow(a, b)) % mod_number;
    return result & (mod_number - 1);
}

uint64_t pow_mod(uint64_t a, int32_t b, uint64_t mod_number)
{
    uint64_t result = 1;
    for (size_t i = 0; i < b; i++)
    {
        result = static_cast<uint128_t>(result) * a % mod_number;
    }

    return result;
}


void compute_hexl_rotate_indexes(std::vector<int32_t>& hexl_ntt_index,
                                std::vector<int32_t>& rotate_index, const int32_t length)
{
    // the ntts.computeForward outputs w^hexl_ntt_index[0], w^hexl_ntt_index[1], \cdots
    // hexl_ntt_index = [1, 4097, 2049, 6145, 1025, 5121, 3073, \cdots]
    // std::vector<uint64_t> hexl_ntt_index(length, 0);
    hexl_ntt_index[0] = 1;
    for (size_t i = 0; i < log2(length); i++)
    {
        int32_t current_fill = 0x01 << i;
        int32_t interval = length / current_fill;
        for (size_t j = 0; j < current_fill; j++)
        {
            hexl_ntt_index[j * interval + interval/2] = hexl_ntt_index[j * interval] *
                                         (2 * current_fill + 1) % (4 * current_fill);
        }
    }

    // rotate_index = [1, 5, 25, 125, 625, \cdots, q-1, q-5, \cdots]
    // std::vector<uint64_t> rotate_index(length, 0);
    /**
    int32_t half_length = length / 2;
    rotate_index[0] = 1;
    rotate_index[half_length] = (2 * length) - 1;
    for (size_t i = 1; i < length/2; i++)
    {
        rotate_index[i] = rotate_index[i - 1] * 5 % (2 * length);
        rotate_index[half_length + i] = rotate_index[half_length + i - 1] * 5 % (2 * length);
    }
    **/
    // int32_t half_length = length / 2;
    rotate_index[0] = 1;
    rotate_index[length - 1] = (2 * length) - 1;
    for (size_t i = 1; i < length/2; i++)
    {
        rotate_index[i] = rotate_index[i - 1] * 5 % (2 * length);
        // rotate_index[half_length + i] = rotate_index[half_length + i - 1] * 5 % (2 * length);
        rotate_index[length - 1 - i] = rotate_index[length - 1 - i + 1] * 5 % (2 * length);
    }
}


void compute_find_index(std::vector<int32_t>& hexl_ntt_index, std::vector<int32_t>& find_index, const int32_t length)
{
    // std::vector<int32_t> hexl_ntt_index(length, 0);
    std::vector<int32_t> rotate_index(length, 0);
    compute_hexl_rotate_indexes(hexl_ntt_index, rotate_index, length);

    // [0, 1,  2, 3, 4, \cdots, 2048, \cdots]
    // [1, 3,  5, 7, 9, \cdots, 4097, \cdots]
    // [0, *,  *, *, *, \cdots,    1, \cdots]
    for (size_t i = 0; i < length; i++)
    {
        find_index[hexl_ntt_index[i] >> 0x01] = i;
    }
}

inline void compute_inter_permutation(std::vector<int32_t>& inter_permutation,
                               const std::vector<int32_t>& rotate_index, int32_t index, int32_t length = N)
{
    // permutation implicts the number of rotatition
    // [0, 1,  2, 3, 4, \cdots]
    // [1, 3,  5, 7, 9, \cdots]
    // [5, *, 25, \cdots] for one rotation
    // int32_t index = 1;
    // std::vector<uint64_t> inter_permutation(length, 0);
    // int32_t half_length = length / 2;
    for (size_t i = 0; i < length/2 - index; i++)
    {
        inter_permutation[rotate_index[i] >> 0x01] = rotate_index[i + index];
        // inter_permutation[rotate_index[half_length + i] >> 0x01] = rotate_index[half_length + i + index];
        inter_permutation[rotate_index[length - 1 - i] >> 0x01] = rotate_index[length - 1 - i - index];
    }
    for (size_t i = length/2 - index; i < length/2; i++)
    {
        inter_permutation[rotate_index[i] >> 0x01] = rotate_index[i + index - length/2];
        // inter_permutation[rotate_index[half_length + i] >> 0x01] = rotate_index[half_length + i + index - length/2];
        inter_permutation[rotate_index[length - 1 - i] >> 0x01] = rotate_index[length - 1 - i - index + length/2];
    }

    // sort(permutation.data(), permutation.data() + length);
}

/**
 * @brief compuate permutaion for a single index
 * 
*/
void compute_permutation(std::vector<int32_t>& permutation,
                         const int32_t index, const int32_t length)
{
    std::vector<int32_t> hexl_ntt_index(length, 0);
    std::vector<int32_t> rotate_index(length, 0);
    compute_hexl_rotate_indexes(hexl_ntt_index, rotate_index, length);

    // [0, 1,  2, 3, 4, \cdots, 2048, \cdots]
    // [1, 3,  5, 7, 9, \cdots, 4097, \cdots]
    // [0, *,  *, *, *, \cdots,    1, \cdots]
    std::vector<int32_t> find_index(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        find_index[hexl_ntt_index[i] >> 0x01] = i;
    }

    // compute permutation matrix
    std::vector<int32_t> inter_permutation(length, 0);
    compute_inter_permutation(inter_permutation, rotate_index, index, length);

    //  permutation = hexl_ntt_index \circ inter_permutation \circ find_index
    for (size_t i = 0; i < length; i++)
    {
        permutation[i] = find_index[inter_permutation[hexl_ntt_index[i] >> 0x01] >> 0x01];
    }
}

/**
 * @brief compuate permutaion matrix for a max index
 * 
*/
void compute_permutation_matrix(std::vector<std::vector<int32_t> >& permutations, 
                                const int32_t max_index, const int32_t length)
{
    std::vector<int32_t> hexl_ntt_index(length, 0);
    std::vector<int32_t> rotate_index(length, 0);
    compute_hexl_rotate_indexes(hexl_ntt_index, rotate_index, length);

    // [0, 1, 2, 3, 4, \cdots, 2048, \cdots]
    // [1, 3, 5, 7, 9, \cdots, 4097, \cdots]
    // [0, *, *, *, *, \cdots,    1, \cdots]
    std::vector<int32_t> find_index(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        find_index[hexl_ntt_index[i] >> 0x01] = i;
    }

    // compute permutation matrix
    std::vector<int32_t> inter_permutation(length, 0);
    for (size_t index = 0; index < max_index; index++)
    {
        compute_inter_permutation(inter_permutation, rotate_index, index, length);

        //  permutation = hexl_ntt_index \circ inter_permutation \circ find_index
        for (size_t i = 0; i < length; i++)
        {
            permutations[index][i] = find_index[inter_permutation[hexl_ntt_index[i] >> 0x01] >> 0x01];
        }
    }
}

/**
 * @brief compuate interval permutaion matrix for a max index
 * return indexes: 0, N1, 2 * N1, 3 * N1, \cdots, (max_index - 1) * N1
 * 
 */
void compute_interval_permutation_matrix(std::vector<std::vector<int32_t> >& permutations, 
                                const int32_t max_index, const int32_t N1, const int32_t length)
{
    std::vector<int32_t> hexl_ntt_index(length, 0);
    std::vector<int32_t> rotate_index(length, 0);
    compute_hexl_rotate_indexes(hexl_ntt_index, rotate_index, length);

    // [0, 1, 2, 3, 4, \cdots, 2048, \cdots]
    // [1, 3, 5, 7, 9, \cdots, 4097, \cdots]
    // [0, *, *, *, *, \cdots,    1, \cdots]
    std::vector<int32_t> find_index(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        find_index[hexl_ntt_index[i] >> 0x01] = i;
    }

    // compute permutation matrix
    std::vector<int32_t> inter_permutation(length, 0);
    for (size_t index = 0; index < max_index; index++)
    {
        // compute permutation matrix for each index * N1
        compute_inter_permutation(inter_permutation, rotate_index, index * N1, length);

        //  permutation = hexl_ntt_index \circ inter_permutation \circ find_index
        for (size_t i = 0; i < length; i++)
        {
            permutations[index][i] = find_index[inter_permutation[hexl_ntt_index[i] >> 0x01] >> 0x01];
        }
    }
}

void compute_query_encode(std::vector<int32_t>& query_encode, const int32_t length)
{
    // compute hexl_ntt_index and rotate_index
    std::vector<int32_t> hexl_ntt_index(length, 0);
    std::vector<int32_t> rotate_index(length, 0);
    compute_hexl_rotate_indexes(hexl_ntt_index, rotate_index, length);

    std::vector<int32_t> find_index(length, 0);
    for (size_t i = 0; i < length; i++)
    {
        find_index[hexl_ntt_index[i] >> 0x01] = i;
    }

    // std::vector<int32_t> query_encode(length, 0); // find encode positon
    for (size_t i = 0; i < length; i++)
    {
        query_encode[i] = find_index[rotate_index[i] / 2];
    }
}

void compute_query_decode(std::vector<int32_t>& query_decode, const int32_t length)
{
    std::vector<int32_t> query_encode(length, 0);
    compute_query_encode(query_encode, length);
    for (size_t i = 0; i < length; i++)
    {
        query_decode[query_encode[i]] = i;
    }
}
