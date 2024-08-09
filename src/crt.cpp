#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include <chrono> // benchmark
#include <iostream>

// include __m512i
#if defined(__AVX512F__)
#include <immintrin.h>
#include <avx512fintrin.h>
#endif

#if defined(__AVX2__)
#define USE_AVX2
#include <immintrin.h>
#endif

#include "crt.h"
#include "constants.h"


constexpr size_t pol_bytes = 2 * poly_len * sizeof(uint64_t);

constexpr __uint128_t b_inv_pa_i = 163640210UL * crtq2;
constexpr __uint128_t pa_inv_b_i = 97389680UL * crtq1;
constexpr uint64_t cr0_Q = 7906011006380390721UL;
constexpr uint64_t cr1_Q = 275UL; 

void ntt_forward(uint64_t *operand_overall) {
    for (size_t coeff_mod = 0; coeff_mod < 2; coeff_mod++) {
        size_t n = poly_len;

        const uint64_t *forward_tables = &tables[poly_len * 2 * 2 + coeff_mod * poly_len * 2];
        uint64_t *operand = &operand_overall[coeff_mod * poly_len];
        uint32_t modulus_small = coeff_mod == 0 ? crtq1 : crtq2;
        uint32_t two_times_modulus_small = 2 * modulus_small;
        
        for (size_t mm = 0; mm < coeff_count_pow_of_2; mm++) {
            size_t m = 1 << mm;
            size_t t = n >> (mm+1);
            // m*t is always 1024

            for (size_t i = 0; i < m; i++) {
                const uint64_t W = forward_tables[m + i];
                const uint64_t Wprime = forward_tables[poly_len + m + i];

                #ifdef USE_AVX2
                if (t < 4) {
                #endif

                    for (size_t j = 0; j < t; j++) {
                        // The Harvey butterfly: assume X, Y in [0, 2p), and return
                        // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                        
                        uint64_t *p_x = &operand[2 * i * t + j];
                        uint64_t *p_y = &operand[2 * i * t + t + j];
                        uint32_t x = *p_x;
                        uint32_t y = *p_y;

                        uint32_t currX = x - (two_times_modulus_small * static_cast<int>(x >= two_times_modulus_small));

                        uint64_t Q = ((y) * (uint64_t)(Wprime)) >> 32;

                        Q = W * y - Q * modulus_small;
                        *p_x = currX + Q;
                        *p_y = currX + (two_times_modulus_small - Q);
                    }
                
                #ifdef USE_AVX2
                } else { // t >= 4
                    uint64_t *p_base = &operand[2 * i * t];
                    for (size_t j = 0; j < t; j+= 4) {
                        // The Harvey butterfly: assume X, Y in [0, 2p), and return
                        // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                        
                        uint64_t *p_x = &p_base[j];
                        uint64_t *p_y = &p_base[t + j];
                        __m256i x = _mm256_loadu_si256((__m256i const *)p_x);
                        __m256i y = _mm256_loadu_si256((__m256i const *)p_y);

                        // uint32_t currX = x - (two_times_modulus_small * static_cast<int>(x >= two_times_modulus_small));
                        __m256i cmp_val = _mm256_set1_epi64x(two_times_modulus_small);
                        __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val);
                        __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val);
                        __m256i currX = _mm256_sub_epi64(x, to_subtract);

                        // uint32_t Q = ((y) * (uint64_t)(Wprime)) >> 32;
                        __m256i Wprime_vec = _mm256_set1_epi64x(Wprime);
                        __m256i product = _mm256_mul_epu32(y, Wprime_vec);
                        __m256i Q = _mm256_srli_epi64(product, 32);

                        // Q = W * y - Q * modulus_small;
                        __m256i W_vec = _mm256_set1_epi64x(W);
                        __m256i WY = _mm256_mul_epu32(y, W_vec);
                        __m256i modulus_small_vec = _mm256_set1_epi64x(modulus_small);
                        __m256i Q_scaled = _mm256_mul_epu32(Q, modulus_small_vec);
                        __m256i Q_final = _mm256_sub_epi64(WY, Q_scaled);

                        __m256i new_x = _mm256_add_epi64(currX, Q_final);
                        __m256i Q_final_inverted = _mm256_sub_epi64(cmp_val, Q_final);
                        __m256i new_y = _mm256_add_epi64(currX, Q_final_inverted);

                        // *p_x = currX + Q;
                        // *p_y = currX + (two_times_modulus_small - Q);
                        _mm256_storeu_si256 ((__m256i *)p_x, new_x);
                        _mm256_storeu_si256 ((__m256i *)p_y, new_y);
                    }
                }
                #endif
            }
        }

        #ifdef USE_AVX2
            for (size_t i = 0; i < poly_len; i+=4) {
                    __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus_small);
                    __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
                    __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
                    __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                    x = _mm256_sub_epi64(x, to_subtract);

                    __m256i cmp_val2 = _mm256_set1_epi64x(modulus_small);
                    gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
                    to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                    x = _mm256_sub_epi64(x, to_subtract);
                    _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);
            }
        #else
            for (size_t i = 0; i < poly_len; i++) {
                operand[i] -= static_cast<int>(operand[i] >= two_times_modulus_small) * two_times_modulus_small;
                operand[i] -= static_cast<int>(operand[i] >= modulus_small) * modulus_small;
            }
        #endif
        
        // forward_tables = &tables[poly_len * 2 * 2 + poly_len * 2];
        // operand = &operand_overall[poly_len];
        // uint64_t modulus = a_i;
        // uint64_t two_times_modulus = 2 * modulus;

        // for (size_t mm = 0; mm < coeff_count_pow_of_2; mm++) {
        //     size_t m = 1 << mm;
        //     size_t t = n >> (mm+1);
        //     // m*t is always 1024

        //     for (size_t i = 0; i < m; i++) {
        //         const uint64_t W = forward_tables[m + i];
        //         const uint64_t Wprime = forward_tables[poly_len + m + i];

        //         for (size_t j = 0; j < t; j++) {
        //             // The Harvey butterfly: assume X, Y in [0, 2p), and return
        //             // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                    
        //             uint64_t *p_x = &operand[2 * i * t + j];
        //             uint64_t *p_y = &operand[2 * i * t + t + j];
        //             uint64_t x = *p_x;
        //             uint64_t y = *p_y;

        //             uint64_t currX = x - (two_times_modulus * static_cast<int>(x >= two_times_modulus));

        //             uint64_t Q = ((y) * (unsigned __int128)(Wprime)) >> 64;

        //             Q = W * y - Q * modulus;
        //             *p_x = currX + Q;
        //             *p_y = currX + (two_times_modulus - Q);
        //         }
        //     }
        // }

        // #if USE_AVX2
        //     for (size_t i = 0; i < poly_len; i+=4) {
        //         // operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
        //         // operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
        //         __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus);
        //         __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
        //         __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
        //         __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
        //         x = _mm256_sub_epi64(x, to_subtract);

        //         __m256i cmp_val2 = _mm256_set1_epi64x(modulus);
        //         gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
        //         to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
        //         x = _mm256_sub_epi64(x, to_subtract);
        //         _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);
        //     }
        // #else
        //     for (size_t i = 0; i < poly_len; i++) {
        //         operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
        //         operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
        //     }
        // #endif
    }
}

void ntt_inverse(uint64_t *operand_overall) {
    for (size_t coeff_mod = 0; coeff_mod < 2; coeff_mod++) {
        const uint64_t *inverse_tables = &tables[coeff_mod * poly_len * 2];
        uint64_t *operand = &operand_overall[coeff_mod == 0 ? 0 : poly_len];
        uint64_t modulus = coeff_mod == 0 ? crtq1 : crtq2;

        uint64_t two_times_modulus = 2 * modulus;
        size_t n = poly_len;
        size_t t = 1;

        for (size_t m = n; m > 1; m >>= 1) {
            size_t j1 = 0;
            size_t h = m >> 1;
            for (size_t i = 0; i < h; i++) {
                size_t j2 = j1 + t;
                // Need the powers of  phi^{-1} in bit-reversed order
                const uint64_t W = inverse_tables[h + i];
                const uint64_t Wprime = inverse_tables[poly_len + h + i];

                uint64_t *U = operand + j1;
                uint64_t *V = U + t;
                uint64_t currU;
                uint64_t T;
                uint64_t H;
                for (size_t j = j1; j < j2; j++) {
                    // U = x[i], V = x[i+m]

                    // Compute U - V + 2q
                    T = two_times_modulus - *V + *U;

                    // Cleverly check whether currU + currV >= two_times_modulus
                    currU = *U + *V - (two_times_modulus * static_cast<int>((*U << 1) >= T));

                    // Need to make it so that div2_uint_mod takes values that
                    // are > q.
                    // div2_uint_mod(U, modulusptr, coeff_uint64_count, U);
                    // We use also the fact that parity of currU is same as
                    // parity of T. Since our modulus is always so small that
                    // currU + masked_modulus < 2^64, we never need to worry
                    // about wrapping around when adding masked_modulus.
                    // uint64_t masked_modulus = modulus &
                    // static_cast<uint64_t>(-static_cast<int64_t>(T & 1));
                    // uint64_t carry = add_uint64(currU, masked_modulus, 0,
                    // &currU); currU += modulus &
                    // static_cast<uint64_t>(-static_cast<int64_t>(T & 1));
                    *U++ = (currU + (modulus * static_cast<int>(T & 1))) >> 1;

                    //multiply_uint64_hw64(Wprime, T, &H);
                    H = ((T) * (uint64_t)(Wprime)) >> 32;
                    
                    // effectively, the next two multiply perform multiply
                    // modulo beta = 2**wordsize.
                    *V++ = W * T - H * modulus;
                }
                j1 += (t << 1);
            }
            t <<= 1;
        }

        #ifdef USE_AVX2
            for (size_t i = 0; i < poly_len; i+=4) {
                // operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
                // operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;

                __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus);
                __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
                __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
                __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                x = _mm256_sub_epi64(x, to_subtract);

                __m256i cmp_val2 = _mm256_set1_epi64x(modulus);
                gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
                to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                x = _mm256_sub_epi64(x, to_subtract);
                _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);

                // operand[i] = operand[i] % modulus;
                // if (operand[i] >= two_times_modulus) operand[i] -=
                // two_times_modulus; if (operand[i] >= modulus) operand[i] -=
                // modulus;
            }
        #else
            for (size_t i = 0; i < poly_len; i++) {
                operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
                operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
            }
        #endif
    }
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

inline uint64_t barrett_coeff_q1(uint64_t val)
{
    // return val % moduli[n];
    return barrett_raw_u64(val, cr1_p, crtq1);
}

inline uint64_t barrett_coeff_q2(uint64_t val)
{
    return barrett_raw_u64(val, cr1_b, crtq2);
}

void computeForward(uint64_t* result, uint64_t* input)
{
    for (size_t n = 0; n < 2; n++)
    {
        for (size_t z = 0; z < poly_len; z++)
        {
            uint64_t val = input[z];
            result[poly_len * n + z] = barrett_coeff(val, n);//val % moduli[n];
        }
    }
    ntt_forward(result);
}

typedef struct __attribute__((packed, aligned(16))) {
    uint64_t x, y;
} ulonglong2_h;

inline ulonglong2_h umul64wide(uint64_t a, uint64_t b) {
    ulonglong2_h z;
    __uint128_t val = ((__uint128_t)a) * ((__uint128_t)b);
    z.x = (uint64_t)val;
    z.y = (uint64_t)(val >> 64);
    return z;
}

inline uint64_t cpu_add_u64(uint64_t operand1, uint64_t operand2,
                                uint64_t *result) {
    *result = operand1 + operand2;
    return (*result < operand1) ? 0x01 : 0x00; // detect overflow
}

inline uint64_t barrett_raw_u128(__uint128_t val, uint64_t const_ratio_0, uint64_t const_ratio_1, uint64_t modulus) {
    uint64_t zx = val & (((__uint128_t)1 << 64) - 1);
    uint64_t zy = val >> 64;

    uint64_t tmp1, tmp3, carry;
    ulonglong2_h prod = umul64wide(zx, const_ratio_0);
    carry = prod.y;
    ulonglong2_h tmp2 = umul64wide(zx, const_ratio_1);
    tmp3 = tmp2.y + cpu_add_u64(tmp2.x, carry, &tmp1);
    tmp2 = umul64wide(zy, const_ratio_0);
    carry = tmp2.y + cpu_add_u64(tmp1, tmp2.x, &tmp1);
    tmp1 = zy * const_ratio_1 + tmp3 + carry;
    tmp3 = zx - tmp1 * modulus;

    return tmp3;
}

inline uint64_t barrett_reduction_u128(__uint128_t val) {
    uint64_t reduced_val = barrett_raw_u128(val, cr0_Q, cr1_Q, crtMod);
    reduced_val -= (crtMod)*(static_cast<int>(reduced_val >= crtMod));
    return reduced_val;
}

uint64_t crt_compose(uint64_t x, uint64_t y) {
    // uint64_t val_ap = x * a_inv_p_i;
    // val_ap += y * p_inv_a_i;
    uint64_t val_ap_u64 = x;
    __uint128_t val = val_ap_u64 * b_inv_pa_i;
    val += y * pa_inv_b_i;
    // cout << "(" << x << ", " << y << ", " << z << ") -> " << (uint64_t)(val % Q_i) << endl;
    // return val % Q_i;
    return barrett_reduction_u128(val);
}

// uint64_t *scratch;
uint64_t* scratch = (uint64_t *)malloc(2 * poly_len * sizeof(uint64_t));

void computeInverse(uint64_t* result, uint64_t* input)
{
    // memcpy(scratch, &input, pol_bytes);
    memcpy(scratch, input, pol_bytes);
    ntt_inverse(scratch);
    for (size_t z = 0; z < poly_len; z++) {
        uint64_t v_x = scratch[z];
        uint64_t v_y = scratch[poly_len + z];
        // uint64_t v_z = scratch[2*coeff_count + z];
        uint64_t val = crt_compose(v_x, v_y);
        result[z] = val;
    }
}

constexpr size_t crt_count = 2;

void fastMultiplyQueryByDatabaseDim1(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per)
{
    // v_firstdim: 
    //     poly_len, dim0, ct.rows, ct.cols
    // db: 
    //     poly_len, num_per, pt.cols, dim0, pt.rows
    size_t base_dim = 2;

    size_t ct_rows = base_dim;
    size_t ct_cols = 1;
    size_t pt_rows = 1;
    size_t pt_cols = 1;

    // size_t count = 0;

    #if defined(__AVX512F__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                            // count++;
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[8];
                        alignas(64) uint64_t sums_out_n2_u64[8];
                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);

                        // TODO: check
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx] + sums_out_n0_u64[4 + idx];
                            // sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                            sums_out_n0_u64_acc[idx] += val;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx] + sums_out_n2_u64[4 + idx];
                            // sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                            sums_out_n2_u64_acc[idx] += val;
                        }
                    }
                    // TODO: check
                    /**
                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }
                    **/

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                }
            }
        }
        // std::cout << "origin count = " << count << std::endl;
    #elif defined(__AVX2__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m256i sums_out_n0 = _mm256_setzero_si256();
                        __m256i sums_out_n2 = _mm256_setzero_si256();

                        #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/4; i_jm++) {
                            size_t jm = o_jm * inner_limit + (4*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            __m256i b  = _mm256_set_epi64x(b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m256i a = _mm256_load_si256((const __m256i *)v_a);
                            __m256i a_lo = a;
                            __m256i a_hi_hi = _mm256_srli_epi64(a, packed_offset_2);
                            __m256i b_lo = b;
                            __m256i b_hi_hi = _mm256_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm256_add_epi64(sums_out_n0, _mm256_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm256_add_epi64(sums_out_n2, _mm256_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[4];
                        alignas(64) uint64_t sums_out_n2_u64[4];
                        _mm256_store_si256((__m256i *)sums_out_n0_u64, sums_out_n0);
                        _mm256_store_si256((__m256i *)sums_out_n2_u64, sums_out_n2);

                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx];
                            sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx];
                            sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                        }
                    }

                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                }
            }
        }
    #else
        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    // size_t idx_b = idx_b_base + i * (pt_cols * dim0 * n0) + c * (dim0 * n0);
                    
                    __uint128_t sums_out_n0_0 = 0, sums_out_n0_1 = 0;
                    __uint128_t sums_out_n1_0 = 0, sums_out_n1_1 = 0;

                    // #pragma unroll 16
                    for (size_t jm = 0; jm < dim0 * pt_rows; jm++) { // 0... (128 * 1) 
                        uint64_t b = db[idx_b_base++];
                        
                        const uint64_t *v_a = &v_firstdim[idx_a_base + jm * ct_rows]; // ct_row = 2;
                        uint64_t v_a0 = v_a[0];
                        uint64_t v_a1 = v_a[1];

                        uint32_t b_lo = b;
                        uint32_t b_hi = b >> 32L;

                        uint32_t v_a0_lo = v_a0;
                        uint32_t v_a0_hi = v_a0 >> 32L;

                        uint32_t v_a1_lo = v_a1;
                        uint32_t v_a1_hi = v_a1 >> 32L;
                        
                        // do n0
                        sums_out_n0_0 += (v_a0_lo) * (uint64_t)b_lo;
                        sums_out_n0_1 += (v_a1_lo) * (uint64_t)b_lo;

                        // do n1
                        sums_out_n1_0 += ((uint64_t)v_a0_hi) * b_hi;
                        sums_out_n1_1 += ((uint64_t)v_a1_hi) * b_hi;
                    }

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z; // c = 0 always
                    out[i][idx_c] = sums_out_n0_0 % crtq1;
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = sums_out_n0_1 % crtq1;

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = sums_out_n1_0 % crtq2;
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = sums_out_n1_1 % crtq2;

                    // out[i]: cta, cta, ctb, ctb
                }
            }
        }
    #endif
}


void fastMultiplyQueryByDatabaseDim1InvCRT(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per)
{
    // v_firstdim: 
    //     poly_len, dim0, ct.rows, ct.cols
    // db: 
    //     poly_len, num_per, pt.cols, dim0, pt.rows
    size_t base_dim = 2;

    size_t ct_rows = base_dim;
    size_t ct_cols = 1;
    size_t pt_rows = 1;
    size_t pt_cols = 1;

    #if defined(__AVX512F__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            // auto start  = std::chrono::high_resolution_clock::now();
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                // for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[8];
                        alignas(64) uint64_t sums_out_n2_u64[8];
                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);

                        // TODO: check
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx] + sums_out_n0_u64[4 + idx];
                            // sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                            sums_out_n0_u64_acc[idx] += val;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx] + sums_out_n2_u64[4 + idx];
                            // sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                            sums_out_n2_u64_acc[idx] += val;
                        }
                    }
                    // TODO: check
                    /**
                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }
                    **/

                    /**
                    // output n0
                    size_t n = 0;
                    size_t idx_c =  n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);

                    idx_c += (crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);


                    // output n1
                    n = 1;
                    idx_c = n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);

                    idx_c += (crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                    **/

                    /**
                    uint64_t temp1 = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    uint64_t temp2 = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);
                    uint64_t temp3 = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    uint64_t temp4 = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);

                    out[i][z] = crt_compose(temp1, temp3, 0);
                    out[i][poly_len + z] = crt_compose(temp2, temp4, 0);;
                    **/
                    out[i][z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i][poly_len + z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                // }
            }
            // auto stop  = std::chrono::high_resolution_clock::now();

            // auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            // std::cout << " (in total) database multiplication total costs " << glapsed.count() << " us." << std::endl;
        }
        #elif defined(__AVX2__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m256i sums_out_n0 = _mm256_setzero_si256();
                        __m256i sums_out_n2 = _mm256_setzero_si256();

                        #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/4; i_jm++) {
                            size_t jm = o_jm * inner_limit + (4*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            __m256i b  = _mm256_set_epi64x(b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m256i a = _mm256_load_si256((const __m256i *)v_a);
                            __m256i a_lo = a;
                            __m256i a_hi_hi = _mm256_srli_epi64(a, packed_offset_2);
                            __m256i b_lo = b;
                            __m256i b_hi_hi = _mm256_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm256_add_epi64(sums_out_n0, _mm256_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm256_add_epi64(sums_out_n2, _mm256_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[4];
                        alignas(64) uint64_t sums_out_n2_u64[4];
                        _mm256_store_si256((__m256i *)sums_out_n0_u64, sums_out_n0);
                        _mm256_store_si256((__m256i *)sums_out_n2_u64, sums_out_n2);

                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx];
                            // sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                            sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val);
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx];
                            // sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                            sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val);
                        }
                    }
                /**
                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }
                **/
                /**
                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                **/ 
                    out[i][z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i][poly_len + z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                }
            }
        }
    #endif
}

// encode to crt form
void database_tocrt(uint64_t* datacrt, std::vector<std::vector<uint64_t> >& data_ntt, int32_t N1)
{
    int32_t N2 = N / 2 / N1;
    
    /** datacrt = [a, e, i, m, b, f, j, \cdots] 
    data_ntt = [a, b, c, d
                e, f, g, h
                i, j, k, l
                m, n, o, p]
    **/

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N2; j++)
        {
            for (size_t k = 0; k < N1; k++)
            {
                uint64_t temp = data_ntt[j * N1 + k][i];
                // TODO: check
                datacrt[i * N1 * N2 + j * N1 + k] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
                // uint64_t temp1 = barrett_coeff(temp, 0);
                // uint64_t temp2 = barrett_coeff(temp, 1);
                // datacrt[i * N1 * N2 + j * N1 + k] = temp1 | temp2 << 32;
            }
        }
    }
}
