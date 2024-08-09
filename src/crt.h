#ifndef CRT_H
#define CRT_H

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "params.h"

/**
 * we use another ntt and multiplication implementation that is same as Spiral.
 * 
*/


#define crtq1 268369921 // 2^28 - 2^16 + 1
#define crtq2 249561089UL // 266334209; // 2^28 - 2^21 - 2^12 + 1
#define crtMod (crtq1 * crtq2)

constexpr uint64_t cr1_p = 68736257792UL;
constexpr uint64_t cr1_b = 73916747789UL;

constexpr size_t max_summed_pa_or_b_in_u64 = 1 << 6;
constexpr uint64_t packed_offset_1 = 32;
constexpr uint64_t packed_offset_diff = 0; // log (a_i)
constexpr uint64_t packed_offset_2 = packed_offset_1 + packed_offset_diff;

// mini root of unity for intel::hexl::NTT
constexpr uint64_t mini_root_of_unity = 10297991595;


#define bsMod 16760833 // 2^24 - 2^14 + 1
#define bsModInv 51309461009346113UL // bsMod^-1 (mod crtMod)
#define crtModInv 3920321 // crtMod^-1 (mod bsMod)

// #define auxMod 16736257UL // auxiliary modulus for offline phase, 24 bits here
#define auxMod 268361729UL // 28 bits here
#if auxMod == 16736257UL
    #define crtBaMod (bsMod * auxMod) // bsMod * auxMod

    #define crtBaModInv 16509166641422680UL // crtBaMod^-1 (mod crtMod)
    #define crtModInvBa 211367399390890UL // crtMod^-1 (mod crtBaMod)

    #define auxModInv 47918398092285591UL // auxMod^-1 (mod crtMod)
    #define auxModInvBs 682 // auxMod^-1 (mod bsMod)

    // two128_mod = 2^128 % (crtMod * crtBaMod) = 16922090791964026345071048261795
    #define two128_mod_hignbits_bsaux 1834696759964UL
    #define two128_mod_lowbits_bsaux 4009972933838110883ULL
#elif auxMod == 268361729UL
    #define crtBaMod (bsMod * auxMod) // bsMod * auxMod

    #define crtBaModInv 26456327859403970UL // crtBaMod^-1 (mod crtMod)
    #define crtModInvBa 2721180491014976UL // crtMod^-1 (mod crtBaMod)

    #define auxModInv 48735363086513670UL // auxMod^-1 (mod crtMod)
    #define auxModInvBs 4280662 // auxMod^-1 (mod bsMod)

    // two128_mod = 2^128 % (crtMod * crtBaMod) = 136175933503495683303827764855712
    #define two128_mod_hignbits_bsaux 14764224294473UL
    #define two128_mod_lowbits_bsaux 1483995998376346528ULL

    #define crtModInvAux 73083388
    // two128_mod = 2^128 % (crtMod * auxMod) = 12725562273183741650962405
    #define two128_mod_hignbits_aux 1379708
    #define two128_mod_lowbits_aux 2086958912630458341ULL
#else
#error "Please add more parameters."
#endif

// two128_mod = 2^128 % (crtMod * bsMod) = 505663912589501352409913
#define two128_mod_hignbits_bs 54824 // 2^128 (mod q1 * q2) >> 63
#define two128_mod_lowbits_bs 1764040975123512121ULL // 2^128 (mod q1 * q2) & (2^63 - 1)



#if N == 4096
    constexpr size_t coeff_count_pow_of_2 = 12;
    constexpr size_t poly_len = 0x01 << coeff_count_pow_of_2;

    constexpr uint64_t root_of_unity_crt = 3375402822066082UL;
    constexpr uint64_t root_of_unity_crt_bamod = 1486445687605966UL;
#elif N == 2048
    constexpr size_t coeff_count_pow_of_2 = 11;
    constexpr size_t poly_len = 0x01 << coeff_count_pow_of_2;

    constexpr uint64_t root_of_unity_crt = 38878761190133527UL;
#elif N == 256
    constexpr size_t coeff_count_pow_of_2 = 8;
    constexpr size_t poly_len = 256;

    constexpr uint64_t root_of_unity_crt = 46801507955698810UL;
#else
    #error "Please add more parameters."
#endif




/**
 * @brief computeForward computes ntt form
 * @param result the length is 2 * N
 * @param input the length is N
 * 
*/
void computeForward(uint64_t* result, uint64_t* input);


/**
 * @brief computeInverse computes ntt form
 * @param result the length is N
 * @param input the length is 2 * N
 * 
*/
void computeInverse(uint64_t* result, uint64_t* input);

inline uint64_t barrett_coeff(uint64_t val, size_t n);

uint64_t crt_compose(uint64_t x, uint64_t y);

void fastMultiplyQueryByDatabaseDim1(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per);

void fastMultiplyQueryByDatabaseDim1InvCRT(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per);

void database_tocrt(uint64_t* datacrt, std::vector<std::vector<uint64_t> >& data_ntt, int32_t N1);

#endif
