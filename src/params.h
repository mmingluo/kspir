#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>

#define INTEL_HEXL

#define SECURITY
#define USE_CRTMOD

#define THREADS

#ifndef SECURITY
#define N 256 // 256
#else
#define N 4096
// #define N 2048
#endif

#ifdef THREADS
#define THREADS_NUM 16
#endif

// #define DATA_INDEPENDENT_HANDLE

#define mod 8380417 // (UINT64_C(0x01) << 32), omega = 1753
//#define bigMod 4398046486529 // (UINT64_C(0x01) << 42), omega = 23254882976, when N = 256
// #define bigMod 281474976694273 // (UINT64_C(0x01) << 42), omega = 384399401, when N = 2048

#define bigMod 1125899906826241 // pow(2, 50)
// #define bigMod 18014398509450241 // pow(2, 54)

#define bigMod2 4398046504961 //TODO


#define mod2inv 569391470430536 // 149000519487509 // q_2^-1 (mod q_1)
#define mod1inv 2173861076666 // 2069913391085 // q_1^-1 (mod q_2)


#if bigMod == 1125899906826241
    #define Delta 17179869183 // floor(q/p)
    #define Bg 0x01 << 17

    // two128_mod = 4874087337207107259134443424
    #define two128_mod_hignbits 528449607 // 2^128 (mod q1 * q2) >> 63
    #define two128_mod_lowbits 9116411467323735968

    #if N == 4096
    #define bNinv 1125625028919301
    #elif N == 256
    #define bNinv 1121501860315201
    #elif N == 2048
    #define bNinv 1125350151012361
    #else
    #error "Please add more parameters."
    #endif
#elif bigMod == 18014398509450241
    #define Delta 1125899906826241 // floor(q/p)
    #define Bg 0x01 << 18 

    #if N == 4096
    #define bNinv 4499201580851464
    #elif N == 256
    #define bNinv 17944029765272701
    #elif N == 2048
    #define bNinv 8998403161702928
    #else
    #error "Please add more parameters."
    #endif
#endif

#if bigMod2 == 4398046504961
    #if N == 4096
    #define bNinv2 2197949510658
    #elif N == 256
    #define bNinv2 4380866635801
    #elif N == 2048
    #define bNinv2 4395899021316
    #else
    #error "Please add more parameters."
    #endif
#endif


#define Pbits 16

#define ell 3 // 3
#define Base 0 // have to set to zero



// bsgs
#if N == 4096
    #define bsgsp 65537 // 40961 // 40961 % 8192 = 1
    #ifndef USE_CRTMOD
    #define bsgsDelta 27487119621
    #else
    #define bsgsDelta 1021937069741 // 1635084342169 for 40961
    #endif
#elif N == 2048
    #define bsgsp 40961 // 40961 % 8192 = 1
    #ifndef USE_CRTMOD
    #define bsgsDelta 27487119621
    #else
    #define bsgsDelta 1635084342169
    #endif
#elif N == 256
    #define bsgsp 7681  // 7681  % 512  = 1
    #ifndef USE_CRTMOD
    #define bsgsDelta 146582464109
    #else
    #define bsgsDelta 8719527371384
    #endif
#endif



/*
#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
    static intel::hexl::NTT ntts(N, mod);
    static intel::hexl::NTT nttb(N, bigMod);
#endif
*/

#endif
