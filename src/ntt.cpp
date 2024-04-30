#include <stdint.h>
#include <vector>

#include "ntt.h"




#define MONT -4186625 // 2^32 % Q
#define QINV 58728449 // q^(-1) mod 2^32
#define Q mod // dilithium

#define bigMOUT 103075020800 // 2^64 % bigQ
#define bigQINV 148629233667694593 // bigQ^(-1) mod 2^64
#define bigQ bigMod // 4398046486529

int32_t montgomery_reduce(int64_t a) {
  int32_t t;

  t = (int32_t)a*QINV;
  t = (a - (int64_t)t*Q) >> 32;
  return t;
}

int64_t montgomery_reduce(__int128_t a) {
  int64_t t;

  t = (int64_t)a*bigQINV;
  t = (a - (__int128_t)t*bigQ) >> 64;
  return t;
}

int64_t barret_reduce(int64_t a) {
    int64_t t;
    
    // v = round(2^(40+64) /bigQ.);
    // const int64_t v = 4611686044196143104;
    const int64_t v = 4611686044196143248;
    t = (__int128_t)v*a >> 104;

    t *= bigQ;
    return a - t;
}

int64_t barret_mult_reduce(__uint128_t a) {
    int64_t t;
    
    // v = round(2^(40+64) /bigQ.);s
    const int64_t v = 4611686044196143248;
    t = (__int128_t)v*a >> 104;

    t *= bigQ;
    return a - t;
}


/*
// just tests
void ntt(uint64_t a[dim])
{

}
*/

void ntt(int32_t a[dim]) {
  unsigned int len, start, j, k;
  int32_t zeta, t;

  k = 0;
  for(len = dim/2; len > 0; len >>= 1) {
    for(start = 0; start < dim; start = j + len) {
      zeta = zetas[++k];
      for(j = start; j < start + len; ++j) {
        t = montgomery_reduce((int64_t)zeta * a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

void ntt(int64_t a[dim]) {
  unsigned int len, start, j, k;
  int64_t zeta, t;

  k = 0;
  // for(len = 128; len > 0; len >>= 1) {
  for(len = dim/2; len > 0; len >>= 1) {
    for(start = 0; start < dim; start = j + len) {
      zeta = bigZetas[++k];
      for(j = start; j < start + len; ++j) {
        t = montgomery_reduce((__int128_t)zeta * a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

void invntt_tomont(int32_t a[dim]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
#ifndef SECURITY
  const int32_t f = 41978; // mont^2/256
#else
  const int32_t f = 6290560; // mont^2/2048
#endif

  k = dim;
  for(len = 1; len < dim; len <<= 1) {
    for(start = 0; start < dim; start = j + len) {
      zeta = -zetas[--k];
      for(j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = montgomery_reduce((int64_t)zeta * a[j + len]);
      }
    }
  }

  for(j = 0; j < dim; ++j) {
    a[j] = montgomery_reduce((int64_t)f * a[j]);
  }
}

void invntt_tomont(int64_t a[dim]) {
  unsigned int start, len, j, k;
  int64_t t, zeta;
#ifndef SECURITY
  const int64_t f = 300619399936; // mont^2/256
#else
  const int64_t f = 37577424992; // mont^2/2048
#endif

  // k = 256;
  k = dim;
  for(len = 1; len < dim; len <<= 1) {
    for(start = 0; start < dim; start = j + len) {
      zeta = -bigZetas[--k];
      for(j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = montgomery_reduce((__int128_t)zeta * a[j + len]);
      }
    }
  }

  for(j = 0; j < dim; ++j) {
    a[j] = montgomery_reduce((__int128_t)f * a[j]);
  }
}

void hadamard_mult(int32_t result[dim], int32_t a[dim], int32_t s[dim])
{
    for (unsigned int i = 0; i < dim; i++)
    {
      // TODO:
      // result[i] = (a[i] * s[i]) % Q;
      result[i] = montgomery_reduce((int64_t)a[i] * s[i]);
    }
}

void hadamard_mult(int64_t result[dim], int64_t a[dim], int64_t s[dim])
{
    for (unsigned int i = 0; i < dim; i++)
    {
      // TODO:
      // result[i] = (a[i] * s[i]) % Q;
      result[i] = montgomery_reduce((__int128_t)a[i] * s[i]);
    }
}

void hadamard_mult(int32_t result[dim], int32_t a[dim], std::vector<uint64_t> s)
{
    for (unsigned int i = 0; i < dim; i++)
    {
      // TODO:
      // result[i] = (a[i] * s[i]) % Q;
      result[i] = montgomery_reduce((int64_t)a[i] * (int64_t)s[i]);
    }
}

void hadamard_mult(int64_t result[dim], int64_t a[dim], std::vector<uint64_t> s)
{
    for (unsigned int i = 0; i < dim; i++)
    {
      // TODO:
      // result[i] = (a[i] * s[i]) % Q;
      result[i] = montgomery_reduce((__int128_t)a[i] * s[i]);
    }
}

// 0 <= a[i],b[i] < modulus
void subtraction(uint64_t result[dim], uint64_t a[dim], uint64_t b[dim], uint64_t modulus)
{
  for (size_t i = 0; i < dim; i++)
  {
    int64_t sub = a[i] - b[i];
    result[i] = sub < 0 ? sub += modulus : sub; 
  }
}
