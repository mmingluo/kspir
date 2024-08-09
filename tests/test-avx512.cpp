#include <stddef.h>
#include <stdint.h>
#include <vector>
#include <chrono>
#include <iostream>

#include "pir.h"

#if defined(__AVX512F__)
#include <immintrin.h>
#include <avx512fintrin.h>
#endif

void test_avx2()
{
#if defined(__AVX2__)
    std::cout << "have avx2" << std::endl;
#else
    std::cout << "error: no avx2" << std::endl;
#endif
}


void databaseMultiply(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    const size_t N1, const size_t num_per,
    size_t inner_limit = 512)
{
    size_t outer_limit = 1;
    // size_t inner_limit = 512;
    // size_t inner_limit = 1024;

    size_t packed_offset_2 = 32;
    size_t count = 0;

    #if defined(__AVX512F__)
                    // uint64_t sums_out_n0_u64_acc = 0;
                    // uint64_t sums_out_n2_u64_acc = 0;

                    // reduce here, otherwise we will overflow
                    alignas(64) uint64_t sums_out_n0_u64[8];
                    alignas(64) uint64_t sums_out_n2_u64[8];

                    size_t idx_a_base = 0;
                    size_t idx_b_base = 0;
    for (size_t z = 0; z < N; z++)
    {
        idx_a_base = z * N1;
        for (size_t i = 0; i < num_per; i++) {
                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

    // auto start = std::chrono::high_resolution_clock::now();
                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            // __m512i b = _mm512_load_si512(&db[idx_b_base + jm]);
                            // __m512i b = _mm512_load_si512(&db[idx_b_base]);
                            // idx_b_base += 8;
                            /**/
                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            uint64_t b_inp_5 = db[idx_b_base++]; // online
                            uint64_t b_inp_6 = db[idx_b_base++];
                            uint64_t b_inp_7 = db[idx_b_base++];
                            uint64_t b_inp_8 = db[idx_b_base++];
                            // __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);
                            __m512i b  = _mm512_set_epi64(b_inp_8, b_inp_7, b_inp_6, b_inp_5, b_inp_4, b_inp_3, b_inp_2, b_inp_1);
                            /**/
                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                            // sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_add_epi32(a_lo, b_lo));
                            // sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_add_epi32(a_hi_hi, b_hi_hi));
                            count++;
                        }
    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "Time (in) taken by AVX512: " << duration.count() << " microseconds" << std::endl;                        

                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);
                    }

                    for (size_t index = 0; index < 8; ++index) {
                        out[i][z] = sums_out_n0_u64[index] + sums_out_n2_u64[index];
                    }
        }
    }
    std::cout << "count: " << count << std::endl;
    #endif
}

void databaseMultiply2(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    const size_t N1, const size_t num_per,
    size_t inner_limit = 512)
{
    size_t outer_limit = 1;
    // size_t inner_limit = 512;
    // size_t inner_limit = 1024;

    size_t packed_offset_2 = 32;
    size_t count = 0;

    #if defined(__AVX512F__)
                    // uint64_t sums_out_n0_u64_acc = 0;
                    // uint64_t sums_out_n2_u64_acc = 0;

                    // reduce here, otherwise we will overflow
                    alignas(64) uint64_t sums_out_n0_u64[8];
                    alignas(64) uint64_t sums_out_n2_u64[8];

                    size_t idx_a_base = 0;
                    size_t idx_b_base = 0;
    for (size_t z = 0; z < N; z++)
    {
        idx_a_base = z * N1;
        for (size_t i = 0; i < num_per; i++) {
                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

    // auto start = std::chrono::high_resolution_clock::now();
                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            // uint64_t b_inp_5 = db[idx_b_base++]; // online
                            // uint64_t b_inp_6 = db[idx_b_base++];
                            // uint64_t b_inp_7 = db[idx_b_base++];
                            // uint64_t b_inp_8 = db[idx_b_base++];
                            __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);
                            // __m512i b  = _mm512_set_epi64(b_inp_8, b_inp_7, b_inp_6, b_inp_5, b_inp_4, b_inp_3, b_inp_2, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                            // sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_add_epi32(a_lo, b_lo));
                            // sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_add_epi32(a_hi_hi, b_hi_hi));
                            count++;
                        }
    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "Time (in) taken by AVX512: " << duration.count() << " microseconds" << std::endl;                        

                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);
                    }

                    for (size_t index = 0; index < 8; ++index) {
                        out[i][z] = sums_out_n0_u64[index] + sums_out_n2_u64[index];
                    }
        }
    }
    std::cout << "count: " << count << std::endl;
    #endif
}

// same as databaseMultiply, but with two same entries
void databaseMultiply3(
    std::vector<std::vector<uint64_t> >& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    const size_t N1, const size_t num_per,
    size_t inner_limit = 512)
{
    size_t outer_limit = 1;
    // size_t inner_limit = 512;
    // size_t inner_limit = 1024;

    size_t packed_offset_2 = 32;
    size_t count = 0;

    #if defined(__AVX512F__)
                    // uint64_t sums_out_n0_u64_acc = 0;
                    // uint64_t sums_out_n2_u64_acc = 0;

                    // reduce here, otherwise we will overflow
                    alignas(64) uint64_t sums_out_n0_u64[8];
                    alignas(64) uint64_t sums_out_n2_u64[8];

                    size_t idx_a_base = 0;
                    size_t idx_b_base = 0;
    for (size_t z = 0; z < N; z++)
    {
        idx_a_base = z * N1;
        idx_b_base = z * N1 * num_per;
        for (size_t i = 0; i < num_per; i++) {
                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

    // auto start = std::chrono::high_resolution_clock::now();
                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            // __m512i b = _mm512_load_si512(&db[idx_b_base + jm]);
                            // __m512i b = _mm512_load_si512(&db[idx_b_base + i * N1 + jm]);
                            // __m512i b = _mm512_load_si512(&db[idx_b_base]);
                            // idx_b_base += 8;
                            /**/
                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            // idx_b_base += 4;
                            uint64_t b_inp_5 = db[idx_b_base++]; // online
                            uint64_t b_inp_6 = db[idx_b_base++];
                            uint64_t b_inp_7 = db[idx_b_base++];
                            uint64_t b_inp_8 = db[idx_b_base++];
                            // __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);
                            // __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_2, b_inp_3, b_inp_1, b_inp_2, b_inp_4, b_inp_1, b_inp_3);
                            __m512i b  = _mm512_set_epi64(b_inp_8, b_inp_7, b_inp_6, b_inp_5, b_inp_4, b_inp_3, b_inp_2, b_inp_1);
                            /**/
                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                            // sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_add_epi32(a_lo, b_lo));
                            // sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_add_epi32(a_hi_hi, b_hi_hi));
                            count++;
                        }
    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "Time (in) taken by AVX512: " << duration.count() << " microseconds" << std::endl;                        

                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);
                    }
                    /**
                    for (size_t index = 0; index < 8; ++index) {
                       out[i][z] = sums_out_n0_u64[index] + sums_out_n2_u64[index];
                    }
                    **/
                    /**/
                    out[i][z] = crt_compose(sums_out_n0_u64[0]+sums_out_n0_u64[1]+
                            sums_out_n0_u64[2]+sums_out_n0_u64[3]+
                            sums_out_n0_u64[4]+sums_out_n0_u64[5]+
                            sums_out_n0_u64[6]+sums_out_n0_u64[7],
                            sums_out_n2_u64[0]+sums_out_n2_u64[1]+
                            sums_out_n2_u64[2]+sums_out_n2_u64[3]+
                            sums_out_n2_u64[4]+sums_out_n2_u64[5]+
                            sums_out_n2_u64[6]+sums_out_n2_u64[7]);
                    /**/
        }
    }
    std::cout << "count: " << count << std::endl;
    #endif
}

void test_avx512()
{
    size_t N1 = 512;
    size_t N2 = N / 2 / N1;

    std::vector<std::vector<uint64_t> > out(N2, std::vector<uint64_t>(N, 0));

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N);
    uint64_t* database = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N2 * N);
    for (size_t i = 0; i < N1 * N; ++i) {
        uint64_t temp1 = (((uint64_t)rand() & 0xfffffff) << 32) | ((uint64_t)rand() & 0xfffffff);
        cipherbuf[i] = temp1;
    }
    for (size_t i = 0; i < N1 * N2 * N; ++i) {
        uint64_t temp2 = (((uint64_t)rand() & 0xfffffff) << 32) | ((uint64_t)rand() & 0xfffffff);
        database[i] = temp2;
    }

int32_t ntimes = 1;
    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; ++i)
{
    databaseMultiply3(out, database, cipherbuf, N1, N2);
    // databaseMultiply(out, database, cipherbuf, N1, N2, 256);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by AVX512: " << duration.count() << " microseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; ++i)
{
    // databaseMultiply(out, database, cipherbuf, N1, N2, 256);
    databaseMultiply2(out, database, cipherbuf, N1, N2);
}
    stop = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time (half) taken by AVX512: " << duration.count() << " microseconds" << std::endl;

}

void test_multidatabase_compare()
{
    std::vector<std::vector<uint64_t> > out(1, std::vector<uint64_t>(2, 0));

    size_t N1 = 512;
    // size_t N2 = N / 2 / N1;

    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N);
    uint64_t* database = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * N1 * N);
    for (size_t i = 0; i < N1 * N; ++i) {
        uint64_t temp1 = (((uint64_t)rand() << 32) & 0xfffffff) || (rand() & 0xffffffff);
        uint64_t temp2 = (((uint64_t)rand() << 32) & 0xfffffff) || (rand() & 0xffffffff);
        cipherbuf[i] = temp1;
        database[i] = temp2;
    }

    auto start = std::chrono::high_resolution_clock::now();
    // databaseMultiply(out, database, cipherbuf);
    // fastMultiplyQueryByDatabaseDim1InvCRTOnline(output, database, cipherbuf_on, N1, N2);

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by AVX512: " << duration.count() << " microseconds" << std::endl;   
}

int main(int argc, char** argv)
{
    // test_avx2();

    test_avx512();

    return 0;
}