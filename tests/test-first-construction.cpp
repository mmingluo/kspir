#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"

// we test the first construction

// 32 MB
void test_first_dimension()
{   
    Secret queryKey(bigMod), answerKey(bigMod);

    uint64_t row = rand() & (N - 1);
    uint64_t col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N, 0));
    sample_database(data);

    // for (size_t i = 0; i < N; i++)
    //    data[i][col] = 0;
    std::cout << "the wanted message is " << data[row][col] << std::endl;

    database_to_signed(data, Pbits, bigMod);
    database_tontt(data);

    auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, bigMod);
    query(query1, queryKey, col);

    auto stop_qu = std::chrono::high_resolution_clock::now();
    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << "query time costs " << glapsed_qu.count() << " us" << std::endl;

    // build keyswitch key
    std::vector<std::vector<uint64_t> > pk(2 * ell * N, std::vector<uint64_t>(N, 0));
    pkKeyGen(pk, queryKey, answerKey);

    // RlweCiphertext kskHint(N, bigMod);
    RlweCiphertext result(N, bigMod);

    // answer
    auto start_ksk = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < 1; i++)
{
    ksKey_hint_first_construction(result.b, result.a, query1.a, pk, data);
}
    // ksKey_hint_variant(result.b, result.a, query1.a, ks);
    auto stop_ksk = std::chrono::high_resolution_clock::now();

    auto start_ans = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < 1; i++)
{
    // answer_first_dimension(result.b, result.a, query1.b, query1.a, data);
    answer_first_dimension(result.b, result.a, query1.b, data);
}
    auto stop_ans = std::chrono::high_resolution_clock::now();


    auto glapsed_ksk = std::chrono::duration_cast<std::chrono::microseconds>(stop_ksk - start_ksk);
    auto glapsed_ans = std::chrono::duration_cast<std::chrono::microseconds>(stop_ans - start_ans);
    std::cout << "ksk hint time costs " << glapsed_ksk.count() << " us" << std::endl;
    std::cout << "answer time costs " << glapsed_ans.count() << " us" << std::endl;

    // recover
    auto start_re = std::chrono::high_resolution_clock::now();
    std::vector<uint64_t> decrypted(N);
    recover(decrypted, result, answerKey);

    auto stop_re = std::chrono::high_resolution_clock::now();
    auto glapsed_re = std::chrono::duration_cast<std::chrono::microseconds>(stop_re - start_re);
    std::cout << "recover costs " << glapsed_re.count() << " us" << std::endl;

    std::cout << "the recovered value is " << decrypted[row] << std::endl;
}


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_first_dimension();

    return 0;
}
