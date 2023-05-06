#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"


// 256 MB / 8
void test_two_dimensions()
{   
    Secret queryKey(bigMod), answerKey(bigMod);
    uint64_t row = 123, col = 123;

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N, 0));
    sample_database(data);

    // for (size_t i = 0 ; i < N; i++)
    //    data[i][col] = 0;
    std::cout << "the wanted message is " << data[row][col] << std::endl;
    
    database_to_signed(data, Pbits, bigMod);
    database_tontt(data);

    auto start_qu = std::chrono::high_resolution_clock::now();
    RlweCiphertext query1(N, bigMod);
    query(query1, queryKey, col);
    RGSWCiphertext query2(N, bigMod);
    query2.keyGen(answerKey, row, true); // note there is reverses
    auto stop_qu = std::chrono::high_resolution_clock::now();
    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query time costs " << glapsed_qu.count() << " us" << std::endl;

    // build keyswitch key
    std::vector<std::vector<uint64_t> > ks(2 * ell * N, std::vector<uint64_t>(N, 0));
    dummy_ksKeyGen(ks, queryKey, answerKey, data);

    RlweCiphertext kskHint(N, bigMod);
    RlweCiphertext result(N, bigMod);

int ntimes = 1;

    // answer
    auto start_ksk = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    ksKey_hint(kskHint.b, kskHint.a, query1.a, ks);
}
    // ksKey_hint_variant(kskHint.b, kskHint.a, query1.a, ks);
    auto stop_ksk = std::chrono::high_resolution_clock::now();
    auto glapsed_ksk = std::chrono::duration_cast<std::chrono::microseconds>(stop_ksk - start_ksk);
    std::cout << " ksk hint time costs " << glapsed_ksk.count() << " us" << std::endl;

    auto start_ans = std::chrono::high_resolution_clock::now();
    // online answer
for (size_t i = 0; i < ntimes; i++)
{
    answer_two_dimensions(result, kskHint.b, kskHint.a, query1.b, data, query2);
}
    auto stop_ans = std::chrono::high_resolution_clock::now();
    auto glapsed_ans = std::chrono::duration_cast<std::chrono::microseconds>(stop_ans - start_ans);
    std::cout << " answer time costs " << glapsed_ans.count() << " us" << std::endl;

    // recover
    auto start_rec = std::chrono::high_resolution_clock::now();
    std::vector<uint64_t> decrypted(N);
    recover(decrypted, result, answerKey);
    auto stop_rec = std::chrono::high_resolution_clock::now();
    auto glapsed_rec = std::chrono::duration_cast<std::chrono::microseconds>(stop_rec - start_rec);
    std::cout << " recover time costs " << glapsed_rec.count() << " us" << std::endl;


    std::cout << "the recovered value is " << decrypted[0] << std::endl;
}


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_two_dimensions();

    return 0;
}
