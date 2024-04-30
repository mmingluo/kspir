#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"


void database_init(std::vector<std::vector<uint64_t> >& data, uint64_t row = 123, uint64_t col = 123)
{
    sample_database(data);

    std::cout << "the wanted message is " << data[row][col] << std::endl;
    database_tontt(data);
}


// TODO: finish this function
void test_larger_database()
{   
    Secret queryKey(bigMod), answerKey(bigMod);

    uint64_t r = (0x01 << 3); // 8

    // sample database
    std::vector<std::vector<uint64_t> > data(r * N, std::vector<uint64_t>(N, 0));

    uint64_t row = 123, col = 123;

    database_init(data, row, col)

    build_query()
    RlweCiphertext query1(N, bigMod);
    query(query1, queryKey, col);
    RGSWCiphertext query2(N, bigMod);
    query2.keyGen(answerKey, row, true); // note there is reverse

    // build keyswitch key
    std::vector<std::vector<uint64_t> > ks(2 * ell * N, std::vector<uint64_t>(N, 0));
    dummy_ksKeyGen(ks, queryKey, answerKey, data);

    RlweCiphertext kskHint(N, bigMod);
    RlweCiphertext result(N, bigMod);

    // answer
    auto start_ksk = std::chrono::high_resolution_clock::now();
    ksKey_hint(kskHint.b, kskHint.a, query1.a, ks);
    auto stop_ksk = std::chrono::high_resolution_clock::now();
    auto glapsed_ksk = std::chrono::duration_cast<std::chrono::milliseconds>(stop_ksk - start_ksk);
    std::cout << " ksk hint time costs " << glapsed_ksk.count() << " ms" << std::endl;

    auto start_ans = std::chrono::high_resolution_clock::now();
    // TODO: change result
    answer_two_dimensions(result, kskHint.b, kskHint.a, query1.b, data, query2);

    auto stop_ans = std::chrono::high_resolution_clock::now();
    auto glapsed_ans = std::chrono::duration_cast<std::chrono::milliseconds>(stop_ans - start_ans);
    std::cout << " answer time costs " << glapsed_ans.count() << " ms" << std::endl;

    // recover
    std::vector<uint64_t> decrypted(N);
    recover(decrypted, result, answerKey);

    std::cout << "the recovered value is " << decrypted[0] << std::endl;
}


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_larger_database();

    return 0;
}
