#include <iostream>
#include <chrono>
#include <time.h>

// #include "answer.h"
#include "pir.h"

/**
 * @brief this function just for test in dummy database mode.
 * 
 */
void test_dummy_data()
{
    // database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N, 0));

    auto start = std::chrono::high_resolution_clock::now();
    sample_database(data);
    database_tontt(data);
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "database initialization time costs " << glapsed.count() << " us." << std::endl;

    // query
    RlweCiphertext query(N, bigMod);
    sample_random(query.a, bigMod);
    sample_random(query.b, bigMod);

    // the results (null value)
    RlweCiphertext result(N, bigMod);
    sample_random(result.b, bigMod);

    // keySwitch key (null value)
    std::vector<std::vector<uint64_t> > ks(2 * N * ell, std::vector<uint64_t>(N, 0));
    /*std::vector<RlweCiphertext> ks(N * ell);
    for (size_t i = 0; i < N * ell; i++)
    {
        sample_random62_vector(ks[i].getA().data(), ks[i].getLength());
        sample_random62_vector(ks[i].getA().data(), ks[i].getLength());
    }
    */
    /*
    for (size_t i = 0; i < 2 * N * ell; i++)
    {
        // sample_random62_vector(ks[i].data(), ks[i].size());
        sample_random(ks[i], bigMod);
    }
    */
    start = std::chrono::high_resolution_clock::now();    
    
int32_t ntimes = 1;
for (size_t i = 0; i < ntimes; i++)
{
    // sample_random(query.a, bigMod);
    // sample_random(query.b, bigMod);

#ifndef INTEL_HEXL
    // answer(result_b.data(), result_a.data(), query_b.data(), query_a.data(), data, ks);
#else
    answer_first_dimension(result.b, result.a, query.b, data);
#endif
    std::cout << "(  " << result.b[0] << "   )   " << std::endl;
}
    stop = std::chrono::high_resolution_clock::now();

    glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " time costs " << glapsed.count() << " us." << std::endl;

}



int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_dummy_data();

    return 0;
}
