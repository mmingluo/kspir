#include <stdint.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>

#include "pir.h"


void test_setup()
{   
    Secret queryKey(bigMod), answerKey(bigMod);

    // sample database
    std::vector<std::vector<uint64_t> > data(N, std::vector<uint64_t>(N, 0));
    sample_database(data);

    std::vector<std::vector<uint64_t> > setup_data(2 * N, std::vector<uint64_t>(N, 0));
    data_to_setupdata(setup_data, data, Pbits, bigMod, bigMod2);

    /************* dummy variant **************************************/
    database_to_signed(data, Pbits, bigMod);
    database_tontt(data);
    std::vector<std::vector<uint64_t> > ks2(2 * ell * N, std::vector<uint64_t>(N, 0));
    dummy_ksKeyGen(ks2, queryKey, answerKey, data);
    /******************************************************************/

    std::vector<std::vector<uint64_t> > pk(2 * ell * N * 2, std::vector<uint64_t>(N, 0));
    publicKeyGen(pk, queryKey, answerKey);

    // build keyswitch keys
    std::vector<std::vector<uint64_t> > ks(2 * ell * N * 2, std::vector<uint64_t>(N, 0));

    auto start = std::chrono::high_resolution_clock::now();

    evalKsKeyGen(ks, pk, setup_data);

    std::vector<std::vector<uint64_t> > result(2 * ell * N, std::vector<uint64_t>(N, 0));
    modSwitch(result, ks);
    // modSwitch(result, pk);

    auto stop = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << " evaluate kskeygen time costs " << glapsed.count() << " ms" << std::endl;

    std::vector<uint64_t> decrypted(N), decrypted_crt2(N), decrypted_crt1(N);
    // uint64_t crt_interval = 2 * ell * N;
    decrypt(decrypted.data(), result[0].data(), result[1].data(), answerKey);
    showLargeVector(decrypted, "the decrypted value is ");

    intel::hexl::NTT ntts = answerKey.getNTT();
    ntts.ComputeForward(ks2[0].data(), ks2[0].data(), 1, 1);
    ntts.ComputeForward(ks2[1].data(), ks2[1].data(), 1, 1);
    decrypt(decrypted.data(), ks2[0].data(), ks2[1].data(), answerKey);
    showLargeVector(decrypted, "the dummy is ");
    /*
    decrypt(decrypted_crt1.data(), pk[0].data(), pk[1].data(), answerKey);
    showLargeVector(decrypted_crt1, "crt1 is ");

    decrypt(decrypted_crt2.data(), pk[0 + crt_interval].data(), pk[1 + crt_interval].data(), answerKey, bigMod2);
    showLargeVector(decrypted_crt2, "crt2 is ");
    */
}


int main(int argc, char** argv)
{
    srand(time(NULL));
    
    test_setup();

    return 0;
}
