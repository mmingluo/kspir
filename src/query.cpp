#include "query.h"

#include "params.h"
#include "samples.h"
#include "encrypt.h"

void query(RlweCiphertext& cipher, Secret& queryKey, uint64_t row)
{
    uint64_t length = queryKey.getLength();

    // uint64_t* message = new uint64_t[length];
    std::vector<uint64_t> message(length);

    // TODO: check it
    if (row == 0)
        message[0] = Delta;
    else
        message[length-row] = bigMod - Delta;

    // store a in coefficient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();

    ntts.ComputeForward(message.data(), message.data(), 1, 1);
    encrypt(cipher.b.data(), cipher.a.data(), queryKey, message.data());

    ntts.ComputeInverse(cipher.a.data(), cipher.a.data(), 1, 1);
    cipher.setIsNtt(false); // note! b is ntt form, a is coefficient form.
#endif

    // delete[] message;
}

// the keyswitch key has size of [2N][N]
void build_keyswitch_key(uint64_t** ks, Secret& secret, std::vector<Secret>& encrypted_keys)
{
    uint64_t modulus = secret.getModulus();

    for (size_t i = 0; i < 2 * ell * N; i += 2)
    {
        sample_random(ks[i], modulus, N);
    }
    
    uint64_t* message = new uint64_t[N];

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < ell; j++)
        {
            for (size_t k = 0; k < N; k++)
            {
                message[k] = 256 * j * encrypted_keys[i].getData(k);
            }
            int index = 2 * ell * i + 2 * j;
            encrypt(ks[index + 1], ks[index], secret, message);
        }
    }

    delete[] message;
}
