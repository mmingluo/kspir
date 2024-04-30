#include <stdint.h>
#include <cassert>
#include <string.h>

#include "encrypt.h"
#include "params.h"
#include "samples.h"
#include "ntt.h"
#include "utils.h"

void native_mult(uint64_t* b, uint64_t* a, Secret& secret);

// encrypt an LWE cipertext
void encrypt(uint64_t* b, uint64_t* a, Secret& secret, uint64_t message)
{
    uint64_t modulus = secret.getModulus();
    LWE_TYPE lweType = secret.getLweType();
    // assert(lweType == LWE, "this secret should be LWE type");
    assert(lweType == LWE);

    int32_t size = secret.getLength();
    // assert(N == size, "the size should be equal to that of ntt");
    assert(N == size);

    sample_random(a, modulus, size);

    int64_t result = 0;
    for (size_t i = 0; i < N; i++)
    {
        uint64_t value = secret.getData(i);
        // It would be faster to set data in public field.
        // uint64_t value = secret.data[i];
        if (value == 1)
        {
            result += a[i];
        } else if (value == modulus - 1)
        {
            result -= a[i];
        } else {}
        // *b += a[i] * secret.getData(i);
    }

    uint64_t err = sample_guass(modulus);
    result += message;
    result += err;

    *b = (result % (int64_t)modulus + modulus) % modulus;
    // *b = barret_reduce((int64_t)(*b));
}

// encrypt an RLWE ciphertext
// a, b are in ntt form
void encrypt(uint64_t* b, uint64_t* a, Secret& secret, uint64_t* message)
{
    uint64_t modulus = secret.getModulus();
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    sample_random(a, modulus, length);

    // native_mult(b, a, secret);
#ifndef INTEL_HEXL
    static_assert(N == dim, "the dimension should be equal to that of ntt");

    // TODO: adjust uint64_t ntt
    int64_t* tempa = new int64_t[N];
    int64_t* tempb = new int64_t[N];
    copy_vector(tempa, a); // a[i] < 2^63

    ntt(tempa);

    hadamard_mult(tempb, tempa, secret.getData());
    invntt_tomont(tempb);

    copy_vector(b, tempb, modulus);
    
    uint64_t* err = new uint64_t[N];
    sample_guass(err, modulus);
    
    for (size_t i = 0; i < N; i++)
    {
        b[i] += message[i];
        b[i] += err[i];
    }

    delete[] tempa;
    delete[] tempb;
    delete[] err;
#endif

#ifdef INTEL_HEXL
    if (!secret.isNttForm())
    {
        secret.toNttForm();
    }

    intel::hexl::EltwiseMultMod(b, a, secret.getData().data(), length,
                    modulus, 1);

    intel::hexl::NTT ntts = secret.getNTT();

    std::vector<uint64_t> err(length);
    // TODO:
    sample_guass(err.data(), modulus);

    ntts.ComputeForward(err.data(), err.data(), 1, 1);

    intel::hexl::EltwiseAddMod(b, b, err.data(), length, modulus);
    intel::hexl::EltwiseAddMod(b, b, message, length, modulus);
#endif
}

// TODO: another modulus
// encrypt an CRT RLWE ciphertext (CRT num = 2)
// a, b are in ntt form
#ifdef INTEL_HEXL
void encrypt(uint64_t* b1, uint64_t* a1, uint64_t* b2, uint64_t* a2,
            Secret& secret, uint64_t* message,
            uint64_t* secret1, uint64_t* secret2,
            intel::hexl::NTT& ntts1, intel::hexl::NTT& ntts2,
            uint64_t modulus1, uint64_t modulus2)
{
    uint64_t modulus = secret.getModulus();
    assert(modulus1 == modulus);
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    sample_random(a1, modulus1, length);
    sample_random(a2, modulus2, length);

    std::vector<uint64_t> temp(length);

    intel::hexl::EltwiseMultMod(b1, a1, secret1, length,
                    modulus1, 1);
    intel::hexl::EltwiseMultMod(b2, a2, secret2, length,
                    modulus2, 1);

    std::vector<uint64_t> err(length);
    sample_guass(err.data(), modulus1);

    // the two RNS part share the sample meessage and error
    // intel::hexl::EltwiseFMAMod(temp.data(), message, modulus2, nullptr, length, modulus1, 1);
    // intel::hexl::EltwiseAddMod(temp.data(), temp.data(), err.data(), length, modulus1);
    intel::hexl::EltwiseFMAMod(temp.data(), message, modulus2, err.data(), length, modulus1, 1);
    ntts1.ComputeForward(temp.data(), temp.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b1, b1, temp.data(), length, modulus1);

    guass_to_modulus(err.data(), modulus1, modulus2);
    // encode: q2*m % q2 = 0, so encrypt 0 to second part
    ntts2.ComputeForward(temp.data(), err.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b2, b2, temp.data(), length, modulus2);
}
#endif

// rns encrypt for two messages
#ifdef INTEL_HEXL
void encrypt(uint64_t* b1, uint64_t* a1, uint64_t* b2, uint64_t* a2,
            Secret& secret, uint64_t* message1, uint64_t* message2,
            uint64_t* secret1, uint64_t* secret2,
            intel::hexl::NTT& ntts1, intel::hexl::NTT& ntts2,
            uint64_t modulus1, uint64_t modulus2, int32_t num = 0)
{
    uint64_t modulus = secret.getModulus();
    assert(modulus1 == modulus);
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    sample_random(a1, modulus1, length);
    sample_random(a2, modulus2, length);

    std::vector<uint64_t> temp(length);

    intel::hexl::EltwiseMultMod(b1, a1, secret1, length,
                    modulus1, 1);
    intel::hexl::EltwiseMultMod(b2, a2, secret2, length,
                    modulus2, 1);

    std::vector<uint64_t> err(length), err_t(length);
    sample_guass(err.data(), modulus1);

    /************* multipy error with N = 2^l in this variant  *********/
    intel::hexl::EltwiseFMAMod(err_t.data(), err.data(), (uint64_t)num,
                        nullptr, length, modulus, 1);
    /*******************************************************************/


    // the two RNS part share the sample meessage and error
    intel::hexl::EltwiseAddMod(temp.data(), err_t.data(), message1, length, modulus1);
    ntts1.ComputeForward(temp.data(), temp.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b1, b1, temp.data(), length, modulus1);

    guass_to_modulus(err.data(), modulus1, modulus2);
    /************* multipy error with N = 2^l in this variant  *********/
    intel::hexl::EltwiseFMAMod(err_t.data(), err.data(), (uint64_t)num,
                        nullptr, length, modulus, 1);
    /*******************************************************************/

    // TODO: we modify it
    // intel::hexl::EltwiseAddMod(temp.data(), err_t.data(), message2, length, modulus1);
    intel::hexl::EltwiseAddMod(temp.data(), err_t.data(), message2, length, modulus2);
    ntts2.ComputeForward(temp.data(), temp.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b2, b2, temp.data(), length, modulus2);
}
#endif

void encrypt_rns(RlweCiphertext& cipher1, RlweCiphertext& cipher2,
                Secret& secret, std::vector<uint64_t>& message)
{
    uint64_t length = secret.getLength();
    uint64_t modulus1 = secret.getModulus();
    uint64_t modulus2 = bigMod2;

    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = secret.getNTT();
    intel::hexl::NTT ntts2(length, modulus2);

    std::vector<uint64_t> secret1(length);
    std::vector<uint64_t> secret2(length);

    for (size_t i = 0; i < secret.getLength(); i++)
    {
        secret1[i] = secret.getData(i);
        secret2[i] = secret.getData(i);
    }
    guass_to_modulus(secret2.data(), modulus1, modulus2);
    ntts1.ComputeForward(secret1.data(), secret1.data(), 1, 1);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    // encode message * bigMod2
    encode_crt(message);

    // CRT encryption
    encrypt(cipher1.b.data(), cipher1.a.data(),
            cipher2.b.data(), cipher2.a.data(),
            secret, message.data(),
            secret1.data(), secret2.data(), // should be ntt form
            ntts1, ntts2, modulus1, modulus2);
    cipher1.setIsNtt(true);
    cipher2.setIsNtt(true);

#endif
}

/**
 * @brief encrypt rns ciphertext for message
 * 
 * @param b1, a1 rns part 1 
 * @param b2, a2 rns part 2
 * @param secret 
 * @param message1 coefficient form
 * @param message2 coefficient form
 */
void encrypt_rns(std::vector<uint64_t>& b1, std::vector<uint64_t>& a1,
                    std::vector<uint64_t>& b2, std::vector<uint64_t>& a2,
                    Secret& secret, std::vector<uint64_t>& message1, std::vector<uint64_t>& message2,
                    int32_t num, uint64_t modulus2)
{
    uint64_t length = secret.getLength();
    uint64_t modulus1 = secret.getModulus();
    // uint64_t modulus2 = bigMod2;

    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = secret.getNTT();
    intel::hexl::NTT ntts2(length, modulus2);

    std::vector<uint64_t> secret1(length);
    std::vector<uint64_t> secret2(length);

    for (size_t i = 0; i < secret.getLength(); i++)
    {
        secret1[i] = secret.getData(i);
        secret2[i] = secret.getData(i);
    }
    guass_to_modulus(secret2.data(), modulus1, modulus2);
    ntts1.ComputeForward(secret1.data(), secret1.data(), 1, 1);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    // CRT encryption
    encrypt(b1.data(), a1.data(), b2.data(), a2.data(),
            secret, message1.data(), message2.data(),
            secret1.data(), secret2.data(), // should be ntt form
            ntts1, ntts2, modulus1, modulus2, num);

#endif
}


// rns encrypt for two messages
#ifdef INTEL_HEXL
void encrypt_bsgs_autokey(uint64_t* b1, uint64_t* a1, uint64_t* b2, uint64_t* a2,
            Secret& secret, uint64_t* message1, uint64_t* message2,
            uint64_t* secret1, uint64_t* secret2,
            intel::hexl::NTT& ntts1, intel::hexl::NTT& ntts2,
            uint64_t modulus1, uint64_t modulus2)
{
    uint64_t modulus = secret.getModulus();
    assert(modulus1 == modulus);
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    sample_random(a1, modulus1, length);
    sample_random(a2, modulus2, length);

    std::vector<uint64_t> temp(length);

    intel::hexl::EltwiseMultMod(b1, a1, secret1, length,
                    modulus1, 1);
    intel::hexl::EltwiseMultMod(b2, a2, secret2, length,
                    modulus2, 1);

    std::vector<uint64_t> err(length);
    sample_guass(err.data(), modulus1);

    // the two RNS part share the sample meessage and error
    intel::hexl::EltwiseAddMod(temp.data(), err.data(), message1, length, modulus1);
    ntts1.ComputeForward(temp.data(), temp.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b1, b1, temp.data(), length, modulus1);

    guass_to_modulus(err.data(), modulus1, modulus2);

    // the second part
    intel::hexl::EltwiseAddMod(temp.data(), err.data(), message2, length, modulus2);
    ntts2.ComputeForward(temp.data(), temp.data(), 1, 1);
    intel::hexl::EltwiseAddMod(b2, b2, temp.data(), length, modulus2);
}
#endif

/**
 * @brief encrypt rns ciphertext for message
 * 
 * @param b1, a1 rns part 1 
 * @param b2, a2 rns part 2
 * @param secret 
 * @param message1 coefficient form
 * @param message2 coefficient form
 */
void encrypt_rns_bsgs_autokey(std::vector<uint64_t>& b1, std::vector<uint64_t>& a1,
                    std::vector<uint64_t>& b2, std::vector<uint64_t>& a2,
                    Secret& secret, std::vector<uint64_t>& message1, std::vector<uint64_t>& message2,
                    uint64_t modulus2)
{
    uint64_t length = secret.getLength();
    uint64_t modulus1 = secret.getModulus();
    // uint64_t modulus2 = bigMod2;

    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = secret.getNTT();
    intel::hexl::NTT ntts2(length, modulus2);

    std::vector<uint64_t> secret1(length);
    std::vector<uint64_t> secret2(length);

    /**
    for (size_t i = 0; i < secret.getLength(); i++)
    {
        secret1[i] = secret.getData(i);
        secret2[i] = secret.getData(i);
    }
    **/
    memcpy(secret1.data(), secret.getData().data(), sizeof(uint64_t) * length);
    memcpy(secret2.data(), secret.getData().data(), sizeof(uint64_t) * length);
    guass_to_modulus(secret2.data(), modulus1, modulus2);
    ntts1.ComputeForward(secret1.data(), secret1.data(), 1, 1);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    // CRT encryption
    encrypt_bsgs_autokey(b1.data(), a1.data(), b2.data(), a2.data(),
            secret, message1.data(), message2.data(),
            secret1.data(), secret2.data(), // should be ntt form
            ntts1, ntts2, modulus1, modulus2);
#endif
}

// encrypt an RLWE ciphertext only for packing. refer to 
// a, b are in form
void encrypt_special_rlwe(uint64_t* b, uint64_t* a, Secret& secret, uint64_t* message, int32_t lwesnum)
{
    uint64_t modulus = secret.getModulus();
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    sample_random(a, modulus, length);

#ifdef INTEL_HEXL
    if (!secret.isNttForm())
    {
        secret.toNttForm();
    }

    intel::hexl::EltwiseMultMod(b, a, secret.getData().data(), length,
                    modulus, 1);

    intel::hexl::NTT ntts = secret.getNTT();

    std::vector<uint64_t> err(length);
    sample_guass(err.data(), modulus);

    /************** multipy error with 2^l in this variant  ************/
    intel::hexl::EltwiseFMAMod(err.data(), err.data(), (uint64_t)lwesnum,
                        nullptr, length, modulus, 1);
    /*******************************************************************/

    ntts.ComputeForward(err.data(), err.data(), 1, 1);

    intel::hexl::EltwiseAddMod(b, b, err.data(), length, modulus);
    intel::hexl::EltwiseAddMod(b, b, message, length, modulus);
#endif
}

void encrypt(RlweCiphertext& cipher, Secret& secret, std::vector<uint64_t>& message)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = secret.getNTT();
    ntts.ComputeForward(message.data(), message.data(), 1, 1);
#endif
    // the message should be ntt form in the following `encrypt' function
    encrypt(cipher.b.data(), cipher.a.data(), secret, message.data());
    cipher.setIsNtt(true);
}

/**
 * @brief encrypt for rns bsgs algorithm
 * @param cipher encrypted cipher
 * @param secret secret
 * @param message stored in coefficients form
 * 
*/
void encrypt_rns_bsgs(std::vector<RlweCiphertext>& cipher, Secret& secret, std::vector<uint64_t>& message)
{
    uint64_t modulus = secret.getModulus();
    uint64_t bsModulus = bsMod;
    uint64_t length = secret.getLength();

    // a is sampling in ntt form.
    // sample_random(cipher[0].a.data(), modulus, length);
    // sample_random(cipher[1].a.data(), bsModulus, length);

#ifdef INTEL_HEXL
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    intel::hexl::NTT ntts1 = secret.getNTT();
    intel::hexl::NTT ntts2(length, bsModulus);

    std::vector<uint64_t> secret1(length);
    std::vector<uint64_t> secret2(length);

    /**
    for (size_t i = 0; i < secret.getLength(); i++)
    {
        secret1[i] = secret.getData(i);
        secret2[i] = secret.getData(i);
    }
    **/
    memcpy(secret1.data(), secret.getData().data(), sizeof(uint64_t) * secret.getLength());
    memcpy(secret2.data(), secret.getData().data(), sizeof(uint64_t) * secret.getLength());
    guass_to_modulus(secret2.data(), modulus, bsModulus);
    ntts1.ComputeForward(secret1.data(), secret1.data(), 1, 1);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    // encode message * bigMod2
    // TODO: check
    // encode_crt(message);

    // CRT encryption
    encrypt(cipher[0].b.data(), cipher[0].a.data(),
            cipher[1].b.data(), cipher[1].a.data(),
            secret, message.data(),
            secret1.data(), secret2.data(), // should be ntt form
            ntts1, ntts2, modulus, bsModulus);
    ntts1.ComputeInverse(cipher[0].b.data(), cipher[0].b.data(), 1, 1);
    ntts1.ComputeInverse(cipher[0].a.data(), cipher[0].a.data(), 1, 1);
    ntts2.ComputeInverse(cipher[1].b.data(), cipher[1].b.data(), 1, 1);
    ntts2.ComputeInverse(cipher[1].a.data(), cipher[1].a.data(), 1, 1);
    cipher[0].setIsNtt(false);
    cipher[1].setIsNtt(false);
#endif
}

// decrypt an LWE cipertext
uint64_t decrypt(uint64_t* b, uint64_t* a, Secret& secret)
{
    uint64_t modulus = secret.getModulus();
    LWE_TYPE lweType = secret.getLweType();
    // assert(lweType == LWE, "this secret should be LWE type");
    assert(lweType == LWE);

    int32_t size = secret.getLength();
    assert(N == size);

    // TODO: if Q > 2^56, please change the following code
    int64_t as = 0;
    for (size_t i = 0; i < N; i++)
    {
        int64_t value = secret.getData(i);
        if (value == modulus - 1)
        {
            value = -1;
        }
        else {}
        as += a[i] * value;
    }
    int64_t reduce = ((int64_t)*b - as) % (int64_t)modulus;
    uint64_t message = (reduce + modulus) % modulus;

    return message;
}

// decrypt an RLWE ciphertext
// a, b are in ntt form
void decrypt(uint64_t* message, uint64_t* b, uint64_t* a, Secret& secret)
{
    uint64_t modulus = secret.getModulus();
    uint64_t length = secret.getLength();

#ifndef INTEL_HEXL
    uint64_t* as = new uint64_t[N];
    // native_mult(as, a, secret);
        // TODO: adjust uint64_t ntt
    int64_t* tempa = new int64_t[N];
    int64_t* tempb = new int64_t[N];
    copy_vector(tempa, a); // a[i] < 2^63

    ntt(tempa);

    hadamard_mult(tempb, tempa, secret.getData());
    invntt_tomont(tempb);

    copy_vector(as, tempb, modulus);

    for (size_t i = 0; i < N; i++)
    {
        message[i] = b[i] - as[i];
    }

    delete[] as;
    delete[] tempa;
    delete[] tempb;
#endif

#ifdef INTEL_HEXL
    if (!secret.isNttForm())
    {
        secret.toNttForm();
    }

    intel::hexl::EltwiseMultMod(message, a, secret.getData().data(), length,
                    modulus, 1);

    intel::hexl::EltwiseSubMod(message, b, message, length, modulus);

    intel::hexl::NTT ntts = secret.getNTT();
    ntts.ComputeInverse(message, message, 1, 1);
#endif

}

void decrypt(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& secret)
{
    if (!cipher.getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = secret.getNTT();
        ntts.ComputeForward(cipher.a.data(), cipher.a.data(), 1, 1);
        ntts.ComputeForward(cipher.b.data(), cipher.b.data(), 1, 1);
        cipher.setIsNtt(true);
#endif
    }
    
    decrypt(message.data(), cipher.b.data(), cipher.a.data(), secret);
}

void decrypt(uint64_t* message, uint64_t* b, uint64_t* a, Secret& secret, uint64_t modulus2)
{
    uint64_t modulus = secret.getModulus();
    uint64_t length = secret.getLength();

#ifdef INTEL_HEXL
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    intel::hexl::NTT ntts2(length, modulus2);
    std::vector<uint64_t> secret2(length);
    for (size_t i = 0; i < secret.getLength(); i++)
    {
        secret2[i] = secret.getData(i);
    }
    guass_to_modulus(secret2.data(), modulus, modulus2);
    ntts2.ComputeForward(secret2.data(), secret2.data(), 1, 1);

    intel::hexl::EltwiseMultMod(message, a, secret2.data(), length,
                    modulus2, 1);

    intel::hexl::EltwiseSubMod(message, b, message, length, modulus2);

    ntts2.ComputeInverse(message, message, 1, 1);
#endif

}


void native_mult(uint64_t* b, uint64_t* a, Secret& secret)
{
    uint64_t modulus = secret.getModulus();
    //TODO: add modulus
    for (size_t i = 0; i < N; i++)
    {
        b[i] = 0;
        for (size_t j = 0; j <= i; j++)
        {
            b[i] += (a[j] * secret.getData(i-j));
            b[i] &= (modulus-1);
        }
        for (size_t j = i+1; j < N; j++)
        {
            b[i] -= a[j] * secret.getData(N-j+i);
            b[i] &= (modulus-1);
        }
    }
}
