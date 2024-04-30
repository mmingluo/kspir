#include <iostream>
#include <math.h>
#include <assert.h>

#include "lwe.h"
#include "samples.h"
#include "utils.h"
#include "encrypt.h"

#include "crt.h"

#ifdef INTEL_HEXL
    #include <hexl/hexl.hpp>
#endif

/**
 * LweCiphertext class
*/
    
LweCiphertext::LweCiphertext(/* args */){}

LweCiphertext::LweCiphertext(int32_t len, uint64_t module) : length(len), modulus(module) 
{
    this->a.resize(len);
}

std::vector<uint64_t> LweCiphertext::getA()
{
    return a;
}

uint64_t LweCiphertext::getA(int32_t index)
{
    return a[index];
}

uint64_t LweCiphertext::getB()
{
    return b;
}

uint64_t LweCiphertext::getModulus()
{
    return modulus;
}

int32_t LweCiphertext::getLength()
{
    return length;
}


/**
 * RlweCiphertext class
*/

RlweCiphertext::RlweCiphertext(/* args */){
    this->a.resize(length);
    this->b.resize(length);
}

RlweCiphertext::RlweCiphertext(uint64_t len, uint64_t module) : length(len), modulus(module) 
{
    this->a.resize(len);
    this->b.resize(len);
    // sample_random(this->a, module);
    // TODO: should be removed
    // sample_random(this->b, module);
}

RlweCiphertext::RlweCiphertext(std::vector<uint64_t>& a, std::vector<uint64_t>& b, bool isn): isntt(isn)
{
    assert(N == a.size());
    assert(N == b.size());

    copy(a.begin(), a.end(), this->a.begin());
    copy(b.begin(), b.end(), this->b.begin());
}

void RlweCiphertext::addAndEqual(const RlweCiphertext& cipher)
{

#ifdef INTEL_HEXL
    intel::hexl::EltwiseAddMod(a.data(), a.data(), cipher.a.data(),
                this->getLength(), this->getModulus());
    intel::hexl::EltwiseAddMod(b.data(), b.data(), cipher.b.data(),
                this->getLength(), this->getModulus());
#endif
}

void RlweCiphertext::subAndEqual(const RlweCiphertext& cipher)
{

#ifdef INTEL_HEXL
    intel::hexl::EltwiseSubMod(a.data(), a.data(), cipher.a.data(),
                this->getLength(), this->getModulus());
    intel::hexl::EltwiseSubMod(b.data(), b.data(), cipher.b.data(),
                this->getLength(), this->getModulus());
#endif
}

void RlweCiphertext::multConst(int32_t number)
{
    std::vector<uint64_t> temp_a(length), temp_b(length);
    std::copy(a.begin(), a.end(), temp_a.begin());
    std::copy(b.begin(), b.end(), temp_b.begin());

    for (size_t i = 0; i < length - number; i++)
    {
        a[i + number] = temp_a[i];
        b[i + number] = temp_b[i];
    }
    for (size_t i = length - number; i < length; i++)
    {
        int32_t index = i + number - length;
        a[index] = modulus - temp_a[i];
        b[index] = modulus - temp_b[i];
    }
}

/**
 * @brief mult negate of input const. This function is used in the expand function.
 *   input ciphertext should be in coefficient form.
 * @param neg_factor 
 */
void RlweCiphertext::multNegConst(int32_t neg_factor)
{
    uint64_t number = length - neg_factor;

    std::vector<uint64_t> temp_a(length), temp_b(length);
    std::copy(a.begin(), a.end(), temp_a.begin());
    std::copy(b.begin(), b.end(), temp_b.begin());

    for (size_t i = 0; i < length - number; i++)
    {
        a[i + number] = modulus - temp_a[i];
        b[i + number] = modulus - temp_b[i];
    }
    for (size_t i = length - number; i < length; i++)
    {
        int32_t index = i + number - length;
        a[index] = temp_a[i];
        b[index] = temp_b[i];
    }
}


std::vector<uint64_t> RlweCiphertext::getA() const
{
    return this->a;
}

uint64_t RlweCiphertext::getA(int32_t index)
{
    return this->a[index];
}

std::vector<uint64_t> RlweCiphertext::getB() const
{
    return this->b;
}

uint64_t RlweCiphertext::getB(int32_t index)
{
    return this->b[index];
}

uint64_t RlweCiphertext::getModulus() const
{
    return this->modulus;
}

void RlweCiphertext::setModulus(uint64_t module)
{
    this->modulus = module;
}

int32_t RlweCiphertext::getLength() const
{
    return this->length;
}

bool RlweCiphertext::getIsNtt() const
{
    return this->isntt;
}

void RlweCiphertext::setIsNtt(bool is)
{
    this->isntt = is;
}

RlweCiphertext& RlweCiphertext::operator=(RlweCiphertext& cipher)
{
    copy(cipher.a.cbegin(), cipher.a.cend(), this->a.begin());
    copy(cipher.b.cbegin(), cipher.b.cend(), this->b.begin());

    this->isntt = cipher.isntt;
    this->length = cipher.length;
    this->modulus = cipher.modulus;

    return *this;
}
/*
void RlweCiphertext::toCoeffForm()
{
     // TODO
#ifdef INTEL_HEXL
    ntts.ComputeInverse(this->a.data(), this->a.data(), 1, 1);
    ntts.ComputeInverse(this->b.data(), this->b.data(), 1, 1);
#endif
    this->nttform = false;  
}

void RlweCiphertext::toNttForm()
{
    // TODO
#ifdef INTEL_HEXL
    ntts.ComputeForward(this->a.data(), this->a.data(), 1, 1);
    ntts.ComputeForward(this->b.data(), this->b.data(), 1, 1);
#endif
    this->nttform = true;
}
*/


/**
 * RGSWCiphertext class
*/

RGSWCiphertext::RGSWCiphertext(/* args */)
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        uint64_t length = this->getLength();
        this->a[i].resize(length);
        this->b[i].resize(length);
    }
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    this->ntts = ntts;
#endif

}

RGSWCiphertext::RGSWCiphertext(uint64_t len, uint64_t module) : length(len), modulus(module) 
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        this->a[i].resize(len);
        this->b[i].resize(len);
    }
    // sample_random(this->a, module);
    // TODO: should be removed
    // sample_random(this->b, module);


#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }
#endif
}

RGSWCiphertext::RGSWCiphertext(uint64_t len, uint64_t module, uint64_t elln, uint64_t base, uint64_t bg):
    length(len), modulus(module), ellnum(elln), PP(base), BBg(bg)
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        this->a[i].resize(len);
        this->b[i].resize(len);
    }

#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }
#endif
} 

void RGSWCiphertext::keyGen(Secret& secret, const uint64_t index, bool is_reverse)
{
    // the index should belong to [0, N)
    assert(0 <= index);
    assert(index < length);
    
    uint64_t length = secret.getLength();
    uint64_t modulus = secret.getModulus();

    std::vector<uint64_t> message(length);
    if (!is_reverse)
    {
        message[index] = 1;
    } else {
        if (index == 0)
            message[0] = 1;
        else
            message[length - index] = modulus - 1;
    }

    // uint64_t mult_factor;
    std::vector<uint64_t> temp(length);

    // m(s, 1)
    if(!secret.isNttForm())
    {
        secret.toNttForm();
    }
    intel::hexl::NTT ntts = secret.getNTT();

    std::vector<uint64_t> mult_factor = powerOfBg(this->PP, this->BBg, this->ellnum);
    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        // faster by single barret reduce
        intel::hexl::EltwiseFMAMod(temp.data(), message.data(),
                mult_factor[i], nullptr, this->length, this->modulus, 1);

        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);

        // rlwe encryption, these rlwes is in ntt form.
        // encrypt(m)
        encrypt(this->b[ellnum + i].data(), this->a[ellnum + i].data(), secret, temp.data());

        // encrypt(-ms)
        // intel::hexl::EltwiseSubMod(temp.data(), this->b[ellnum - 1].data(), temp.data(),
        //        this->length, this->modulus);
        negate(temp.data(), this->length, this->modulus);

        intel::hexl::EltwiseMultMod(temp.data(), temp.data(), secret.getData().data(),
                this->length, this->modulus, 1);
        encrypt(this->b[i].data(), this->a[i].data(), secret, temp.data());
#endif
        this->isntt = true;
    }
}

/*
std::vector<uint64_t> RGSWCiphertext::getA() const
{
    return this->a;
}

uint64_t RGSWCiphertext::getA(int32_t index)
{
    return this->a[index];
}

std::vector<uint64_t> RGSWCiphertext::getB() const
{
    return this->b;
}

uint64_t RGSWCiphertext::getB(int32_t index)
{
    return this->b[index];
}
*/
uint64_t RGSWCiphertext::getModulus() const
{
    return this->modulus;
}

uint64_t RGSWCiphertext::getLength() const
{
    return this->length;
}

uint64_t RGSWCiphertext::getEllnum() const
{
    return this->ellnum;
}

uint64_t RGSWCiphertext::getBg() const
{
    return this->BBg;
}

uint64_t RGSWCiphertext::getBase() const
{
    return this->PP;
}

bool RGSWCiphertext::getIsNtt() const
{
    return this->isntt;
}

#ifdef INTEL_HEXL
intel::hexl::NTT RGSWCiphertext::getNTT() const
{
    return ntts;
}
#endif


/**
 * AutoKey class
*/

AutoKey::AutoKey(/* args */)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    this->ntts = ntts;
#endif  
}

AutoKey::AutoKey(int32_t len, uint64_t module, int32_t elln, uint64_t base, uint64_t bg) :
         length(len), modulus(module), ellnum(elln), PP(base), BBg(bg)
{

#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }
#endif
}

void AutoKey::generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                                Secret& secret)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp(this->length), autosec(this->length);
    std::vector<uint64_t> mult_factor = powerOfBg(this->PP, this->BBg, this->ellnum);

    automorphic(autosec, secret.getData(), index, modulus);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        intel::hexl::EltwiseFMAMod(temp.data(), autosec.data(),
                mult_factor[i], nullptr, this->length, this->modulus, 1);
        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
#endif
        // rlwe encryption, these rlwes is in ntt form.
        encrypt(result[i].b.data(), result[i].a.data(), secret, temp.data());
        // encrypt_special_rlwe(result[i].b.data(), result[i].a.data(), secret, temp.data(), num);
        result[i].setIsNtt(true);
        this->isntt = true;
    }
}

void AutoKey::generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                        int32_t num, Secret& secret)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp(this->length), autosec(this->length);
    std::vector<uint64_t> mult_factor = powerOfBg(this->PP, this->BBg, this->ellnum);

    automorphic(autosec, secret.getData(), index, modulus);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        intel::hexl::EltwiseFMAMod(temp.data(), autosec.data(),
                mult_factor[i], nullptr, this->length, this->modulus, 1);
        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
#endif
        // rlwe encryption, these rlwes is in ntt form.
        // encrypt(result[i].b.data(), result[i].a.data(), secret, temp.data());
        encrypt_special_rlwe(result[i].b.data(), result[i].a.data(), secret, temp.data(), num);
        result[i].setIsNtt(true);
        this->isntt = true;
    }
}

/**
 * @brief genreate automorphic keys
 * 
 * @param secret 
 * @param num packing numbers
 */
void AutoKey::keyGen(Secret& secret, const int32_t num)
{
    std::vector<int32_t> indexLists;
    // Algorithm 2: PackLWEs
    for (size_t i = 2; i <= num + 1; i <<= 1)
    {
        indexLists.push_back(i + 1);
    }
    int32_t log2Nn = log2(length / num);
    // Algorithm 1: EvalTrNn
    for (size_t k = 1; k <= log2Nn; k++)
    {
        int32_t index = (length >> k << 1) + 1;
        indexLists.push_back(index);
    }

    std::vector<RlweCiphertext> key_for_index(this->ellnum);

    for (auto index : indexLists)
    {
        // encrypt automorphic^()(s)
        this->generateSingleKey(key_for_index, index, num, secret);

        keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
    }
}

/**
 * @brief genreate automorphic keys
 * 
 * @param secret 
 * @param num packing numbers
 */
void AutoKey::keyGen(Secret& secret, const std::vector<int32_t> indexLists)
{
    std::vector<RlweCiphertext> key_for_index(this->ellnum);

    for (auto index : indexLists)
    {
        // encrypt automorphic^()(s)
        this->generateSingleKey(key_for_index, index, secret);

        keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
    }
}

uint64_t AutoKey::getModulus() const
{
    return this->modulus;
}

uint64_t AutoKey::getLength() const
{
    return this->length;
}

uint64_t AutoKey::getEllnum() const
{
    return this->ellnum;
}

uint64_t AutoKey::getBg() const
{
    return this->BBg;
}

uint64_t AutoKey::getBase() const
{
    return this->PP;
}

bool AutoKey::getIsNtt() const
{
    return this->isntt;
}

#ifdef INTEL_HEXL
intel::hexl::NTT AutoKey::getNTT() const
{
    return ntts;
}
#endif


/**
 * AutoKeyRNS class
*/

AutoKeyRNS::AutoKeyRNS(/* args */)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1(length, modulus1);
    this->ntts1 = ntts1;

    intel::hexl::NTT ntts2(length, modulus2);
    this->ntts2 = ntts2;
#endif  
}

AutoKeyRNS::AutoKeyRNS(int32_t len, uint64_t module1, uint64_t module2, int32_t elln) :
         length(len), modulus1(module1), modulus2(module2), ellnum(elln)
{

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1(length, modulus1);
    this->ntts1 = ntts1;

    intel::hexl::NTT ntts2(length, modulus2);
    this->ntts2 = ntts2;
#endif
}

// TODO: fix it
void AutoKeyRNS::generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                        int32_t num, Secret& secret)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp1(this->length), temp2(this->length);
    std::vector<uint64_t> autosec(this->length), autosec2(this->length);
    // note the powerofBg for 128 bits
    std::vector<uint128_t> mult_factor = powerOfBg128(this->PP, this->BBg, this->ellnum);

    automorphic(autosec, secret.getData(), index, modulus1);

    // store in another modulus
    copy(autosec.begin(), autosec.end(), autosec2.begin());
    guass_to_modulus(autosec2.data(), modulus1, modulus2);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        uint64_t temp = mult_factor[i] % this->modulus1;
        intel::hexl::EltwiseFMAMod(temp1.data(), autosec.data(),
                temp, nullptr, this->length, this->modulus1, 1);

        temp = mult_factor[i] % this->modulus2;
        intel::hexl::EltwiseFMAMod(temp2.data(), autosec2.data(),
                temp, nullptr, this->length, this->modulus2, 1);

        // rns encryption, these rlwes is in ntt form.
        encrypt_rns(result[i].b, result[i].a, result[i + ellnum].b, result[i + ellnum].a,
                    secret, temp1, temp2);
        //encrypt_special_rns(result[i].b, result[i].a, result[i + ellnum].b, result[i + ellnum].a,
        //            secret, temp1, temp2, num);
        
        result[i].setIsNtt(true);
        result[i + ellnum].setIsNtt(true);
        result[i + ellnum].setModulus(modulus2);
#endif
    }
}

void AutoKeyRNS::keyGen(Secret& secret, const int32_t num)
{
    std::vector<int32_t> indexLists;

    int32_t log2N = log2(length);
    // Algorithm: Eval expand
    for (size_t k = 1; k <= log2N; k++)
    {
        int32_t index = (length >> k << 1) + 1;
        indexLists.push_back(index);
    }

    std::vector<RlweCiphertext> key_for_index(this->ellnum * 2);

    for (auto index : indexLists)
    {
        this->generateSingleKey(key_for_index, index, num, secret);

        keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
    }
}

uint64_t AutoKeyRNS::getModulus1() const
{
    return this->modulus1;
}

uint64_t AutoKeyRNS::getModulus2() const
{
    return this->modulus2;
}

uint64_t AutoKeyRNS::getLength() const
{
    return this->length;
}

uint64_t AutoKeyRNS::getEllnum() const
{
    return this->ellnum;
}

uint64_t AutoKeyRNS::getBg() const
{
    return this->BBg;
}

uint64_t AutoKeyRNS::getBase() const
{
    return this->PP;
}

bool AutoKeyRNS::getIsNtt() const
{
    return this->isntt;
}

#ifdef INTEL_HEXL
intel::hexl::NTT AutoKeyRNS::getNTT1() const
{
    return ntts1;
}

intel::hexl::NTT AutoKeyRNS::getNTT2() const
{
    return ntts2;
}

#endif

/**
 * @brief 
 * 
 * @param result coefficient forms
 * @param firstDim  form
 * @param query ntt form
 */
void externalProduct(RlweCiphertext& result, RlweCiphertext& firstDim, const RGSWCiphertext& query)
{
    assert(query.getIsNtt());

    uint64_t ellnum = query.getEllnum();
    uint64_t length = query.getLength();
    uint64_t modulus = query.getModulus();

#ifdef INTEL_HEXL
    // keyswitch
    std::vector<std::vector<uint64_t> > decFirstDim(2 * ellnum, std::vector<uint64_t>(length, 0));
    std::vector<uint64_t> temp(length, 0);
    intel::hexl::NTT ntts = query.getNTT();

    if (firstDim.getIsNtt())
    {
        ntts.ComputeInverse(firstDim.a.data(), firstDim.a.data(), 1, 1);
        ntts.ComputeInverse(firstDim.b.data(), firstDim.b.data(), 1, 1);
        firstDim.setIsNtt(false);
    }

    // g^-1(a, b) * key
    decompose_rlwe(decFirstDim, firstDim.a, firstDim.b, ellnum, query.getBase(), query.getBg(), modulus);
    for (size_t i = 0; i <  2 * ellnum; i++)
    {
        ntts.ComputeForward(decFirstDim[i].data(), decFirstDim[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.a[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp.data(), length, modulus);

    }

    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.b[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp.data(), length, modulus);
    }

    // return to coeffs form
    ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
    ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
    result.setIsNtt(false);
#endif
}

void addRnsCiphertext(RlweCiphertext& result1, RlweCiphertext& result2,
                        RlweCiphertext& input1, RlweCiphertext& input2)
{
    result1.addAndEqual(input1);
    result2.addAndEqual(input2);
}

/**
 * AutoKeyBSGS class
*/

AutoKeyBSGS::AutoKeyBSGS(/* args */)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    this->ntts = ntts;
#endif  
}

AutoKeyBSGS::AutoKeyBSGS(int32_t len, uint64_t module) : length(len), modulus(module)
{
#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }
#endif
}

AutoKeyBSGS::AutoKeyBSGS(int32_t len, uint64_t module, int32_t elln_bs, int32_t elln_gs) :
         length(len), modulus(module), ellnum_bs(elln_bs), ellnum_gs(elln_gs)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    this->ntts = ntts;
#endif
}

void AutoKeyBSGS::generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                                Secret& secret, const int32_t ellnum, const uint64_t PP,
                                const uint64_t BBg)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp(this->length), autosec(this->length);
    std::vector<uint64_t> mult_factor = powerOfBg(PP, BBg, ellnum);

    automorphic(autosec, secret.getData(), index, modulus);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        intel::hexl::EltwiseFMAMod(temp.data(), autosec.data(),
                mult_factor[i], nullptr, this->length, this->modulus, 1);
        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
#endif
        // rlwe encryption, these rlwes is in ntt form.
        encrypt(result[i].b.data(), result[i].a.data(), secret, temp.data());
        result[i].setIsNtt(true);
    }
    this->isntt = true;
}

/**
 * @brief genreate automorphic keys
 * 
 * @param secret 
 * @param num packing numbers
 */
void AutoKeyBSGS::keyGen(Secret& secret, const std::vector<int32_t> indexLists, StepName stepname)
{
    int32_t ellnum = this->ellnum_bs;
    uint64_t PP = this->PP_bs;
    uint64_t BBg = this->BBg_bs;

    // if the input index lists are for gaint step, use the second parameters
    if (stepname == GaintStep)
    {
        ellnum = this->ellnum_gs;
        PP = this->PP_gs;
        BBg = this->BBg_gs;
    }
    
    std::vector<RlweCiphertext> key_for_index(ellnum, RlweCiphertext(N, modulus));

    for (auto index : indexLists)
    {
        // encrypt automorphic^()(s)
        this->generateSingleKey(key_for_index, index, secret, ellnum, PP, BBg);

        keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
    }
}

uint64_t AutoKeyBSGS::getModulus() const
{
    return this->modulus;
}

uint64_t AutoKeyBSGS::getLength() const
{
    return this->length;
}

uint64_t AutoKeyBSGS::getEllnumBS() const
{
    return this->ellnum_bs;
}

uint64_t AutoKeyBSGS::getBgBS() const
{
    return this->BBg_bs;
}

uint64_t AutoKeyBSGS::getBaseBS() const
{
    return this->PP_bs;
}

uint64_t AutoKeyBSGS::getEllnumGS() const
{
    return this->ellnum_gs;
}

uint64_t AutoKeyBSGS::getBgGS() const
{
    return this->BBg_gs;
}

uint64_t AutoKeyBSGS::getBaseGS() const
{
    return this->PP_gs;
}

bool AutoKeyBSGS::getIsNtt() const
{
    return this->isntt;
}

#ifdef INTEL_HEXL
intel::hexl::NTT AutoKeyBSGS::getNTT() const
{
    return ntts;
}
#endif

/**
 * AutoKeyBSGSRNS class
*/

AutoKeyBSGSRNS::AutoKeyBSGSRNS(/* args */)
{
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(length, modulus);
    this->ntts = ntts;

    intel::hexl::NTT bsNtts(length, bsModulus);
    this->bsNtts = bsNtts;
#endif  
}

AutoKeyBSGSRNS::AutoKeyBSGSRNS(int32_t len, uint64_t module, uint64_t bsModulue) : length(len),
                                modulus(module), bsModulus(bsModulue)
{
#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }

    intel::hexl::NTT bsNtts(length, bsModulus);
    this->bsNtts = bsNtts;    
#endif
}

AutoKeyBSGSRNS::AutoKeyBSGSRNS(int32_t len, uint64_t module, uint64_t bsModule,
         int32_t elln_bs, uint64_t pp_bs, uint64_t bbg_bs) :
         length(len), modulus(module), bsModulus(bsModule),
         ellnum_bs(elln_bs), PP_bs(pp_bs), BBg_bs(bbg_bs)
{
#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }

    if (bsModule == crtBaMod)
    {
        intel::hexl::NTT bsNtts(N, crtBaMod, root_of_unity_crt_bamod);
        this->bsNtts = bsNtts;
    } else {
        intel::hexl::NTT bsNtts(length, bsModulus);
        this->bsNtts = bsNtts;
    }
#endif
}

void AutoKeyBSGSRNS::generateSingleKeyBS(std::vector<RlweCiphertext>& result, int32_t index,
                                Secret& secret, const int32_t ellnum, const uint64_t PP,
                                const uint64_t BBg)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp1(this->length), temp2(this->length);
    std::vector<uint64_t> autosec(this->length), autosec2(this->length);

    // note the powerofBg for 128 bits
    std::vector<uint128_t> mult_factor = powerOfBg128(PP, BBg, ellnum);

    automorphic(autosec, secret.getData(), index, modulus);

    // store in another modulus
    copy(autosec.begin(), autosec.end(), autosec2.begin());
    guass_to_modulus(autosec2.data(), modulus, bsModulus);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        uint64_t temp = mult_factor[i] % this->modulus;
        intel::hexl::EltwiseFMAMod(temp1.data(), autosec.data(),
                temp, nullptr, this->length, this->modulus, 1);

        temp = mult_factor[i] % this->bsModulus;
        intel::hexl::EltwiseFMAMod(temp2.data(), autosec2.data(),
                temp, nullptr, this->length, this->bsModulus, 1);

        // rns encryption, these rlwes is in ntt form.
        encrypt_rns_bsgs_autokey(result[i].b, result[i].a, result[i + ellnum].b, result[i + ellnum].a,
                    secret, temp1, temp2, bsModulus);
        
        result[i].setIsNtt(true);
        result[i + ellnum].setIsNtt(true);
        // result[i + ellnum].setModulus(bsModulus);
#endif
    }
    this->isntt = true;
}

void AutoKeyBSGSRNS::generateSingleKeyGS(std::vector<RlweCiphertext>& result, int32_t index,
                                Secret& secret, const int32_t ellnum, const uint64_t PP,
                                const uint64_t BBg)
{
    // The automorphic transform should be done in coefficient form, instead of ntt form.
    // Also, adding both forms in Secret class can save some computation.
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }

    // compute automorphic^()(s)
    std::vector<uint64_t> temp(this->length), autosec(this->length);
    std::vector<uint64_t> mult_factor = powerOfBg(PP, BBg, ellnum);

    automorphic(autosec, secret.getData(), index, modulus);

    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        intel::hexl::EltwiseFMAMod(temp.data(), autosec.data(),
                mult_factor[i], nullptr, this->length, this->modulus, 1);
        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);
#endif
        // rlwe encryption, these rlwes is in ntt form.
        encrypt(result[i].b.data(), result[i].a.data(), secret, temp.data());
        result[i].setIsNtt(true);
    }
    this->isntt = true;
}

/**
 * @brief genreate automorphic keys
 * 
 * @param secret 
 * @param num packing numbers
 */
void AutoKeyBSGSRNS::keyGen(Secret& secret, const std::vector<int32_t> indexLists, StepName stepname)
{
    // if the input index lists are for gaint step, use the second parameters
    if (stepname == GaintStep)
    {
        std::vector<RlweCiphertext> key_for_index(this->ellnum_gs, RlweCiphertext(N, modulus));

        for (auto index : indexLists)
        {
            // encrypt automorphic^()(s)
            this->generateSingleKeyGS(key_for_index, index, secret, this->ellnum_gs, this->PP_gs, this->BBg_gs);

            keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
        }

        return;
    }

    // there are rns rlwe ciphertexts in baby step
    std::vector<RlweCiphertext> key_for_index(this->ellnum_bs, RlweCiphertext(N, modulus));
    for (size_t i = 0; i < this->ellnum_bs; i++)
    {
        key_for_index.push_back(RlweCiphertext(N, bsModulus));
    }

    for (auto index : indexLists)
    {
        // encrypt automorphic^()(s)
        this->generateSingleKeyBS(key_for_index, index, secret, this->ellnum_bs, this->PP_bs, this->BBg_bs);

        keyMap.insert(std::pair<int32_t, std::vector<RlweCiphertext> >(index, key_for_index));
    }
}

uint64_t AutoKeyBSGSRNS::getModulus() const
{
    return this->modulus;
}

uint64_t AutoKeyBSGSRNS::getBSModulus() const
{
    return this->bsModulus;
}

uint64_t AutoKeyBSGSRNS::getLength() const
{
    return this->length;
}

uint64_t AutoKeyBSGSRNS::getEllnumBS() const
{
    return this->ellnum_bs;
}

uint64_t AutoKeyBSGSRNS::getBgBS() const
{
    return this->BBg_bs;
}

uint64_t AutoKeyBSGSRNS::getBaseBS() const
{
    return this->PP_bs;
}

uint64_t AutoKeyBSGSRNS::getEllnumGS() const
{
    return this->ellnum_gs;
}

uint64_t AutoKeyBSGSRNS::getBgGS() const
{
    return this->BBg_gs;
}

uint64_t AutoKeyBSGSRNS::getBaseGS() const
{
    return this->PP_gs;
}

bool AutoKeyBSGSRNS::getIsNtt() const
{
    return this->isntt;
}

#ifdef INTEL_HEXL
intel::hexl::NTT AutoKeyBSGSRNS::getNTT() const
{
    return ntts;
}

intel::hexl::NTT AutoKeyBSGSRNS::getBSNTT() const
{
    return bsNtts;
}

// used for offline keys only
void AutoKeyBSGSRNS::setBsModules(uint64_t modules)
{
    this->bsModulus = modules;
}

void AutoKeyBSGSRNS::setBSNTT(intel::hexl::NTT& ntts)
{
    this->bsNtts = ntts;
}

#endif
