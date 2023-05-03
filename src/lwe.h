#ifndef LWE_CIPHERTEXT_H
#define LWE_CIPHERTEXT_H

#include <stdint.h>
#include <vector>
#include <map>

#include "params.h"
#include "secret.h"

class LweCiphertext
{
private:
    /* data */
    int32_t length = N;
    uint64_t modulus = bigMod;
    std::vector<uint64_t> a;
    uint64_t b;

public:
    LweCiphertext(/* args */);

    LweCiphertext(int32_t len, uint64_t module);

    std::vector<uint64_t> getA();

    uint64_t getA(int32_t index);

    uint64_t getB();

    uint64_t getModulus();

    int32_t getLength();

    ~LweCiphertext();
};


/**
 * @brief Rlwe Ciphertext class
 * 
 */
class RlweCiphertext
{
private:
    /* data */
    uint64_t length = N;
    uint64_t modulus = bigMod;

    bool isntt = false;

public:
    std::vector<uint64_t> a;
    std::vector<uint64_t> b;

    RlweCiphertext(/* args */);

    RlweCiphertext(uint64_t len, uint64_t module); 

    RlweCiphertext(std::vector<uint64_t>& a, std::vector<uint64_t>& b, bool isn = false);

    // inline void addAndEqual(const RlweCiphertext& cipher);
    void addAndEqual(const RlweCiphertext& cipher);

    // inline void subAndEqual(const RlweCiphertext& cipher);
    void subAndEqual(const RlweCiphertext& cipher);

    /**
     * @brief 
     *  ensure: 0 < number < length
     * 
     * @param number 
     * @return RlweCiphertext& 
     */
    // inline void multConst(int32_t number);
    void multConst(int32_t number);

    void multNegConst(int32_t neg_factor);

    std::vector<uint64_t> getA() const;

    uint64_t getA(int32_t index);

    std::vector<uint64_t> getB() const;

    uint64_t getB(int32_t index);

    uint64_t getModulus() const;

    void setModulus(uint64_t module);

    int32_t getLength() const;

    bool getIsNtt() const;

    void setIsNtt(bool is);
    
    RlweCiphertext& operator=(RlweCiphertext& cipher);

    // void toCoeffForm();

    // void toNttForm();

    // ~RlweCiphertext();
};


class RGSWCiphertext
{
private:
    /* data */
    uint64_t length = N;
    uint64_t modulus = bigMod;

    uint64_t ellnum = 2;
    uint64_t PP = 0x01 << 20; // 18, the base
    uint64_t BBg = 0x01 << 17; // 16

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
#endif

public:
    std::vector<std::vector<uint64_t> > a;
    std::vector<std::vector<uint64_t> > b;

    RGSWCiphertext(/* args */);

    RGSWCiphertext(uint64_t len, uint64_t module); 

    void keyGen(Secret& secret, const uint64_t index, bool is_reverse = false);

/*
    std::vector<uint64_t> getA() const;

    uint64_t getA(int32_t index);

    std::vector<uint64_t> getB() const;

    uint64_t getB(int32_t index);
*/
    uint64_t getModulus() const;

    uint64_t getLength() const;

    uint64_t getEllnum() const;

    uint64_t getBg() const;

    uint64_t getBase() const;

    bool getIsNtt() const;

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT() const;
#endif

    // ~RGSWCiphertext();
};

/**
 * @brief public key for automorphic transform
 * 
 */
class AutoKey
{
private:
    /* data */
    uint64_t length = N;
    uint64_t modulus = bigMod;

    // TODO: change parameter
    int32_t ellnum = 4; // 3; // 4
    uint64_t PP = 0; // 13;// 0x01 << 11; // the base
    uint64_t BBg = 0x01 << 13;

    // std::vector<int32_t> indexLists;

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
#endif

public:
    // std::vector<uint64_t> a;
    // std::vector<uint64_t> b;

    // the map between index and RLWE ciphertext vector
	std::map<int32_t, std::vector<RlweCiphertext> > keyMap;

    AutoKey(/* args */);

    AutoKey(int32_t len, uint64_t module, int32_t elln);

    void generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                    int32_t num, Secret& secret);

    void keyGen(Secret& secret, const int32_t num);

    uint64_t getModulus() const;

    uint64_t getLength() const;

    uint64_t getEllnum() const;

    uint64_t getBg() const;

    uint64_t getBase() const;

    bool getIsNtt() const;

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT() const;
#endif
    
    // ~RlweCiphertext();
};

/**
 * @brief public key for automorphic transform
 * 
 */
class AutoKeyRNS
{
private:
    uint64_t length = N;
    uint64_t modulus1 = bigMod;
    uint64_t modulus2 = bigMod2;

    int32_t ellnum = 9;
    uint64_t PP = 0;// 0x01 << 11; // the base
    uint16_t BBg = 0x01 << 12;

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1;
    intel::hexl::NTT ntts2;
#endif

public:
    // the map between index and RLWE ciphertext vector
	std::map<int32_t, std::vector<RlweCiphertext> > keyMap;

    AutoKeyRNS(/* args */);

    AutoKeyRNS(int32_t len, uint64_t module1, uint64_t module2, int32_t elln);

    void generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                    int32_t num, Secret& secret);

    void keyGen(Secret& secret, const int32_t num);

    uint64_t getModulus1() const;

    uint64_t getModulus2() const;

    uint64_t getLength() const;

    uint64_t getEllnum() const;

    uint64_t getBg() const;

    uint64_t getBase() const;

    bool getIsNtt() const;

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT1() const;
    intel::hexl::NTT getNTT2() const;
#endif
};

void externalProduct(RlweCiphertext& result, RlweCiphertext& firstDim, const RGSWCiphertext& query);

void addRnsCiphertext(RlweCiphertext& result1, RlweCiphertext& result2,
                        RlweCiphertext& input1, RlweCiphertext& input2);

#endif
