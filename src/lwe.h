#ifndef LWE_CIPHERTEXT_H
#define LWE_CIPHERTEXT_H

#include <stdint.h>
#include <vector>
#include <map>

#include "params.h"
#include "secret.h"
#include "crt.h"

#ifdef INTEL_HEXL
    #include <hexl/hexl.hpp>
#endif

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

    RGSWCiphertext(uint64_t len, uint64_t module, uint64_t elln, uint64_t base, uint64_t bg); 

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

    AutoKey(int32_t len, uint64_t module, int32_t elln, uint64_t base, uint64_t bg);

    void generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index, 
                    Secret& secret);

    void generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                    int32_t num, Secret& secret);

    void keyGen(Secret& secret, const int32_t num, bool packing_rlwe = false);

    void keyGen(Secret& secret, const std::vector<int32_t> indexLists);

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


enum StepName {BabyStep, GaintStep};

/**
 * @brief public key for automorphic transform in bsgs algorithm
 *        the difference from AutoKey is that we have two ellnums, where the one is for baby step
 *        and the other one is for giant step
 * 
 */
class AutoKeyBSGS
{
private:
    /* data */
    uint64_t length = N;
    uint64_t modulus = bigMod;

    // parameters for baby step
    // TODO: change parameter
    int32_t ellnum_bs = 8; // 5; // 3; // 4
    uint64_t PP_bs = 0; // 0x01 << 8; // 0x01 << 10; // 13;// 0x01 << 11; // the base
    uint64_t BBg_bs = 0x01 << 7; // 0x01 << 7; // 0x01 << 8;

    // std::vector<int32_t> indexLists;

    // parameters for giant step
    int32_t ellnum_gs = 3; // 3; // 6; // 3;
    uint64_t PP_gs = 0x01 << 11; // 0x01 << 11; // 0x01 << 8; // 0x01 << 12;
    uint64_t BBg_gs = 0x01 << 15; // 0x01 << 7; // 0x01 << 14;

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
#endif

public:
    // the map between index and RLWE ciphertext vector
	std::map<int32_t, std::vector<RlweCiphertext> > keyMap;

    AutoKeyBSGS(/* args */);

    AutoKeyBSGS(int32_t len, uint64_t module);

    AutoKeyBSGS(int32_t len, uint64_t module, int32_t elln_bs, int32_t elln_gs);

    void generateSingleKey(std::vector<RlweCiphertext>& result, int32_t index,
                           Secret& secret, const int32_t ellnum, const uint64_t PP,
                           const uint64_t BBg);

    void keyGen(Secret& secret, const std::vector<int32_t> indexLists, StepName stepname = BabyStep);

    uint64_t getModulus() const;

    uint64_t getLength() const;

    uint64_t getEllnumBS() const;

    uint64_t getBgBS() const;

    uint64_t getBaseBS() const;

    uint64_t getEllnumGS() const;

    uint64_t getBgGS() const;

    uint64_t getBaseGS() const;

    bool getIsNtt() const;

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT() const;
#endif
    
    // ~RlweCiphertext();
};


/**
 * @brief public key for automorphic transform in bsgs algorithm
 *        the difference from AutoKeyBSGS is that we use RNS for the baby step.
 *        And then we do a modswitch.
 * 
 */
class AutoKeyBSGSRNS
{
private:
    /* data */
    uint64_t length = N;
    uint64_t modulus = crtMod;
    uint64_t bsModulus = bsMod; // this modulu is for baby step

    // parameters for baby step
    // TODO: change parameter
    int32_t ellnum_bs = 3; // 5; // 3; // 4
    uint64_t PP_bs = 0x01 << 20; // 0x01 << 8; // 0x01 << 10; // 13;// 0x01 << 11; // the base
    uint64_t BBg_bs = 0x01 << 20; // 0x01 << 7; // 0x01 << 8;

    // std::vector<int32_t> indexLists;

    // parameters for giant step
    // int32_t ellnum_gs = 3;
    // uint64_t PP_gs = 0x01 << 11;
    // uint64_t BBg_gs = 0x01 << 15;

    int32_t ellnum_gs = 2;
    uint64_t PP_gs = 0x01 << 20;
    uint64_t BBg_gs = 0x01 << 18;

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
    intel::hexl::NTT bsNtts;
#endif

public:
    // the map between index and RLWE ciphertext vector
	std::map<int32_t, std::vector<RlweCiphertext> > keyMap;

    AutoKeyBSGSRNS(/* args */);

    AutoKeyBSGSRNS(int32_t len, uint64_t module, uint64_t bsModule);

    AutoKeyBSGSRNS(int32_t len, uint64_t module, uint64_t bsModule, int32_t elln_bs, uint64_t pp_bs, uint64_t bbg_bs);

    void generateSingleKeyBS(std::vector<RlweCiphertext>& result, int32_t index,
                           Secret& secret, const int32_t ellnum, const uint64_t PP,
                           const uint64_t BBg);

    void generateSingleKeyGS(std::vector<RlweCiphertext>& result, int32_t index,
                            Secret& secret, const int32_t ellnum, const uint64_t PP,
                            const uint64_t BBg);

    void keyGen(Secret& secret, const std::vector<int32_t> indexLists, StepName stepname = BabyStep);

    void bsgsKeyGen(Secret& secret, int32_t N1);

    uint64_t getModulus() const;
    
    uint64_t getBSModulus() const;

    uint64_t getLength() const;

    uint64_t getEllnumBS() const;

    uint64_t getBgBS() const;

    uint64_t getBaseBS() const;

    uint64_t getEllnumGS() const;

    uint64_t getBgGS() const;

    uint64_t getBaseGS() const;

    bool getIsNtt() const;

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT() const;
    intel::hexl::NTT getBSNTT() const;

    // used for offline keys only
    void setBsModules(uint64_t modules);
    void setBSNTT(intel::hexl::NTT& ntts);
#endif
    
    // ~RlweCiphertext();
};

#endif
