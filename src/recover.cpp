#include <math.h>

#include "recover.h"
#include "encrypt.h"

void recover(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& new_secret)
{
    decrypt(message, cipher, new_secret);

    uint64_t modulus = cipher.getModulus();
    uint64_t halfMolulus = modulus >> 1;

    int64_t temp;
    for (auto iter = message.begin(); iter != message.end(); iter++)
    {
        // TODO: check it
        if (*iter > halfMolulus)
            temp = *iter - modulus;
        else
            temp = *iter;
        temp = (int64_t)roundl((long double)temp / Delta);

        if (temp < 0)
            temp += (0x01 << Pbits);
        *iter = temp;
    }
}
