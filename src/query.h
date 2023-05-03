#ifndef QUERY_H
#define QUERY_H

#include <stdint.h>
#include <vector>

#include "secret.h"
#include "lwe.h"


void query(RlweCiphertext& cipher, Secret& queryKey, uint64_t row);

// the keyswitch key has size of [2N][N]
void build_keyswitch_key(uint64_t** ks, Secret& secret, std::vector<Secret>& encrypted_keys);

#endif
