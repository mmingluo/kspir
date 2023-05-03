#ifndef RECOVER_H
#define RECOVER_H

#include <stdint.h>
#include <stddef.h>
#include <vector>

#include "secret.h"
#include "lwe.h"

void recover(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& new_secret);

#endif
