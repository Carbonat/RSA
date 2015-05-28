#ifndef TEST_RSA_H
#define TEST_RSA_H


#include "BigNum.h"

class RSA {
private:
    BigNum privateKey;
    BigNum publicKey;
    BigNum modulus;

public:
    RSA(unsigned int bitSize);
    BigNum* encrypt(string message, int& ciphersSize);
    string decrypt(BigNum* cyphers, int ciphersSize);
};


#endif //TEST_RSA_H
