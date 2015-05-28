
#include "RSA.h"

RSA::RSA(unsigned int bitSize) {
    BigNum p = BigNum::getPrime(bitSize);
    BigNum q = BigNum::getPrime(bitSize);
    modulus = p * q;
    BigNum eulerFunc = (p - BigNum::ONE) * (q - BigNum::ONE);
    publicKey = (BigNum::ONE << 16) + BigNum::ONE;
    privateKey = BigNum::ZERO;
    BigNum temp = BigNum::ZERO;
    publicKey.binaryGcd(eulerFunc, privateKey, temp);

}

BigNum* RSA::encrypt(string message, int& cyphersSize) {
    cyphersSize = message.size();
    BigNum* encryptedMessage = new BigNum[cyphersSize];
    for (int i = 0; i < cyphersSize; ++i) {
        encryptedMessage[i] = (BigNum::valueOf(message[i])).modPow(publicKey, modulus);
    }
    return encryptedMessage;
}

string RSA::decrypt(BigNum* cyphers, int cyphersSize) {
    string message(cyphersSize, 0);
    for (int i = 0; i < cyphersSize; ++i) {
        message[i] = cyphers[i].modPow(privateKey, modulus).getSmallNum();
    }
    return message;
}