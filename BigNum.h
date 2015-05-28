#ifndef BIGNUM_H
#define BIGNUM_H
#include <algorithm>
#include <cmath>
#include <iostream>
#include <climits>
#include <cstring>
#include <stack>
#include <vector>
#include <memory>
#include <bitset>
#include <fstream>
#include <cstdio>
#include <utility>
#include <cstdlib>
#include <ctime>
using namespace std;

const int NEGATIVE = -1;
const int POSITIVE = 1;

class BigNum
{
public:
    BigNum ();
    BigNum (const unsigned int& bits);
    BigNum (const BigNum& bigNum);
    BigNum (BigNum&& bigNum);
    ~BigNum ();
    bool equals(const BigNum& bigNum) const;
    int compareTo(const BigNum& bigNum) const;
    int signum();

    unsigned int bitLength();
    long getSmallNum();

    BigNum& operator= (BigNum b);
    BigNum operator- ();
    BigNum operator+ (const BigNum& b);
    BigNum operator- (const BigNum& b);
    BigNum operator* (const BigNum& b) const;
    BigNum operator* (int val) const;
    BigNum operator/ (int val);
    BigNum operator/ (const BigNum& b);
    BigNum operator% (const BigNum& b);
    BigNum operator& (const BigNum& b);
    BigNum operator>> (int n);
    BigNum operator<< (int n);

    bool isEven();

    static BigNum valueOf(int val);
    static BigNum ZERO;
    static BigNum ONE;
    static BigNum TWO;



    BigNum binaryGcd (const BigNum& bigNum, BigNum& _x, BigNum& _y) const;

    BigNum modPow (const BigNum &power, const BigNum &mod);
    bool bitAt(unsigned int index) const;
    void setRandom();
    bool isPrime();
    static BigNum getPrime(unsigned int bitSize);
    void print() const;

    void print(ofstream& f) const;

    void read(ifstream& f);

private:
    unsigned long* _num;

    int _signum;

    unsigned int _bits;
    size_t _size;
    const int _baseBits = 32;
    const long _base = 4294967296;
    shared_ptr<BigNum> add (const BigNum& bigNum) const;



    shared_ptr<BigNum> naiveMultiplication(const BigNum &bigNum) const;

    shared_ptr<BigNum> splitToomCook3(const int&start, const int&end) const;
    shared_ptr<BigNum> toomCook3(const BigNum& bigNum) const;
    shared_ptr<BigNum> subtract (const BigNum& bigNum) const;
    shared_ptr<BigNum> multiply (unsigned int val) const;
    shared_ptr<BigNum> multiply (const BigNum& bigNum) const;
    shared_ptr<BigNum>* divideAndRemainder (unsigned int val);
    shared_ptr<BigNum>* divideAndRemainder (const BigNum& bigNum);
    shared_ptr<BigNum> getBits (unsigned int start, unsigned int length) const;
    void breakIntoParts(pair<shared_ptr<BigNum>, int>* C, int&sizeC, shared_ptr<BigNum>* U, int&sizeU, shared_ptr<BigNum>* V, int&sizeV, int r, int q, int p) const;


    void recurse(pair<shared_ptr<BigNum>, int>* C, int&sizeC, shared_ptr<BigNum>* U, int&sizeU, shared_ptr<BigNum>* V, int&sizeV, int r) const;
    void findCoefs(shared_ptr<BigNum>* W, int&sizeW, int r) const;
    shared_ptr<BigNum> setAnswer(shared_ptr<BigNum>* W, int&sizeW, int r, int q, int resultBits) const;
    shared_ptr<BigNum> extend (unsigned int bits) const;

    void halveT3(shared_ptr<BigNum>* t, const BigNum& x, const BigNum& y) const;
    shared_ptr<BigNum> montgomeryProduct(const BigNum& c, const BigNum& d, const BigNum& mod, const BigNum& expressMod, int s, const BigNum& Rminus1);
    bool MillerRabin(const int &r);
    void setCorrectBitSize ();
    shared_ptr<BigNum> abs () const;

    shared_ptr<BigNum> negate () const;
    friend void swap(BigNum& first, BigNum& second);
};


#endif
