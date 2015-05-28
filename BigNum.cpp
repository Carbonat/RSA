#include <tgmath.h>
#include "BigNum.h"

BigNum BigNum::ZERO = BigNum::valueOf(0);
BigNum BigNum::ONE = BigNum::valueOf(1);
BigNum BigNum::TWO = BigNum::valueOf(2);

BigNum::BigNum() : BigNum(BigNum::ONE) { }

BigNum::BigNum(const unsigned int &bits) {
    _bits = bits;
    if (bits != 0) {
        _size = ceil(((double) bits) / _baseBits);
    }
    else {
        _size = 1;
    }
    _signum = 0;
    _num = new unsigned long[_size];
    memset(_num, 0, _size * sizeof(_num[0]));
}

BigNum::BigNum(const BigNum& bigNum) {
    _bits = bigNum._bits;
    _size = bigNum._size;
    _signum = bigNum._signum;
    _num = new unsigned long[_size];
    copy(bigNum._num, bigNum._num + bigNum._size, _num);
}

BigNum::BigNum(BigNum&& bigNum)
        : BigNum(bigNum._bits) {
    swap(*this, bigNum);
}

BigNum::~BigNum() {
    if (_size != 0) {
        delete[] _num;
    }
    _bits = 0;
    _size = 0;
    _signum = 0;
}

bool BigNum::equals(const BigNum& bigNum) const {
    int res = this->compareTo(bigNum);
    return this->compareTo(bigNum) == 0 ? true : false;
}

int BigNum::compareTo(const BigNum &bigNum) const {
    if (_signum != bigNum._signum) {
        return _signum > bigNum._signum ? 1 : -1;
    }
    if (_size != bigNum._size)
        return _size > bigNum._size ? 1 : -1;
    for (int i = _size - 1; i >= 0; --i) {
        if (_num[i] != bigNum._num[i])
            return _num[i] > bigNum._num[i] ? 1 : -1;
    }
    return 0;
}

int BigNum::signum() {
    return _signum;
}

unsigned int BigNum::bitLength() {
    return _bits;
}

long BigNum::getSmallNum() {
    return _signum * _num[0];
}

BigNum& BigNum::operator=(BigNum b) {
    swap(*this, b);
    return *this;
}

BigNum BigNum::operator-() {
    return *(this->negate());
}

BigNum BigNum::operator+(const BigNum &b) {
    if (this->_signum == 0)
        return BigNum(b);
    if (b._signum == 0)
        return BigNum(*this);

    shared_ptr<BigNum> result;

    shared_ptr<BigNum> first = this->abs();
    shared_ptr<BigNum> second = b.abs();

    if (this->_signum * b._signum == -1) {
        int compare = first->compareTo(*second);

        if (compare == 1) {
            result = first->subtract(*second);
            if (this->_signum == NEGATIVE)
                result->_signum = NEGATIVE;
            else // b._signum == NEGATIVE
                result->_signum = POSITIVE;
        }
        else if (compare == 0)
            result = make_shared<BigNum>(0);
        else {
            result = second->subtract(*first);
            if (this->_signum == NEGATIVE)
                result->_signum = POSITIVE;
            else // b._sinum == NEGATIVE
                result->_signum = NEGATIVE;
        }

        return *result;
    }

    result = first->add(*second);
    if (this->_signum == POSITIVE && b._signum == POSITIVE)
        result->_signum = POSITIVE;
    else // this->_signum == NEGATIVE && b._signum == NEGATIVE
        result->_signum = NEGATIVE;

    return *result;
}

BigNum BigNum::operator-(const BigNum &b) {
    if (b._bits == 0 || b._signum == 0)
        return BigNum(*this);

    if (this->_signum == 0)
        return *(b.negate());

    shared_ptr<BigNum> result;

    shared_ptr<BigNum> first = this->abs();
    shared_ptr<BigNum> second = b.abs();

    if (this->_signum * b._signum == -1) {
        result = first->add(*second);
        if (this->_signum == NEGATIVE)
            result->_signum = NEGATIVE;
        else
            result->_signum = POSITIVE;
    }
    else {
        int compare = first->compareTo(*second);
        if (compare == 1) {
            result = first->subtract(*second);
            if (this->_signum == POSITIVE && b._signum == POSITIVE)
                result->_signum = POSITIVE;
            else // this._signum == NEGATIVE && b._signum == NEGATIVE
                result->_signum = NEGATIVE;

        }
        else if (compare == 0) {
            return BigNum(0);
        }
        else {
            result = second->subtract(*first);
            if (this->_signum == POSITIVE && b._signum == POSITIVE)
                result->_signum = NEGATIVE;
            else // this._signum == NEGATIVE && b._signum == NEGATIVE
                result->_signum = POSITIVE;
        }
    }
    return *result;
}

BigNum BigNum::operator*(const BigNum &b) const {
    if (this->equals(BigNum::ZERO) || b.equals(BigNum::ZERO)) {
        return BigNum::ZERO;
    }
    shared_ptr<BigNum> result;
    if (this->_size < 150 || b._size < 150) {
        result = this->naiveMultiplication(b);
    }
    else {
        result = this->toomCook3(b);
    }
    if (this->_signum == NEGATIVE ^ b._signum == NEGATIVE)
        result->_signum = NEGATIVE;
    else
        result->_signum = POSITIVE;
    return *result;
}

BigNum BigNum::operator* (int val) const {
    if (this->equals(BigNum::ZERO) || val == 0) {
        return BigNum::ZERO;
    }
    unsigned int absVal = std::abs(val);
    shared_ptr<BigNum> result = this->multiply(absVal);
    if (this->_signum == NEGATIVE ^ val < 0)
        result->_signum = NEGATIVE;
    else
        result->_signum = POSITIVE;
    return *result;
}

BigNum BigNum::operator/ (int val) {
    if (this->equals(BigNum::ZERO) || val == 0) {
        return BigNum::ZERO;
    }
    shared_ptr<BigNum> *quotientAndRemainder;
    shared_ptr<BigNum> quotient;
    shared_ptr<BigNum> dividend = this->abs();
    unsigned int divisor = std::abs(val);
    if (this->_size == 1 && this->_num[0] < divisor) {
        quotient = make_shared<BigNum>(BigNum::ZERO);
    }
    else {
        quotientAndRemainder = dividend->divideAndRemainder(divisor);
        quotient = quotientAndRemainder[0];
        delete[] quotientAndRemainder;
        if (this->_signum == NEGATIVE ^ val < 0)
            quotient->_signum = NEGATIVE;
        else
            quotient->_signum = POSITIVE;
    }
    return *quotient;
}

BigNum BigNum::operator/(const BigNum &b) {
    if (b.equals(BigNum::ZERO)) {
        return BigNum::ZERO;
    }
    shared_ptr<BigNum> *quotientAndRemainder;
    shared_ptr<BigNum> quotient;
    shared_ptr<BigNum> dividend = this->abs();
    shared_ptr<BigNum> divisor = b.abs();
    int compare = dividend->compareTo(*divisor);
    if (compare >= 0) {
        quotientAndRemainder = dividend->divideAndRemainder(*divisor);
        quotient = quotientAndRemainder[0];
        delete[] quotientAndRemainder;
        if (this->_signum == NEGATIVE ^ b._signum == NEGATIVE)
            quotient->_signum = NEGATIVE;
        else
            quotient->_signum = POSITIVE;
    }
    else {
        quotient = make_shared<BigNum>(BigNum::ZERO);
    }
    return *quotient;
}

BigNum BigNum::operator%(const BigNum &b) {
    shared_ptr<BigNum> *quotientAndRemainder;
    shared_ptr<BigNum> quotient, remainder;
    shared_ptr<BigNum> dividend = this->abs();
    shared_ptr<BigNum> divisor = b.abs();
    int compare = dividend->compareTo(*divisor);

    if (compare == 0)
        remainder = make_shared<BigNum>(0);
    else {
        if (this->_signum * b._signum == 1) {
            if (compare == 1) {
                quotientAndRemainder = dividend->divideAndRemainder(*divisor);
                remainder = quotientAndRemainder[1];
                delete[] quotientAndRemainder;
            }
            else // compare == -1
                remainder = make_shared<BigNum>(*dividend);

            if (this->_signum == POSITIVE && b._signum == POSITIVE)
                remainder->_signum = POSITIVE;
            else // this->_signum == NEGATIVE && b._signum == NEGATIVE
                remainder->_signum = NEGATIVE;
        }
        else {
            if (compare == 1) {
                quotientAndRemainder = dividend->divideAndRemainder(*divisor);
                quotient = quotientAndRemainder[0];
                delete[] quotientAndRemainder;

                BigNum one(BigNum::valueOf(1));
                one._signum = 0;                                                // attention
                shared_ptr<BigNum> temp = (quotient->add(one));
                temp = temp->multiply(*divisor);                                // attention
                remainder = temp->subtract(*dividend);
                if (this->_signum == NEGATIVE)
                    remainder->_signum = POSITIVE;
                else  // b->_signum == NEGATIVE
                    remainder->_signum = NEGATIVE;
            }
            else { // compare == -1
                remainder = divisor->subtract(*dividend);
                if (this->_signum == NEGATIVE)
                    remainder->_signum = POSITIVE;
                else // b._signum == NEGATIVE
                    remainder->_signum = NEGATIVE;
            }
        }
    }

    return *remainder;
}


BigNum BigNum::operator&(const BigNum &b) {
    if (this->_signum == 0 || b._signum == 0)
        return BigNum(0);
    shared_ptr<BigNum> first, second, result;
    if (this->_bits > b._bits) {
        first = make_shared<BigNum>(*this);
        second = b.extend(this->_bits);
        result = make_shared<BigNum>(this->_bits);
    }
    else {
        first = this->extend(b._bits);
        second = make_shared<BigNum>(b);
        result = make_shared<BigNum>(b._bits);
    }
    for (int i = 0; i < result->_size; ++i) {
        result->_num[i] = first->_num[i] & second->_num[i];
    }
    result->_signum = POSITIVE;
    return *result;
}


BigNum BigNum::operator>>(int n) {
    if (n >= _bits || _signum == 0) {
        return BigNum::ZERO;
    }

    shared_ptr<BigNum> result;

    int bigShift = n / _baseBits;
    int smallShift = n % _baseBits;

    if (bigShift != 0) {
        result = make_shared<BigNum>(_bits - bigShift * _baseBits);
        for (int j = 0; j + bigShift < _size; ++j) {
            result->_num[j] = this->_num[j + bigShift];
        }
    }
    else {
        result = make_shared<BigNum>(*this);
    }

    if (smallShift != 0) {
        unsigned int prevK = 0, k = 0;
        unsigned int bitsSelector = (1 << smallShift) - 1;
        for (int i = result->_size - 1; i >= 0; --i) {
            k = result->_num[i] & bitsSelector;
            result->_num[i] >>= smallShift;
            result->_num[i] |= prevK;
            prevK = k << (result->_baseBits - smallShift);
        }
    }
    result->_signum = this->_signum;
    result->setCorrectBitSize();
    return *result;
}

BigNum BigNum::operator<<(int n) {
    if (_signum == 0) {
        return BigNum::ZERO;
    }
    int bigShift = n / _baseBits;
    int smallShift = n % _baseBits;

    shared_ptr<BigNum> result = make_shared<BigNum>(_bits + n);
    if (bigShift != 0) {
        for (int i = 0; i < this->_size; ++i) {
            result->_num[i + bigShift] = this->_num[i];
        }
    }
    else {
        for (int i = 0; i < this->_size; ++i) {
            result->_num[i] = this->_num[i];
        }
    }

    if (smallShift != 0) {
        unsigned int prevK = 0, k = 0;
        unsigned int bitsSelector = UINT_MAX << (_baseBits - smallShift);
        unsigned int bitsDelete = UINT_MAX >> smallShift;
        for (int i = 0; i < result->_size; ++i) {
            k = result->_num[i] & bitsSelector;
            k >>= (result->_baseBits - smallShift);
            result->_num[i] &= bitsDelete;
            result->_num[i] <<= smallShift;
            result->_num[i] |= prevK;
            prevK = k;
        }
    }
    result->_signum = this->_signum;
    result->setCorrectBitSize();
    return *result;
}

shared_ptr<BigNum> BigNum::naiveMultiplication(const BigNum &bigNum) const {
    unsigned int m = this->_size;
    unsigned int n = bigNum._size;
    int i = 0;
    unsigned long k = 0;
    unsigned long t = 0;
    shared_ptr<BigNum> w = make_shared<BigNum>((m + n + 1) * this->_baseBits);
    for (int j = 0; j < bigNum._size; ++j) {
        if(bigNum._num[j] == 0) {
            w->_num[j + m] = 0;
        }
        else {
            for (k = 0, i = 0; i < m; ++i) {
                t = this->_num[i] * bigNum._num[j] + w->_num[i + j] + k;
                w->_num[i + j] = t % w->_base;
                k = t >> w->_baseBits;
            }
            w->_num[j + m] = k;
        }
    }

    w->setCorrectBitSize();
    return w;
}

bool BigNum::isEven() {
    return _num[0] & 1 == 1 ? false : true;
}


BigNum BigNum::valueOf(int val) {
    if (val == 0)
        return BigNum(0);

    shared_ptr<BigNum> result = make_shared<BigNum>(32);        // magic number
    result->_num[0] = std::abs(val);

    if (val > 0)
        result->_signum = POSITIVE;
    else
        result->_signum = NEGATIVE;

    result->setCorrectBitSize();

    return *result;
}

shared_ptr<BigNum> BigNum::abs() const {
    shared_ptr<BigNum> result = make_shared<BigNum>(*this);
    if (result->_signum == NEGATIVE) {
        result->_signum = POSITIVE;
    }
    return result;
}


shared_ptr<BigNum> BigNum::negate() const {
    shared_ptr<BigNum> result = make_shared<BigNum>(*this);
    result->_signum = -(result->_signum);
    return result;
}


// change types of k.
shared_ptr<BigNum> BigNum::add(const BigNum &bigNum) const {
    shared_ptr<BigNum> result;
    shared_ptr<BigNum> augend;
    shared_ptr<BigNum> addend;
    int compare = this->compareTo(bigNum);
    if (compare == 1) {
        augend = make_shared<BigNum>(*this);
        addend = make_shared<BigNum>(bigNum);
    }
    else {
        augend = make_shared<BigNum>(bigNum);
        addend = make_shared<BigNum>(*this);
    }
    result = make_shared<BigNum>(augend->_bits);

    int k = 0;
    size_t j = 0;
    while (j < addend->_size) {
        result->_num[j] = (augend->_num[j] + addend->_num[j] + k) % _base;
        k = (augend->_num[j] + addend->_num[j] + k) / _base;
        ++j;
    }

    while (j < augend->_size) {
        result->_num[j] = (augend->_num[j] + k) % _base;
        k = (augend->_num[j] + k) / _base;
        ++j;
    }

    if (k != 0) {
        shared_ptr<BigNum> formatResult = make_shared<BigNum>(result->_bits + _baseBits);
        size_t i = 0;
        while (i < result->_size) {
            formatResult->_num[i] = result->_num[i];
            ++i;
        }
        formatResult->_num[i] = k;
        formatResult->setCorrectBitSize();
        return formatResult;
    }

    result->setCorrectBitSize();
    return result;
}

shared_ptr<BigNum> BigNum::subtract(const BigNum &bigNum) const {
    shared_ptr<BigNum> result = make_shared<BigNum>(_bits);
    long k = 0;
    size_t j = 0;
    long difference = 0;
    while (j < bigNum._size) {
        difference = _num[j] - bigNum._num[j] + k;

        if (difference < 0) {
            result->_num[j] = (_base + difference) % _base;
            k = -1;
        }
        else {
            result->_num[j] = difference % _base;
            k = 0;
        }
        ++j;
    }
    while (j < _size) {
        difference = _num[j] + k;

        if (difference < 0) {
            result->_num[j] = (_base + difference) % _base;
            k = -1;
        }
        else {
            result->_num[j] = difference % _base;
            k = 0;
        }
        ++j;
    }

    result->setCorrectBitSize();
    return result;
}


// delete k
void BigNum::breakIntoParts(pair<shared_ptr<BigNum>, int>* C, int& sizeC,
                            shared_ptr<BigNum>* U, int&sizeU,
                            shared_ptr<BigNum>* V, int&sizeV,
                            int r, int q, int p) const {
    shared_ptr<BigNum> u = C[--sizeC].first;
    shared_ptr<BigNum> v = C[--sizeC].first;

    shared_ptr<BigNum> u_split[r + 1];
    shared_ptr<BigNum> v_split[r + 1];
    int nextStart = 0;
    for (int i = 0; i <= r; ++i) {
        nextStart = i * q;
        u_split[i] = u->getBits(nextStart, q);
        v_split[i] = v->getBits(nextStart, q);
    }

    int doubleR = r * 2;
    U[sizeU++] = u_split[0];
    V[sizeV++] = v_split[0];
    for (int j = 1; j <= doubleR; ++j) {
        shared_ptr<BigNum> sum_v = make_shared<BigNum>(p);
        shared_ptr<BigNum> sum_u = make_shared<BigNum>(p);
        for (int i = r; i > 0; --i) {
            sum_u = (sum_u->add(*(u_split[i])))->multiply(j);
            sum_v = (sum_v->add(*(v_split[i])))->multiply(j);
        }
        sum_u = sum_u->add(*(u_split[0]));
        sum_v = sum_v->add(*(v_split[0]));
        U[sizeU++] = sum_u;
        V[sizeV++] = sum_v;
    }
}

void BigNum::recurse(pair<shared_ptr<BigNum>, int>* C, int&sizeC,
                     shared_ptr<BigNum>* U, int&sizeU,
                     shared_ptr<BigNum>* V, int&sizeV,
                     int r) const {
    C[sizeC++] = make_pair(nullptr, 2);
    C[sizeC++] = make_pair(V[--sizeV], 0);
    C[sizeC++] = make_pair(U[--sizeU], 0);
    for (int i = r * 2 - 1; i >= 0; --i) {
        C[sizeC++] = make_pair(nullptr, 3);
        C[sizeC++] = make_pair(V[--sizeV], 0);
        C[sizeC++] = make_pair(U[--sizeU], 0);
    }
}


void BigNum::findCoefs(shared_ptr<BigNum>* W, int&sizeW, int r) const {
    int doubleR = 2 * r;
    int indexW0 = sizeW - doubleR - 1;

    for (int j = 1; j <= doubleR; ++j) {
        for (int t = doubleR; t >= j; --t) {
            shared_ptr<BigNum> *res = (W[indexW0 + t]->subtract(*(W[indexW0 + t - 1])))->divideAndRemainder(j);
            W[indexW0 + t] = res[0];
            delete[] res;
        }
    }

    for (int j = doubleR - 1; j > 0; --j) {
        for (int t = j; t < doubleR; ++t) {
            W[indexW0 + t] = (W[indexW0 + t]->subtract(*(W[indexW0 + t + 1]->multiply(j))));
        }
    }
}

shared_ptr<BigNum> BigNum::setAnswer(shared_ptr<BigNum>* W, int&sizeW, int r, int q, int resultBits) const {
    int doubleR = 2 * r;
    unsigned long indexW0 = sizeW - doubleR - 1;
    unsigned int factor = pow(2, q);
    int repeat = 1;
    if (q > 16) {
        repeat = q / 16;
        factor = 1 << 16;
    }


    shared_ptr<BigNum> result = make_shared<BigNum>(resultBits);
    int i = doubleR;
    while (i > 0) {
        result = result->add(*(W[i + indexW0]));
        for (int j = 0; j < repeat; ++j) {
            result = result->multiply(factor);
        }
        --i;
    }
    result = result->add((*W[i + indexW0]));
    sizeW -= doubleR + 1;
    return result;
}

// this and bigNum have the same size
shared_ptr<BigNum> BigNum::multiply(const BigNum &bigNum) const {

    if (this->_bits <= 32 && bigNum._bits <= 32) {
        return multiply(bigNum._num[0]);
    }

    vector<unsigned long long> q_arr(0);
    vector<unsigned long long> r_arr(0);
    q_arr.push_back(16); q_arr.push_back(16);
    r_arr.push_back(4); r_arr.push_back(4);
    unsigned long Q = 4;
    unsigned long R = 2;

    unsigned int sizesW[] = { 9, 17, 25, 41, 57};
    unsigned int sizesC[] = { 28, 53, 78, 12, 176};
    unsigned int sizesUV[] = { 9, 9, 9, 17, 17};

    int k = 1;
    int maxBits = _bits > bigNum._bits ? _bits : bigNum._bits;
    while (q_arr[k - 1] + q_arr[k] < maxBits) {
        Q += R;
        R = floor(sqrt(Q));
        unsigned long long qtemp = pow(2.0, Q);
        q_arr.push_back(qtemp);
        unsigned long long rtemp = pow(2.0, R);
        r_arr.push_back(rtemp);
        ++k;
        cout << "k = " << k << endl;
        cout << "Q = " << Q << endl;
        cout << "R = " << R << endl;
        cout << "q = " << q_arr.back() << endl;
        cout << "r = " << r_arr.back() << endl << endl;        // delete
    }
    unsigned int r = 0, q = 0;
    unsigned int p = q_arr[k - 1] + q_arr[k];
    cout << "p = " << p << endl;

    pair<shared_ptr<BigNum>, int >* C = new pair<shared_ptr<BigNum>, int> [176];
    shared_ptr<BigNum>* U = new shared_ptr<BigNum>[17];
    shared_ptr<BigNum>* V = new shared_ptr<BigNum>[17];
    shared_ptr<BigNum>* W = new shared_ptr<BigNum>[57];
    shared_ptr<BigNum> w;
    int sizeU = 0, sizeV = 0, sizeW = 0, sizeC = 0;
    C[sizeC++] = make_pair(nullptr, 1);
    C[sizeC++] = make_pair(this->extend(p), 0);
    C[sizeC++] = make_pair(bigNum.extend(p), 0);

    int code = 0;

    clock_t start = 0, end = 0;
    clock_t bipSum = 0, rSum = 0, multSum = 0, fcSum = 0, saSum = 0;

    while (code != 1) {
        while (k > 0) {
            --k;
            if (k != 0) {
                r = r_arr[k];
                q = q_arr[k];
                p = q_arr[k - 1] + q_arr[k];
                breakIntoParts(C, sizeC, U, sizeU, V, sizeV, r, q, p);                   // delete k
                recurse(C, sizeC, U, sizeU, V, sizeV, r);                             // delete k
                rSum += clock() - start;

            }
            else {
                shared_ptr<BigNum> u = C[--sizeC].first;
                shared_ptr<BigNum> v = C[--sizeC].first;
                w = u->multiply(*v);
                code = 0;
            }
        }
        while (code != 3 && code != 1) {
            ++k;
            code = C[--sizeC].second;
            W[sizeW++] = w;

            if (code == 2) {
                r = r_arr[k];
                q = q_arr[k];
                p = q_arr[k - 1] + q_arr[k];
                findCoefs(W, sizeW, r);
                w = setAnswer(W, sizeW, r, q, 2 * (q_arr[k] + q_arr[k + 1]));
            }
        }
    }
    start = clock();
    w->setCorrectBitSize();
    end = clock();

    delete[] W;
    delete[] C;
    delete[] U;
    delete[] V;
    return w;
}


shared_ptr<BigNum> BigNum::multiply(unsigned int val) const {
    if (this->equals(BigNum::ZERO) || val == 0) {
        return make_shared<BigNum>(BigNum::ZERO);
    }

    if (val == 1) {
        shared_ptr<BigNum> result = make_shared<BigNum>(_bits);
        for (int i = 0; i < _size; ++i) {
            result->_num[i] = _num[i];
        }
        result->setCorrectBitSize();
        return result;
    }
    else {
        shared_ptr<BigNum> result = make_shared<BigNum>(_bits + _baseBits);
        unsigned long k = 0;
        unsigned long product = 0;
        int i = 0;
        while (i < _size) {
            product = _num[i] * val + k;
            result->_num[i] = product % _base;
            k = product / _base;
            ++i;
        }
        result->_num[i] = k;

        result->setCorrectBitSize();
        return result;
    }
}


shared_ptr<BigNum> *BigNum::divideAndRemainder(unsigned int val) {
    shared_ptr<BigNum> *result = new shared_ptr<BigNum>[2];
    if (_bits == 0) {
        result[0] = make_shared<BigNum>(0);
        result[1] = make_shared<BigNum>(0);
        return result;
    }

    shared_ptr<BigNum> quotient = make_shared<BigNum>(_bits);
    shared_ptr<BigNum> remainder = make_shared<BigNum>(_baseBits);
    unsigned long dividend = 0;
    for (int i = _size - 1; i >= 0; --i) {
        dividend = (remainder->_num[0] << _baseBits) + _num[i];
        quotient->_num[i] = dividend / val;
        remainder->_num[0] = dividend % val;
    }

    quotient->setCorrectBitSize();
    result[0] = quotient;
    result[1] = remainder;
    return result;
}


shared_ptr<BigNum> *BigNum::divideAndRemainder(const BigNum &bigNum) {
    shared_ptr<BigNum> *result = new shared_ptr<BigNum>[2];
    if (_bits == 0) {
        result[0] = make_shared<BigNum>(BigNum::ZERO);
        result[1] = make_shared<BigNum>(BigNum::ZERO);
        return result;
    }
    int n = bigNum._size;
    if (n == 1) {
        delete[] result;
        return this->divideAndRemainder(bigNum._num[0]);
    }

    int m = _size - bigNum._size;

    shared_ptr<BigNum> quotient = make_shared<BigNum>((m + 1) * _baseBits);
    shared_ptr<BigNum> remainder = make_shared<BigNum>(n * _baseBits);

    unsigned int d = _base / (bigNum._num[n - 1] + 1);
    shared_ptr<BigNum> u = make_shared<BigNum>(this->operator*(d));
    if (u->_size == _size) {
        unsigned long *temp = u->_num;
        ++(u->_size);
        u->_bits += _baseBits;
        u->_num = new unsigned long[u->_size];
        memset(u->_num, 0, u->_size * sizeof(u->_num[0]));
        for (int i = 0; i < _size; ++i) {
            u->_num[i] = temp[i];
        }
        delete[] temp;
    }

    shared_ptr<BigNum> v = bigNum.multiply(d);


    for (int j = m; j >= 0; --j) {
        unsigned long q = (u->_num[j + n] * _base + u->_num[j + n - 1]) / v->_num[n - 1];
        unsigned long r = (u->_num[j + n] * _base + u->_num[j + n - 1]) % v->_num[n - 1];

        bool qIsCorrect = q * v->_num[n - 2] > _base * r + u->_num[j + n - 1];

        while (q == _base || qIsCorrect) {
            --q;
            r += v->_num[n - 1];
            qIsCorrect = q * v->_num[n - 2] > _base * r + u->_num[j + n - 1];
            if (r >= _base) {
                break;
            }
        }

        if (q != 0) {
            shared_ptr<BigNum> subtrahend = v->multiply(q);
            shared_ptr<BigNum> splitedU = make_shared<BigNum>((n + 1) * _baseBits);
            for (int i = 0; i <= n; ++i) {
                splitedU->_num[i] = u->_num[i + j];
            }

            int compare = splitedU->compareTo(*subtrahend);
            if (compare == 0) {
                for (int i = j + n; i >= j; --i) {
                    u->_num[i] = 0;
                }
            }
            else {
                if (compare == -1) {
                    --q;
                    subtrahend = make_shared<BigNum>(v->operator*(q));
                }
                splitedU = splitedU->subtract(*subtrahend);

                int i = 0;
                while (i <= splitedU->_size) {
                    u->_num[i + j] = splitedU->_num[i];
                    ++i;
                }
                while (i <= n) {
                    u->_num[i + j] = 0;
                    ++i;
                }
            }
        }

        quotient->_num[j] = q;
    }

    shared_ptr<BigNum> preRemainder = make_shared<BigNum>(n * _baseBits);
    for (int i = 0; i < n; ++i) {
        preRemainder->_num[i] = u->_num[i];
    }
    shared_ptr<BigNum> *temp = preRemainder->divideAndRemainder(d);
    remainder = temp[0];

    delete[] temp;
    quotient->setCorrectBitSize();

    result[0] = quotient;
    result[1] = remainder;

    return result;
}


shared_ptr<BigNum> BigNum::getBits(unsigned int start, unsigned int length) const {
    shared_ptr<BigNum> result = make_shared<BigNum>(length);
    unsigned int k = start / _baseBits; // index of this->_num
    int g = start % _baseBits; // necessary first bit in this->_num[k]
    unsigned long bit = 0;
    int i = 0; // index of result->_num
    int j = 0; // current bit in result->_num[i]
    int h = 0; // bits in result
    while (h < length && k < _size && i < result->_size) {
        while (h < length && g < _baseBits && j < _baseBits) {
            bit = _num[k] & (1 << g);
            if (bit != 0) {
                result->_num[i] |= (unsigned int) 1 << j;
            }
            ++h;
            ++j;
            ++g;
        }
        if (j == _baseBits) {
            j = 0;
            ++i;
        }
        if (g == _baseBits) {
            g = 0;
            ++k;
        }
    }
    return result;
}

void BigNum::halveT3(shared_ptr<BigNum>* t, const BigNum& x, const BigNum& y) const {
    if (t[0]->isEven() && t[1]->isEven()) {
        t[0] = make_shared<BigNum>(t[0]->operator>>(1));
        t[1] = make_shared<BigNum>(t[1]->operator>>(1));
        t[2] = make_shared<BigNum>(t[2]->operator>>(1));
    }
    else {
        t[0] = make_shared<BigNum>(((t[0]->operator+(y)).operator>>(1)));
        t[1] = make_shared<BigNum>(((t[1]->operator-(x)).operator>>(1)));
        t[2] = make_shared<BigNum>(t[2]->operator>>(1));
    }
}

BigNum BigNum::binaryGcd (const BigNum& bigNum, BigNum& _x, BigNum& _y) const {
    if (this->_signum == 0) {
        return *bigNum.abs();
    }
    if (bigNum._signum == 0) {
        return *(this->abs());
    }
    BigNum x = *this;
    BigNum y = bigNum;

    unsigned int beta = 0;
    while (x.isEven() && y.isEven()) {
        x = x >> 1;
        y = y >> 1;
        ++beta;
    }

    BigNum a = BigNum::ONE;
    BigNum b = BigNum::ZERO;
    BigNum h = x;

    shared_ptr<BigNum> *v = new shared_ptr<BigNum>[3];
    v[0] = make_shared<BigNum>(y);
    v[1] = make_shared<BigNum>(BigNum::ONE - x);
    v[2] = make_shared<BigNum>(y);
    shared_ptr<BigNum> *t = new shared_ptr<BigNum>[3];

    if (x.isEven()) {
        t[0] = make_shared<BigNum>(BigNum::ONE);
        t[1] = make_shared<BigNum>(BigNum::ZERO);
        t[2] = make_shared<BigNum>(x);
    }
    else {
        t[0] = make_shared<BigNum>(BigNum::ZERO);
        t[1] = make_shared<BigNum>(BigNum::valueOf(-1));
        t[2] = make_shared<BigNum>(-y);
    }
    do {
        while (t[2]->isEven()) {
            halveT3(t, x, y);
        }
        if (t[2]->_signum > 0) {
            a = *t[0];
            b = *t[1];
            h = *t[2];
        }
        else {
            v[0] = make_shared<BigNum>(y - *t[0]);
            v[1] = make_shared<BigNum>((-x) - *t[1]);
            v[2] = make_shared<BigNum>(*(t[2]->negate()));
        }
        t[0] = make_shared<BigNum>(a - (*v[0]));
        t[1] = make_shared<BigNum>(b - (*v[1]));
        t[2] = make_shared<BigNum>(h - (*v[2]));

        if (t[0]->_signum < 0) {
            t[0] = make_shared<BigNum>(t[0]->operator+(y));
            t[1] = make_shared<BigNum>(t[1]->operator-(x));
        }
    }
    while (t[2]->_signum != 0);

    delete[] t;
    delete[] v;

    BigNum result = h << beta;
    _x = a;
    _y = b;
    return result;
}

shared_ptr<BigNum> BigNum::montgomeryProduct(const BigNum& c, const BigNum& d, const BigNum& mod, const BigNum& expressMod, int s, const BigNum& Rminus1) {
    BigNum x = c * d;
    BigNum z = x * expressMod;
    z = z & Rminus1;
    z = z * mod;
    z = z + x;
    z = z >> s;
    if (z.compareTo(mod) >= 0) {
        z = z - mod;
    }
    return make_shared<BigNum>(z);
}

BigNum BigNum::modPow(const BigNum &power, const BigNum &mod) {
    int s = mod._bits;
    BigNum R = BigNum::ONE << s;
    BigNum invMod = BigNum::ZERO;
    BigNum temp = BigNum::ZERO;
    mod.binaryGcd(R, invMod, temp);
    BigNum expressMod = (-invMod) % R;
    BigNum imX = (*this * R) % mod;
    shared_ptr<BigNum> imP = make_shared<BigNum>(R % mod);
    BigNum Rminus1 = R - BigNum::ONE;
    for (int j = power._bits - 1; j >= 0; --j) {
        imP = montgomeryProduct(*imP, *imP, mod, expressMod, s, Rminus1);

        if (power.bitAt(j)) {
            imP = montgomeryProduct(*imP, imX, mod, expressMod, s, Rminus1);
        }
    }
    return *(montgomeryProduct(*imP, BigNum::ONE, mod, expressMod, s, Rminus1));
}

bool BigNum::bitAt(unsigned int index) const {
    shared_ptr<BigNum> temp = this->getBits(index, 1);
    return temp->_num[0] == 1 ? true : false;
}

bool BigNum::isPrime() {
    if (this->isEven() || this->equals(BigNum::ONE)) {
        return false;
    }

    if (    this->operator%(BigNum::valueOf(3)).equals(BigNum::ZERO)  ||
            this->operator%(BigNum::valueOf(5)).equals(BigNum::ZERO)  ||
            this->operator%(BigNum::valueOf(7)).equals(BigNum::ZERO)  ||
            this->operator%(BigNum::valueOf(11)).equals(BigNum::ZERO) ||
            this->operator%(BigNum::valueOf(13)).equals(BigNum::ZERO) ||
            this->operator%(BigNum::valueOf(17)).equals(BigNum::ZERO) ||
            this->operator%(BigNum::valueOf(19)).equals(BigNum::ZERO) ||
            this->operator%(BigNum::valueOf(23)).equals(BigNum::ZERO) ||
            this->operator%(BigNum::valueOf(29)).equals(BigNum::ZERO) )
        return false;
    int rounds = 2;
    return this->MillerRabin(rounds);
}


bool BigNum::MillerRabin(const int &round) {
    if (this->isEven() || this->equals(BigNum::ONE)) {
        return false;
    }

    int s = 0;
    BigNum thisMinus1 = *this - BigNum::ONE;
    BigNum t = thisMinus1;
    while (t.isEven()) {
        t = t >> 1;
        ++s;
    }

    srand(time(NULL));
    int maxABits = (thisMinus1 - BigNum::ONE)._bits;
    int aBits = 0;
    for (int i = 0; i < round; ++i) {
        aBits = rand() % maxABits;
        BigNum a(aBits);
        a.setRandom();
        a.setCorrectBitSize();
        while (a.compareTo(*this) >= 0) {
            a.setRandom();
            a.setCorrectBitSize();
        }
        BigNum x = a.modPow(t, *this);
        int j = 0;
        while (!((j == 0 && x.compareTo(BigNum::ONE) == 0) || x.compareTo(thisMinus1) == 0)) {
            if ((j > 0 && x.compareTo(BigNum::ONE) == 0) || ++j == s)
                return false;
            x = x.modPow(BigNum::TWO, *this);
        }
    }
    return true;
}

void swap(BigNum& first, BigNum& second)
{
    swap(first._bits, second._bits);
    swap(first._size, second._size);
    swap(first._num, second._num);
    swap(first._signum, second._signum);
}

BigNum BigNum::getPrime(unsigned int bitSize) {
    BigNum result(bitSize);
    result.setRandom();
    if (result.isEven()) {
        result = result - BigNum::ONE;
    }
    int r = 100 * result._bits;
    while (r > 0 && !(result.isPrime())) {
        result = result + BigNum::TWO;
        --r;
    }
    return result;
}

shared_ptr<BigNum> BigNum::extend(unsigned int bits) const {
    if (bits < _bits) {
        return nullptr;
    }

    shared_ptr<BigNum> result = make_shared<BigNum>(bits);
    for (size_t i = 0; i < _size; ++i) {
        result->_num[i] = _num[i];
    }
    return result;
}


void BigNum::setCorrectBitSize() {
    if (_bits != 0) {
        if (_num[_size - 1] == 0) {
            int i = _size - 2;
            while (i >= 0 && _num[i] == 0) {
                --i;
            }
            if (i < 0) {
                _size = 1;
                _bits = 1;
                delete[] _num;
                _num = new unsigned long[_size];
                _num[0] = 0;
                _signum = 0;
                return;
            }
            unsigned long *num = _num;
            _size = i + 1;
            _num = new unsigned long[_size];

            while (i >= 0) {
                _num[i] = num[i];
                --i;
            }
            delete[] num;
        }

        _bits = (_size - 1) * _baseBits + ceil(log2(_num[_size - 1]));
        bitset<32> leadCoefficient(_num[_size - 1]);
        if (leadCoefficient.count() == 1) {
            ++_bits;
        }
    }
}


shared_ptr<BigNum> BigNum::splitToomCook3(const int&start, const int&length) const {
    int sliceSize = length;
    if (start == 0 && length >= this->_size) {
        return make_shared<BigNum>(*this);
    }
    else if (start >= this->_size) {
        return make_shared<BigNum>(BigNum::ZERO);
    }
    else if (this->_size - start < length) {
        sliceSize = this->_size - start;
    }
    shared_ptr<BigNum> result = make_shared<BigNum>(sliceSize * this->_baseBits);
    for (int i = 0; i < sliceSize; ++i) {
        result->_num[i] = this->_num[i + start];
    }
    result->_signum = POSITIVE;
    result->setCorrectBitSize();
    return result;
}


// http://en.wikipedia.org/wiki/Toom%E2%80%93Cook_multiplication
shared_ptr<BigNum> BigNum::toomCook3(const BigNum& bigNum) const {
    int largest = max(this->_size, bigNum._size);
    double largest3 = 0;
    double fractpart = modf(((double) largest) / 3, &largest3);
    int minSize = largest3;
    int maxSize = largest - 2 * minSize;
    if (fractpart > 0.5) {
        maxSize = (largest - minSize) / 2;
    }

//    int size1 = largest - minSize - maxSize;
    int start2 = largest - 2 * maxSize;
    shared_ptr<BigNum> m0 = this->splitToomCook3(0, maxSize);
    shared_ptr<BigNum> m1 = this->splitToomCook3(maxSize, maxSize);
    shared_ptr<BigNum> m2 = this->splitToomCook3(maxSize + maxSize, minSize);

    shared_ptr<BigNum> n0 = bigNum.splitToomCook3(0, maxSize);
    shared_ptr<BigNum> n1 = bigNum.splitToomCook3(maxSize, maxSize);
    shared_ptr<BigNum> n2 = bigNum.splitToomCook3(maxSize + maxSize, minSize);

    BigNum p = *m0 + *m2;
    BigNum p1 = p + *m1;
    BigNum pNeg1 = p - *m1;
    BigNum pNeg2 = ((pNeg1 + *m2) << 1) - *m0;

    BigNum q = *n0 + *n2;
    BigNum q1 = q + *n1;
    BigNum qNeg1 = q - *n1;
    BigNum qNeg2 = ((qNeg1 + *n2) << 1) - *n0;

    BigNum r0 = (*m0) * (*n0);
    BigNum r1 = p1 * q1;
    BigNum rNeg1 = pNeg1 * qNeg1;
    BigNum rNeg2 = pNeg2 * qNeg2;
    BigNum rInf = (*m2) * (*n2);

    BigNum a3 = (rNeg2 - r1) / 3;
    BigNum a1 = (r1 - rNeg1) >> 1;
    BigNum a2 = rNeg1 - r0;
    a3 = ((a2 - a3) >> 1) + (rInf << 1);
    a2 = a2 + a1 - rInf;
    a1 = a1 - a3;

    int base = maxSize * this->_baseBits;
    BigNum result = (((((((rInf << base) + a3) << base) + a2) << base) + a1) << base) + r0;
    result.setCorrectBitSize();
    return make_shared<BigNum>(result);
}

void BigNum::setRandom() {
    ifstream f("/dev/urandom", ifstream::in);
    if (!f) {
        cout << "/dev/urandom does not exist" << endl;
        f.close();
    }

    int randomLen = 0;

    int bytes = _baseBits / 8;

    for (int i = 0; i < _size; ++i) {
        while (randomLen < bytes) {
            f.read(((char *) &_num[i]) + randomLen, bytes - randomLen);
            randomLen += f.gcount();
            if (f.fail()) {
                cout << "Some problems with file" << endl;
                f.close();
            }
        }
        randomLen = 0;
    }
    f.close();
    this->setCorrectBitSize();
    if (this->_signum == 0)
        this->_signum = POSITIVE;
}

void BigNum::print() const {
    cout << _signum << endl;
    if (_signum != 0)
        for (int i = 0; i < _size; ++i) {
            cout << _num[i] << endl;
        }
}

void BigNum::print(ofstream &f) const {
    f << _signum << endl;
    if (_signum != 0)
        for (int i = 0; i < _size; ++i) {
            f << _num[i] << endl;
        }
}

void BigNum::read(ifstream& f) {
    unsigned long temp = 0;
    f >> temp;
    _signum = temp;
    for (int i = 0; i < _size; ++i) {
        f >> temp;
        _num[i] = temp;
    }
    this->setCorrectBitSize();
}
