#include <iostream>
#include <cstdlib>
#include <ctime>
#include "BigNum.h"
#include "RSA.h"

using namespace std;

const long N = 3;
const int bits = 512;

int main() {
    unsigned long sum = 0;
    double sumTime = 0;
    clock_t start, end;
    time_t startTime, endTime;
    unsigned long long max = 0, min = ULLONG_MAX;
    double maxTime = 0, minTime = ULLONG_MAX;
    clock_t different;
    double differentTime;

    for (long i = 0; i < N; ++i) {
        cout << i << endl;

        start = clock();
        startTime = time(NULL);

        RSA rsa(bits);
        string message = "Hello World!";
        int cyphersSize = 0;
        BigNum* cyphers = rsa.encrypt(message, cyphersSize);
        string message1 = rsa.decrypt(cyphers, cyphersSize);
        cout << endl << message1 << endl;
        delete[] cyphers;


        end = clock();
        endTime = time(NULL);
        differentTime = difftime(endTime, startTime);
        sumTime += differentTime;

        different = end - start;

        sum += different;
        if (min > different) {
            min = different;
        }
        if (max < different) {
            max = different;
        }
        if (minTime > differentTime) {
            minTime = differentTime;
        }
        if (maxTime < differentTime) {
            maxTime = differentTime;
        }
    }
    cout << "Ticks" << endl;
    cout << "max =     " << max << endl;
    cout << "min =     " << min << endl;
    cout << "sum =     " << sum << endl;
    cout << "sum / N = " << sum / N << endl;
    cout << "Time (s)" << endl;
    cout << "max =     " << maxTime << endl;
    cout << "min =     " << minTime << endl;
    cout << "sum =     " << sumTime << endl;
    cout << "sum / N = " << sumTime / N << endl;
    return 0;
}

