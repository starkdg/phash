#include "ph_fft.h"
#include <assert.h>


int main(int argc, char **argv){

    const int N = 4;
    double signal[4] = { 1, 4, 3, 2};

	complex double X[N];
    assert(fft(signal, N, X) == 0);

    assert(abs(creal(X[0]) - 10.0) < 0.0000001);
    assert(abs(cimag(X[0]) - 0.0)  < 0.0000001);

    assert(abs(creal(X[1]) + 2.0) < 0.0000001);
    assert(abs(cimag(X[1]) - 2.0) < 0.0000001);

    assert(abs(creal(X[2]) + 2.0) < 0.0000001);
    assert(abs(cimag(X[2]) - 0.0) < 0.0000001);

    assert(abs(creal(X[3]) + 2.0) < 0.0000001);
    assert(abs(cimag(X[3]) + 2.0) < 0.0000001);

    double signal2[4] = { 3, 5, 9, 2};
   
    assert(fft(signal2, N, X) == 0);
    
    assert(abs(creal(X[0]) - 19.0) < 0.0000001);
    assert(abs(cimag(X[0]) - 0.0)  < 0.0000001);

    assert(abs(creal(X[1]) + 6.0) < 0.0000001);
    assert(abs(cimag(X[1]) - 3.0) < 0.0000001);

    assert(abs(creal(X[2]) - 5.0) < 0.0000001);
    assert(abs(cimag(X[2]) - 0.0) < 0.0000001);

    assert(abs(creal(X[3]) + 6.0) < 0.0000001);
    assert(abs(cimag(X[3]) + 3.0) < 0.0000001);

    return 0;
}
