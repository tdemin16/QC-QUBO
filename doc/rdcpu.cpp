#include <cstring>
#include <iostream>
#include <random>

// static bool RDRAND(void) { return CPU_Rep.f_1_ECX_[30]; }
// static bool RDSEED(void) { return CPU_Rep.f_7_EBX_[18]; }
// EAX=7, ECX=0: Extended Features
// <x86intrin.h> OR <immintrin.h>
using namespace std;
// compile with gcc -mrdrnd -mrdseed
#include "cpuid.h"
int main(void) {
    //  const unsigned int flag_RDRAND = (1u << 30);
    //  const unsigned int flag_RDSEED = (1u << 18);
    unsigned int rnd;
    unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
    unsigned int ext = 0x8000000;

    int a = __get_cpuid_max(ext, nullptr);
    cout << hex << "a= " << a << " ext: " << ext << endl;
    int b = __get_cpuid(0x7, &eax, &ebx, &ecx, &edx);
    cout << "leaf 7 = " << b << " ecx: " << ecx << endl;

    __cpuid(0x01, eax, ebx, ecx, edx);
    cout << "ecx: " << hex << ecx << endl;
    if (ecx & 1) cout << "SSE3 supported" << endl;
    if (ecx & bit_RDRND) {
        cout << "Using RDRAND32 CPU intrinsic (Ivy Bridge)" << endl;
        for (int i = 0; i < 10; i++) {
            _rdrand32_step(&rnd);
            cout << dec << rnd << endl;
        }
    }
    __cpuid(0x07, eax, ebx, ecx, edx);
    cout << "ebx: " << hex << ebx << endl;
    //  cout << "Testing RDSEED32 " << bit_RDSEED << " CPU intrinsic (Broadwell, AMD Zen)" << endl;
    if (ebx & bit_RDSEED) {
        cout << "Using RDSEED32 " << bit_RDSEED << " CPU intrinsic (Broadwell, AMD Zen)" << endl;
        for (int i = 0; i < 10; i++) {
            _rdseed32_step(&rnd);
            cout << rnd << endl;
        }
    }
}
