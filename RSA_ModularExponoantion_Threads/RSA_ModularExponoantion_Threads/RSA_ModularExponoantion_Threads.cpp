#include <iostream>
#include <omp.h>
#include<time.h>
using namespace std;
double p = 73;
double q = 31;
double e = 2;

int gcd(int a, int h)// Returns gcd of a and b
{
    int temp;
    while (1) {
        temp = a % h;
        if (temp == 0)
            return h;
        a = h;
        h = temp;
    }
}
double n_find(double p, double q) {
    double n = p * q;
    return n;
}
double phi_find(double p, double q) {
    double phi = (p - 1) * (q - 1);
    return phi;
}
long long key(long long e, long long iphi) {
    for (long long x = 1; x < iphi; x++) {
        if ((e * x) % iphi == 1) {
            return x;
        }
    }
    return -1; // Inverse doesn't exist
}

bool areArraysEqual(long long arr1[], long long arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false;  // Arrays are not equal
        }
    }
    return true;  // Arrays are equal
}
#define size 100000000
int main() {
    time_t start_t, end_t;
    omp_set_num_threads(5);
    long long* Message = new long long[size];
    long long* encrypted = (long long*)malloc(size * sizeof(long long));
    long long* decrypted = (long long*)malloc(size * sizeof(long long));
    srand(time(NULL));
    start_t = clock();
    double start = omp_get_wtime();
    long long in = n_find(p, q);
    long long iphi = phi_find(p, q);
   
    while (e < iphi) {
        // e must be co-prime to phi and
        // smaller than phi.
        if (gcd(e, iphi) == 1)
            break;
        else
            e++;
    }
    long long ik = key(e, iphi);
  // printf("e is %f\n", e);
   // printf("k is %lld\n", ik);
   // printf("n is %lld\n", in);
    int max = 4, min = 2;
    int range = max - min + 1;
    for (int i = 0; i < size;++i) {
        // Initialize random number generator
        Message[i] =  rand() % range + min;
        encrypted[i] = 1;
        decrypted[i] = 1;
        //printf("%lld\n", Message[i]);
    }
   
   
#pragma omp parallel 
    {
        int chunk_e = e / omp_get_num_threads();
        int chunkRemainder_e = int(e) % omp_get_num_threads();
        int start_e = chunk_e * omp_get_thread_num();
        int end_e = start_e + chunk_e;

        if (omp_get_thread_num() == omp_get_num_threads() - 1) {
            end_e = end_e + chunkRemainder_e;
        }

        long long localResult = 1;
        for (int j = 0;j < size;j++) {
            for (int i = start_e; i < end_e; i++) {
                localResult = (localResult * Message[j]) % in;
            }
#pragma omp critical
            encrypted[j] = (encrypted[j] * localResult) % in;
            localResult = 1;
        }
    }
#pragma omp parallel 
    {
        int chunk_k = ik / omp_get_num_threads();
        int chunkRemainder_k = int(ik) % omp_get_num_threads();
        int start_k = chunk_k * omp_get_thread_num();
        int end_k = start_k + chunk_k;

        if (omp_get_thread_num() == omp_get_num_threads() - 1) {
            end_k = end_k + chunkRemainder_k;
        }

        long long localResult_k = 1;
        for (int j = 0;j < size;j++) {
            for (int i = start_k; i < end_k; ++i) {
                localResult_k = (localResult_k * encrypted[j]) % in;
            }
#pragma omp critical
            decrypted[j] = (decrypted[j] * localResult_k)% in;
            localResult_k = 1;
        }
     
    }
    //for (int i = 0; i < size;++i) {
        // Initialize random number generator
        
        //printf("Message is %lld\n", Message[i]);
       // printf("encryptyion is %lld\n", encrypted[i]);
       // printf("decryption is %lld\n", decrypted[i]);
   // }
    bool equal = areArraysEqual(Message, decrypted, size);
    if (equal) {
        printf("Arrays are equal.\n");
    }
    else {
        printf("Arrays are not equal.\n");
    }
    double end = omp_get_wtime();
    end_t = clock();
    double time_taken = (double)(end - start); /// (double)(CLOCKS_PER_SEC);
    double time_taken_t = (double)(end - start) / (double)(CLOCKS_PER_SEC);
    printf("time taken is %f \n", time_taken);
    printf("time from time.h taken is %f \n", time_taken_t);
    free(encrypted);
    free(decrypted);
    delete[] Message;
    return 0;
}