#include <iostream>
#include <omp.h>
#include<time.h>
using namespace std;
double p = 73;
double q = 31;
double e = 2;
int threadNum;
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
long long key(long long e , long long iphi) {
    for (long long x = 1; x < iphi; x++) {
        if ((e * x) % iphi == 1) {
            return x;
        }
    }
    return -1; // Inverse doesn't exist
}
int modularExponentiation(long long base, long long exponent, long long modulus) {
    long long prod = 1;
        for (long long i = 0; i < exponent; i++) {
            
            prod = (prod * base) % modulus;
        } 
    return prod;
}
bool areArraysEqual(long long arr1[], long long arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false;  // Arrays are not equal
        }
    }
    return true;  // Arrays are equal
}
#define size  1000

int main(int argc, char* argv[])
{
    omp_set_num_threads(6);
    time_t start_t, end_t;
    long long* Message = new long long[size]; 


    srand(time(NULL));
    start_t = clock();
    double start = omp_get_wtime();

    double in = n_find(p, q);
    double iphi = phi_find(p, q);
    while (e < iphi) {
        // e must be co-prime to phi and
        // smaller than phi.
        if (gcd(e, iphi) == 1)
            break;
        else
            e++;
    }
    double ik = key(e, iphi);
 
    int max = 4, min = 2;
    int range = max - min + 1;
    for (int i = 0; i < size;++i) {
        // Initialize random number generator
        Message[i] = rand()% range + min;
    }
    long long* encrypted = (long long*)malloc(size * sizeof(long long));
    long long* decrypted = (long long*)malloc(size * sizeof(long long));
    

#pragma omp parallel for
    
        for (int i = 0; i < size; i++) {
             threadNum = omp_get_thread_num();
            // RSA encryption on the block
            encrypted[i] = modularExponentiation(Message[i], e, in);
            // RSA decryption on the block
            decrypted[i] = modularExponentiation(encrypted[i], ik, in);
           // printf("Thread %d:  %lld\n", threadNum, Message[i]);
        //   printf("Thread %d: Encrypted Message: %lld\n", threadNum, encrypted[i]);
        //  printf("Thread %d: Decrypted Message: %lld\n", threadNum, decrypted[i]);
          
        }
       // printf("Thread %d:  ", threadNum);
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
        double time_taken_t = (double)(end - start)/ (double)(CLOCKS_PER_SEC);
        printf("time taken is %f \n", time_taken);
        printf("time from time.h taken is %f \n", time_taken_t);
        free(encrypted);
        free(decrypted);
        delete[] Message;
           return 0;

    }
