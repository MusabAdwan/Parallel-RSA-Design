#include <iostream>
#include<time.h>
using namespace std;

// Function to calculate modular exponentiation (x^y mod m)

long long modularExponentiation(long long x, unsigned int y, int p) {
    int res = 1;     // Initialize result 

    x = x % p; // Update x if it is more than or 
    // equal to p

    if (x == 0) return 0; // In case x is divisible by p;

    while (y > 0)
    {
        // If y is odd, multiply x with result 
        if (y & 1)
            res = (res * x) % p;

        // y must be even now 
        y = y >> 1; // y = y/2 
        x = (x * x) % p;
    }
    return res;
}





// Function to find the modular multiplicative inverse using Extended Euclidean algorithm
long long key(long long a, long long m) {
    for (long long x = 1; x < m; x++) {
        if ((a * x) % m == 1) {
            return x;
        }
    }
    return -1; // Inverse doesn't exist
}bool areArraysEqual(long long arr1[], long long arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false;  // Arrays are not equal
        }
    }
    return true;  // Arrays are equal
}
// Function to calculate the greatest common divisor (GCD) using Euclidean algorithm
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

#define size 10
int main() {
    long long* Message = new long long[size];
    long long* encrypted = new long long[size];
    long long* decrypted = new long long[size];
    srand(time(NULL));
    time_t start_t, end_t;

    // Key generation (In practice, these values are kept secret)

    long long p = 11; // Prime number
    long long q = 13; // Prime number
    long long n = p * q; // Modulus
    long long phi = (p - 1) * (q - 1); // Euler's totient function

    // Choose public exponent (e) - commonly 65537
    long long e = 2;

    // Ensure e and phi are coprime
    while (e < phi) {
        // e must be co-prime to phi and
        // smaller than phi.
        if (gcd(e, phi) == 1)
            break;
        else
            e++;
    }
    // Calculate private exponent (d)
    long long d = key(e, phi);

    // Message to be encrypted (as string)
    int max = 104, min = 2;
    int range = max - min + 1;
    for (int i = 0;i < size;i++)
    {
        Message[i] = rand() % range + min;
    }
    start_t = clock();
    // Encrypt the message block by block (ASCII code)
    for (int i = 0; i < size; i++) {
        // RSA encryption on the block
        encrypted[i] = modularExponentiation(Message[i], e, n);
        // RSA decryption on the block (for illustration purposes)
        decrypted[i] = modularExponentiation(encrypted[i], d, n);
        // Output results for each block
        printf("Original :  %lld)\n", Message[i]);
        printf("Encrypted : %lld\n", encrypted[i]);
        printf("Decrypted :  %lld)\n", decrypted[i]);
    }
    end_t = clock();
    double time_taken = (double)(end_t - start_t) / (double)(CLOCKS_PER_SEC);
    printf("time taken is %f \n", time_taken);
    bool equal = areArraysEqual(Message, decrypted, size);

    if (equal) {
        printf("Arrays are equal.\n");
    }
    else {
        printf("Arrays are not equal.\n");
    }
    delete[]Message;
    delete[]encrypted;
    delete[]decrypted;
    return 0;
}