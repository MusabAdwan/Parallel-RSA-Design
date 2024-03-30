#include <iostream>
#include <mpi.h>
#include<time.h>
using namespace std;
double p = 10259;
double q = 10253;
double e = 5;
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
long long  key(long long e, long long phi) {
    for (long long x = 1; x < phi; x++) {
        if ((e * x) % phi == 1) {
            return x;
        }
    }
    return -1; // Inverse doesn't exist
}
long long modularExponentiation(long long base, long long exponent, long long modulus) {
    if (exponent == 0)
        return 1;

    long long prod = 1;

    while (exponent > 0) {
        if (exponent % 2 == 1)
            prod = (prod * base) % modulus;
        exponent = exponent / 2;
        base = (base * base) % modulus;
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
#define size  100000000
int main(int argc, char** argv) {

    long long* Message = new long long[size];
    long long* encrypted = (long long*)malloc(size * sizeof(long long));;
    long long* dycrypted = (long long*)malloc(size * sizeof(long long));;
    long long* finalarray = (long long*)malloc(size * sizeof(long long));;
    long long ik, in, iphi;
    srand(time(NULL));
    time_t start_t, end_t;
  
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    start_t = clock();

    if (myrank == 0) {
         in = n_find(p, q);
         iphi = phi_find(p, q);
        while (e < iphi) {
            // e must be co-prime to phi and
            // smaller than phi.
            if (gcd(e, iphi) == 1)
                break;
            else
                e++;
        }

         ik = key(e, iphi);
        for (int i = 0; i < size;++i) {
            // Initialize random number generator
           // int range = 4 - 2 + 1;
            Message[i] = rand();// % range + 2;
        }
    }
    
    MPI_Bcast(&Message, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&in, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iphi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int chunkSize = size / npes;
    int modsize = size % npes;
    int startpoint = myrank * chunkSize;
    int endpoint;
    if (myrank == npes - 1) {
        endpoint = size;
    }
    else {
        endpoint = startpoint + chunkSize;
    }
    for (int i = startpoint;i < endpoint;++i) {

        encrypted[i] = modularExponentiation(Message[i], e, in);
        dycrypted[i] = modularExponentiation(encrypted[i], ik, in);
        //  printf("Rank is %d i %d  message. %lld\n", myrank, i, Message[i]);
        // printf("Rank is %d Arrays %d encrypted. %lld\n", myrank,i, encrypted[i]);
      //  printf("Rank is %d Arrays%d decrypted .%lld\n", myrank, i, dycrypted[i]);
    }
    end_t = clock();
 /*   if (myrank == npes - 1) {
        for (int i = 0;i < size;++i) {
            finalarray[i]=0;
        }
    MPI_Gather(dycrypted, chunkSize, MPI_LONG_LONG, finalarray, chunkSize, MPI_LONG_LONG, npes - 1, MPI_COMM_WORLD);
    if (myrank == npes - 1) {
        for (int i = 0;i < size;++i) {
            printf("Rank is %d Arrays  %lld\n", myrank, finalarray[i]);
        }
    }
        //  printf("Rank is %d finalarray  %lld\n", myrank, sizeof(finalarray));
        //  printf("Rank is %d dycrypted  %lld\n", myrank, sizeof(dycrypted));
          bool equal = areArraysEqual(Message, finalarray, size);

          if (equal) {
              printf("Arrays are equal.\n");
          }
          else {
              printf("Arrays are not equal.\n");
          }
      }*/

    double time_taken = (double)(end_t - start_t) / (double)(CLOCKS_PER_SEC);
    printf("My Rank is %d and the time taken is %f \n",myrank, time_taken);
    delete[]Message;
    free(encrypted);
    free(dycrypted);
    free(finalarray);
    MPI_Finalize();
    return 0;

}

