#include <iostream>
#include <mpi.h>
#include<time.h>
using namespace std;

double p = 73;
double q = 31;
long long e = 2;
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

long long  key(long long e, long long iphi) {
    for (long long x = 1; x < iphi; x++) {
        if ((e * x) % iphi == 1) {
            return x;
        }
    }
    return -1; // Inverse doesn't exist
}

int modularExponentiation(long long x, unsigned int y, int p) {
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

bool areArraysEqual(long long  arr1[], long long  arr2[], int size) {
    for (int i = 0; i < size; i++) {
        if (arr1[i] != arr2[i]) {
            return false;  // Arrays are not equal
        }
    }
    return true;  // Arrays are equal
}


# define size 100000000
int main(int argc, char** argv) {

    srand(time(NULL));
    long long* Message = new long long[size];
    long long* encrypted = (long long*)malloc(size * sizeof(long long));
    long long* decrypted = (long long*)malloc(size * sizeof(long long));
    long long* finalencryption = (long long*)malloc(size * sizeof(long long));
    long long* finaldecryption = (long long*)malloc(size * sizeof(long long));
    long long in, ik, iphi;
    int max = 4, min = 2;
    int range = max - min + 1;
    for (int i = 0; i < size;i++) {
        Message[i] = rand() % range + min;
        encrypted[i] = 1;
        decrypted[i] = 1;
        finalencryption[i] = 1;
        finaldecryption[i] = 1;
    }
    time_t start_t, end_t;
   
    int myrank, npes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    start_t = clock();
    if (myrank == 0)
    {
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
        
    }
   // MPI_Bcast(&Message, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&encrypted, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&decrypted, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&finalencryption, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&finaldecryption, size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&in, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iphi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printf("in is %lld\n", in);
    //printf("ik is %lld\n", ik);
    //printf("here is % lld\n", e);
    long long chunk_size_e = (e + npes - 1) / npes;
    long long start_e = myrank * chunk_size_e + 1;
    long long end_e = (myrank + 1) * chunk_size_e;

    long long chunk_size_k = (ik + npes - 1) / npes;
    long long start_ik = myrank * chunk_size_k + 1;
    long long end_ik = (myrank + 1) * chunk_size_k;

    if (end_e > e) {
        end_e = e; // Last process handles the remaining part
    }

    if (end_ik > ik) {
        end_ik = ik; // Last process handles the remaining part
    }


    // Each process performs its part of the modular exponentiation
    long long partialResult = 1;
    for (int j = 0; j < size; j++) {
        for (long long i = start_e; i <= end_e; ++i) {
            partialResult = (partialResult * Message[j]) % in;
            // printf("Process %d: Result: %d %lld\n", rank, j, partialResult);

        }
        encrypted[j] = partialResult;
        // printf(" my rank %d %d,partialResult is is %lld \n", myrank, j,encrypted[j]);
        partialResult = 1;
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(encrypted, finalencryption, size, MPI_LONG_LONG, MPI_PROD, MPI_COMM_WORLD);

    for (int j = 0; j < size; j++)
    {
        finalencryption[j] = (partialResult * finalencryption[j]) % in;
    }
    partialResult = 1;
    for (int j = 0; j < size; j++) {
        for (long long i = start_ik; i <= end_ik; ++i) {
            partialResult = (partialResult * finalencryption[j]) % in;
            // printf("Process %d: Result: %d %lld\n", myrank, j, partialResult);
        }
        decrypted[j] = partialResult;
        partialResult = 1;
        //  printf("Message is %d %lld \n",j,Message[j]);
        //  printf("my ranks %d %d is decryption is is %lld \n", myrank,j, decrypted[j]);
    }
    //  MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(decrypted, finaldecryption, size, MPI_LONG_LONG, MPI_PROD, MPI_COMM_WORLD);
    // printf("my ranks %d is decryption is is %lld \n", myrank,  finaldecryption[0]);
    for (int j = 0; j < size; j++)
    {
        finaldecryption[j] = (partialResult * finaldecryption[j]) % in;
        //   printf("my rank %d finaldecryption %d is %lld \n", myrank,j,finaldecryption[j]);
    }
    // printf("my ranks %d is decryption is is %lld \n", myrank, finaldecryption[0]);
    end_t = clock();

    if (myrank == 0) {
        /*  for (int i = 0;i < size; i++) {
              printf("Message is  %lld finaldecryption %lld is % \n", Message[i], finaldecryption[i]);
          }*/
        bool equal = areArraysEqual(Message, finaldecryption, size);

        if (equal) {
            printf("Arrays are equal.\n");
        }
        else {
            printf("Arrays are not equal.\n");
        }
    }
    double time_taken = (double)(end_t - start_t) / (double)(CLOCKS_PER_SEC);
    printf("time taken is %f \n", time_taken);
    delete[]Message;
    free(encrypted);
    free(decrypted);
    free(finalencryption);
    free(finaldecryption);
    // delete[]finalencryption;
    // delete[]finaldecryption;
    MPI_Finalize();

    return 0;
}