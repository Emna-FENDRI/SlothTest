
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <gmp.h>
#include <string.h>
#include <math.h>  


#define ITERATIONS 1
#define MOD 4
#define MOD_str "4"
#define MASK_str "100000000000000000000000" // 93 in HEXA
#define SAMPLE_SIZE 20000
#define s 11
#define PATH "/Users/emnafendri/Desktop/sloth_core/Data/m4_20000_s11_xor.csv"



void xor_mod(mpz_t result, const mpz_t input1, const mpz_t flip, const mpz_t mod) {
    mpz_xor(result,input1,flip);
    while (mpz_cmp(result, mod) >= 0 || mpz_cmp_ui(result, 0) == 0) {
        mpz_xor(result,result,flip);
    }
}

// The sqrt permutation as specified in the paper (returns a sqrt of either input or -input)
void sqrt_permutation(mpz_t result, const mpz_t input, const mpz_t p, const mpz_t e) {
	mpz_t tmp;
	mpz_init(tmp);
	if (mpz_jacobi(input, p) == 1) {
		mpz_powm(tmp, input, e, p);
		if (mpz_even_p(tmp)) mpz_set(result, tmp);
		else mpz_sub(result, p, tmp);
	}
	else {
		mpz_sub(tmp, p, input);
        mpz_powm(tmp, tmp, e, p);
		if (mpz_odd_p(tmp)) mpz_set(result, tmp);
		else mpz_sub(result, p, tmp);
	}

	mpz_clear(tmp);
}


void sloth_core(mpz_t witness, const mpz_t seed, int iterations, const mpz_t p) {
    mpz_t a, ones, e;
    mpz_init_set(a, seed);
    
    mpz_init_set_ui(ones, 1); 
    mpz_mul_2exp(ones, ones, mpz_sizeinbase(p,2) >> 1);
    mpz_sub_ui(ones, ones, 1);
    
    // compute the exponent for sqrt extraction
    
    mpz_init_set(e, p);
    mpz_add_ui(e, e, 1);
    mpz_tdiv_q_ui(e, e, 4);
    
    for (int i = 0; i < iterations; ++i) {
        xor_mod(a,a,ones,p);
        sqrt_permutation(a, a, p, e);
    }
    
    mpz_set(witness, a);
    
    mpz_clear(a);
    mpz_clear(ones);
}

/**
 * @brief 
 * Generate a prime in the range 0 to 2^2048 
 * @param prime prime generated
 */
void generate_prime(mpz_t prime){
    gmp_randstate_t state; 
    //Initialize state with a default algorithm.
    gmp_randinit_default(state);

    mp_bitcnt_t n; //length prime bits
    n = 2048;
    
    while(!mpz_probab_prime_p(prime, 25)){
        //Generate a uniformly distributed random integer in the range 0 to 2^n-1, inclusive.
        mpz_urandomb(prime, state, n);
    } 
    gmp_randclear(state);
}

/**
 * @brief prints the given 2D array of the given size 
 * 
 * @param tab 
 * @param size 
 */
void print_vect(int** tab , int size){
    if(tab != NULL){
        for(int i =0; i< size; ++i){
            for(int j =0; j< size; ++j){
            printf("Interval (%d,%d): %20d\n", i, j ,tab[i][j]);
            }
        }
    }
}

/**
 * @brief Computes the Chi-squared test statistic
 * 
 * @param obs Contingency table, stores all observations
 * @param X_total row totals
 * @param Y_total column totals
 * @return double  
 */
double Stat( int** obs, const int* X_total , const int* Y_total){
    double T = 0;
    if((obs != NULL) && (X_total != NULL) && (Y_total != NULL)){
        for(int i =0; i<MOD ; ++i){
            if(obs[i] != NULL){
                for(int j=0; j<MOD; ++j){
                    double E_val = ((double)X_total[i]*(double)Y_total[j])/(double)SAMPLE_SIZE;
                    if(E_val == 0){ 
                        fprintf(stderr, "Error: DIV/0 \n");
                        return -1;}
                    T += pow((obs[i][j]-E_val), 2.0)/E_val;
                }
            }
        }
        return T;
    }
    return -1;
}

/**
 * @brief writes in a csv file the value of the Test Statistic and 
 * all the observations as : number of observations per category [i][j]
 * 
 * @param obs Contingency table, stores all observations
 * @param size of the given table
 * @param T Test Statistic
 */
void write_csv(int** obs , int size , double T ){
    FILE *file;
    file = fopen(PATH, "w");
    if (file != NULL){
        fprintf(file, "Modulo, Observations\n");
        for(int i=0; i<size; ++i){
            if(obs[i] != NULL){
                for(int j=0; j<size ; ++j){
                    fprintf(file, "[%d][%d],%d\n",i,j,obs[i][j]);
                }

            }
        }  
        fprintf(file, "Test Statistic = %lf \n" , T);
    }   
    else{
        fprintf(stderr,"Impossible to open file");
    }
    fclose(file);

}

/**
 * @brief test computation
 * 
 */
void test(){
    
    int** vect =(int**)calloc(MOD, sizeof(int*));
    if (vect != NULL){
        for(int i =0 ; i<MOD; ++i ){
            vect[i] = (int*)calloc(MOD, sizeof(int));
            if(vect[i]==NULL){
                fprintf(stderr, "Error: Heap allocation\n");
            }
        }

        int* X_total = (int*)calloc(MOD, sizeof(int));
        int* Y_total = (int*)calloc(MOD, sizeof(int));
        if((X_total == NULL) || (Y_total == NULL) ){
            fprintf(stderr, "Error: Heap allocation\n");
        }

        mpz_t seed_X, seed_Y,mask;
        mpz_t output_sloth_X, output_sloth_Y,prime;
        gmp_randstate_t state;
        mpz_t r_X_mpz, r_Y_mpz;

        
        mpz_init(seed_X);
        mpz_init(seed_Y);
        mpz_init_set_str(mask, MASK_str, 16);
        mpz_init(output_sloth_X);
        mpz_init(output_sloth_Y);
        mpz_init(prime);
        mpz_init(r_X_mpz);
        mpz_init(r_Y_mpz);

        generate_prime(prime);

        unsigned long se = s; 
        gmp_randinit_default(state);
        gmp_randseed_ui(state,se);

    for(int i=0; i<SAMPLE_SIZE; ++i) {
            printf("loop %d\n",i);
            mpz_urandomm (seed_X,state,prime); // generate random input in the range 0 and prime-1
            mpz_xor (seed_Y, seed_X, mask); // create seed_Y = X(xor)MASK

            sloth_core(output_sloth_X,seed_X, ITERATIONS, prime );
            sloth_core(output_sloth_Y,seed_Y, ITERATIONS, prime );

            int r_X = mpz_mod_ui(r_X_mpz,output_sloth_X, MOD); 
            int r_Y = mpz_mod_ui(r_Y_mpz,output_sloth_Y, MOD);

            ++vect[r_X][r_Y];
            ++ X_total[r_X];
            ++ Y_total[r_Y];
    
    }
        print_vect(vect , MOD);
        double stat = Stat(vect, X_total, Y_total);
        printf("Chi-Square Statistic: %f \n", stat);

        //write in file
        write_csv(vect, MOD, stat);

        //free allocated memory
        for(int i =0 ; i<MOD; ++i ){
            free(vect[i]);
        }
        free(vect);

        gmp_randclear(state);
        mpz_clear(seed_X);
        mpz_clear(seed_Y);
        mpz_clear(mask);
        mpz_clear(output_sloth_X);
        mpz_clear(output_sloth_Y);
        mpz_clear(prime);   
        mpz_clear(r_X_mpz);
        mpz_clear(r_Y_mpz);

    }
}

int main() {

    test();
    return 0;
    
}

