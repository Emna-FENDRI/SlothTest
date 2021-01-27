
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <gmp.h>
#include <string.h>
#include <math.h>  


//#include <statistics.h>

//100000
#define ITERATIONS 1
#define TOT_INT 32
#define TOT_INT_str "32"
#define MASK_str "100000000000000000000000" // 93 iéme bit HEXA
#define SAMPLE_SIZE 100000
#define s 41
#define PATH "/Users/emnafendri/Desktop/sloth_core/Data/new_100000_32_s41_xor.csv"
// Remember xor


void xor_mod(mpz_t result, const mpz_t input1, const mpz_t flip, const mpz_t mod) {
    mpz_xor(result,input1,flip);
    while (mpz_cmp(result, mod) >= 0 || mpz_cmp_ui(result, 0) == 0) {
        mpz_xor(result,result,flip);
    }
}

// the sqrt permutation as specified in the paper (returns a sqrt of either input or -input)
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


// computes witness = the sloth witness, for the given seed, number of iterations and prime p
void sloth_core(mpz_t witness, const mpz_t seed, int iterations, const mpz_t p) {
    mpz_t a, ones, e;
    //a est initialisé à seed
    mpz_init_set(a, seed);
    
    mpz_init_set_ui(ones, 1); //init ones = 1 , void mpz_init_set_ui (MP_INT *dest_integer, unsigned long int src_ulong)
    mpz_mul_2exp(ones, ones, mpz_sizeinbase(p,2) >> 1); // flip half the bits (least significant), sizeinbase retourne le nombre de digit dans p (base 2) et >>1 pour diviser par 2 ce nombre (c un size_t)
    mpz_sub_ui(ones, ones, 1);
    
    // compute the exponent for sqrt extraction
    
    mpz_init_set(e, p);
    mpz_add_ui(e, e, 1);
    mpz_tdiv_q_ui(e, e, 4);
    
    for (int i = 0; i < iterations; ++i) {
        xor_mod(a,a,ones,p);  //HERE
        sqrt_permutation(a, a, p, e);
    }
    
    //witness devient a
    mpz_set(witness, a);
    
    mpz_clear(a);
    mpz_clear(ones);
}


void generate_prime(mpz_t prime){

    gmp_randstate_t state; 

    //Initialize state with a default algorithm. This will be a compromise between speed and randomness,
    // and is recommended for applications with no special requirements.
    gmp_randinit_default(state); //  Initi  laze it with default function 

    mp_bitcnt_t n; //length prime bits
    n = 2048;
    
    while(!mpz_probab_prime_p(prime, 25)){
        //mpz_add_ui(prime,prime,1); 
        //void mpz_urandomb (mpz_t rop, gmp_randstate_t state, unsigned long int n)
        //Generate a uniformly distributed random integer in the range 0 to 2^n-1, inclusive.
        mpz_urandomb(prime, state, n);
    } 

    gmp_randclear(state);

}




// VERSION1 O(n)
int get_interval_1(const mpz_t output_sloth, const mpz_t prime ){
    
    mpz_t length; // length of one chunk
    mpz_t n; // number of chunks
    mpz_init(length);
    mpz_init_set_str(n,TOT_INT_str,10);//100

   
   mpz_cdiv_q(length , prime , n); // prime /128
   
    // COMPARE
    //Compare op1 and op2. Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
    mpz_t zero;
    mpz_t one;
    
    mpz_t lower;
    mpz_t higher;
    mpz_init(zero);
    mpz_init_set_str(one, "1", 10);
    
    mpz_init_set (lower, one); //lower = 1
    mpz_init_set (higher, length);
    
    int index = 0;
    while( index < TOT_INT ){

        int check1 = mpz_cmp(output_sloth,lower);
        int check2 = mpz_cmp(higher, output_sloth);

        if((check1>0 )&& (check2>0)){
            return index ;
      
        }
        //Function: void mpz_add (mpz_t rop, const mpz_t op1, const mpz_t op2)
        mpz_add(lower, lower, length); //lower += length;
        mpz_add(higher, lower, length); // higher += length;

       
        ++index;
    }

    mpz_clear(zero);
    mpz_clear(one);
    
    mpz_clear(length);
    mpz_clear(n);

    mpz_clear(lower);
    mpz_clear(higher);
    return -1;
}


void print_vect(int** tab , int size){
    if(tab != NULL){
        for(int i =0; i< size; ++i){
            for(int j =0; j< size; ++j){
            printf("Interval (%d,%d): %20d\n", i, j ,tab[i][j]);
            }
        }
    }
}

double Stat( int** obs, const int* X_total , const int* Y_total){
    double T = 0;
    if((obs != NULL) && (X_total != NULL) && (Y_total != NULL)){

        for(int i =0; i<TOT_INT ; ++i){
            
            if(obs[i] != NULL){
                for(int j=0; j<TOT_INT; ++j){
                    double E_val = ((double)X_total[i]*(double)Y_total[j])/(double)SAMPLE_SIZE;
                    if(E_val == 0){ 
                        fprintf(stderr, "Error: DIV/0 \n");
                        return -1;}

                    T += pow((obs[i][j]-E_val), 2.0)/E_val;
                }
            }
            //return -1; // not sure
        }
        return T;
    }
    return -1;
}

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
        fprintf(stderr,"Cannot open file");
    }
    fclose(file);

}

void get_input(){

    int** vect =(int**)calloc(TOT_INT, sizeof(int*));
    if (vect != NULL){
        for(int i =0 ; i<TOT_INT; ++i ){
            vect[i] = (int*)calloc(TOT_INT, sizeof(int));
            if(vect[i]==NULL){
                fprintf(stderr, "Error: Heap allocation\n");
            }
        }

        int* X_total = (int*)calloc(TOT_INT, sizeof(int));
        int* Y_total = (int*)calloc(TOT_INT, sizeof(int));
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



        //f6ec558f285278de66f9e6717cc79275189328f350e3cf8d4b3235dda6714badc5f8c888723c3f0554fd1ab46ffafbaa7a9305f53c0a2788d71bc07667edd151
        //Define fixed and known prime (512bits)
        // HEX: 89991233b87081c8a95ffbec76352d6c7631ac54e7946696e6f4fbab61d1cee14fa6cb59c3074b079d1c6fdf1d1cd4fe44306dcb5386bd8aecc95784b83bbb73
        // b10 :7206588556671961696914174949518405627123847171625399748521614959002826444849593981836508543000436246142568273303062039539289673298832254224814417525848947
        // mpz_init_set_str(prime,"89991233b87081c8a95ffbec76352d6c7631ac54e7946696e6f4fbab61d1cee14fa6cb59c3074b079d1c6fdf1d1cd4fe44306dcb5386bd8aecc95784b83bbb73", 16);
       
        generate_prime(prime);

        /*printf ("Prime:");
        mpz_out_str(stdout,10,prime);
        printf("\n");*/

        unsigned long se = s; 
        gmp_randinit_default(state);
        gmp_randseed_ui(state,se);

    for(int i=0; i<SAMPLE_SIZE; ++i) {
            printf("loop %d\n",i);
            mpz_urandomm (seed_X,state,prime); // generate random input in the range 0 and prime-1
            mpz_xor (seed_Y, seed_X, mask); // create seed_Y = X(xor)MASK

            // print Seed
            /*printf("Seed X:");
            mpz_out_str(stdout,16,seed_X);
            printf("\n");
            printf("Seed Y:");
            mpz_out_str(stdout,16,seed_Y);
            printf("\n");*/

            sloth_core(output_sloth_X,seed_X, ITERATIONS, prime );
            sloth_core(output_sloth_Y,seed_Y, ITERATIONS, prime );

            // print Output
            /*printf("output X:");
            mpz_out_str(stdout,16,output_sloth_X);
            printf("\n");
            printf("output Y:");
            mpz_out_str(stdout,16,output_sloth_Y);
            printf("\n");*/
            
            int r_X = mpz_mod_ui(r_X_mpz,output_sloth_X, TOT_INT); 
            int r_Y = mpz_mod_ui(r_Y_mpz,output_sloth_Y, TOT_INT);

            ++vect[r_X][r_Y];
            ++ X_total[r_X];
            ++ Y_total[r_Y];

            
            
    }
        print_vect(vect , TOT_INT);
        
        double stat = Stat(vect, X_total, Y_total);
        printf("Chi-Square Statistic: %f \n", stat);
        //write in file
        write_csv(vect, TOT_INT, stat);

       

        for(int i =0 ; i<TOT_INT; ++i ){
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

    get_input();
    return 0;
    

    //alglib::chisquaredistribution(16, stat);
    //double E_val = ((double)100032*(double)100053)/(double)SAMPLE_SIZE;
    //printf("E_VAL : %f \n", E_val);
    //Prime p
    // ou en b10 nv : 12932416320295881761671611278642326573505016881349297659022122875335761010885136553022759118034562056831395799829560255106185681381489452149779058467393873
    //201434053360839760574708769124663506817179517656154357328732411229106995022301843102615358007453227521115601504656579158239306825468604226460172665412541
    // 256 : 87891909009787391741064284376429281279755419278687755710227881169949927172819
    // 512 : 420143405336083976057470876912466350681717951765615435mak7328732411229106995022301843102615358007453227521115601504656579158239306825468604226460172665412541
    //input : 4201434053360839760574708769124604226460172665412500
   // p en hex : 3d896d2159fef9200f7071deee377a2b60579aa8b260c2731d37cb71ec96e6109dde575cdeed74be4f876f7044515839826cff0e664ac070a45339498eedbbd
    //mpz_init_set_str(p,"3d896d2159fef9200f7071deee377a2b60579aa8b260c2731d37cb71ec96e6109dde575cdeed74be4f876f7044515839826cff0e664ac070a45339498eedbbd", 16);
    //char input[1024];


    
}

