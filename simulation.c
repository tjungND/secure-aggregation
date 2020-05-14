#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

#define E_FIXED_POINT 32
#define MAX_USER_NUM 10
#define MAX_DIMENSION 512


void encode(mpz_t x, double y){
    
    mpz_set_d(x, y * pow(2, E_FIXED_POINT));
    
}

void decode(double *x, mpz_t y){
    
    mpf_set_default_prec(E_FIXED_POINT + 32);
    mpf_t temp;
    mpf_init(temp);
    mpf_set_z(temp, y);
    mpf_div_2exp(temp, temp, E_FIXED_POINT);
    *x = mpf_get_d(temp);
    mpf_clear(temp);
}

void simpleHash(mpz_t x, time_t y, mpz_t prime_base, mpz_t modulus){ // Our hash function does not need be cryptographic. All we need is a deterministic function that makes sure the base remains to be large (i.e., the x in H(t)^{p_nu}=x^{p_nu} has at least one large prime factor). This does not affect the security because the input of the hash is known to everyone.
    // This function turns the UNIX time to a large integer deterministically by computing prime_base^y
    
    mpz_mul_ui(x, prime_base, y);
    mpz_mod(x, x, modulus);
}



void generateKeys(mpz_t p_nu[], gmp_randstate_t state){
    mpz_t s_nu_mu[MAX_USER_NUM][MAX_USER_NUM];
    for(int i = 0; i < MAX_USER_NUM; i++){
        for(int j = 0; j < MAX_USER_NUM; j++){
            mpz_init(s_nu_mu[i][j]);
            if(i == j){
                mpz_set_ui(s_nu_mu[i][j], 0);
                continue;
            }
            mpz_urandomb(s_nu_mu[i][j], state, 1000);
        }
    }
    
    mpz_t p_nu_mu[MAX_USER_NUM][MAX_USER_NUM];
    for(int i = 0; i < MAX_USER_NUM; i++){
        for(int j = 0; j < MAX_USER_NUM; j++){
            mpz_init(p_nu_mu[i][j]);
            mpz_sub(p_nu_mu[i][j], s_nu_mu[i][j], s_nu_mu[j][i]);
        }
    }
    
    for(int i = 0; i < MAX_USER_NUM; i++){
        mpz_init(p_nu[i]);
        mpz_set_ui(p_nu[i], 0);
        for(int j = 0; j < MAX_USER_NUM; j++){
            mpz_clear(s_nu_mu[i][j]);
            mpz_add(p_nu[i], p_nu[i], p_nu_mu[i][j]);
            mpz_clear(p_nu_mu[i][j]);
        }
    }
}


void init3Darray(mpz_t ****arr){
    *arr = (mpz_t***)malloc(MAX_USER_NUM * sizeof(mpz_t**));
    //printf("%x\n", *arr);
    for(int i = 0; i < MAX_USER_NUM; i++){
        (*arr)[i] = (mpz_t**)malloc(MAX_DIMENSION * sizeof(mpz_t*));
        for(int j = 0; j < MAX_DIMENSION; j++){
            (*arr)[i][j] = (mpz_t*)malloc(MAX_DIMENSION * sizeof(mpz_t));
            for(int k = 0; k < MAX_DIMENSION; k++){
                mpz_init((*arr)[i][j][k]);
            }
            
        }
    }
}

void clear3Darray(mpz_t ****arr){
    for(int i = 0; i < MAX_USER_NUM; i++){
        for(int j = 0; j < MAX_DIMENSION; j++){
            for(int k = 0; k < MAX_DIMENSION; k++){
                mpz_clear((*arr)[i][j][k]);
            }
            free((*arr)[i][j]);
        }
        free((*arr)[i]);
    }
    free(*arr);
}


int main(int argc, char **argv){
    
    // This program uses a single machine to simulate multi-user secure aggregation
    

    int user_num[4] = {1, 2, 5, MAX_USER_NUM};
    int dimension[5] = {32, 64, 128, 256, MAX_DIMENSION};
    
    srand (time ( NULL));
    
    double aggregate_gradient[MAX_DIMENSION][MAX_DIMENSION] = {0};
    
    //double user_gradient[MAX_USER_NUM][MAX_DIMENSION][MAX_DIMENSION];
    double*** user_gradient;
    user_gradient = (double***)malloc(MAX_USER_NUM * sizeof(double**));
    for(int i = 0; i < MAX_USER_NUM; i++){
        user_gradient[i] = (double**)malloc(MAX_DIMENSION * sizeof(double*));
        for(int j = 0; j < MAX_DIMENSION; j++){
            user_gradient[i][j] = (double*)malloc(MAX_DIMENSION * sizeof(double));
            for(int k = 0; k < MAX_DIMENSION; k++){
                user_gradient[i][j][k] = (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1
            }
        }
    }
    
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    
    
    mpz_t p,q,N,N2,p_inverse_mod_q,q_inverse_mod_p;
    mpz_inits(p,q,N,N2,p_inverse_mod_q,q_inverse_mod_p,NULL);
    mpz_urandomb(p, state, 1024);
    mpz_nextprime(p, p);
    mpz_urandomb(q, state, 1024);
    mpz_nextprime(q, q);
    mpz_invert(p_inverse_mod_q, p, q);
    mpz_invert(q_inverse_mod_p, q, p);
    mpz_mul(N, p, q);
    mpz_mul(N2, N, N);
    
    
    mpz_t p_nu[MAX_USER_NUM];
    generateKeys(p_nu, state);
    
    //mpz_t x_nu[MAX_USER_NUM][MAX_DIMENSION][MAX_DIMENSION];
    mpz_t*** x_nu;
    init3Darray(&x_nu);
    
    
    //mpz_t y_nu[MAX_USER_NUM][MAX_DIMENSION][MAX_DIMENSION];
    mpz_t*** y_nu;
    init3Darray(&y_nu);
    
    
    mpz_t prime_base; // this is the prime base for the hash implementation
    mpz_init(prime_base);
    mpz_urandomb(prime_base, state, 1024);
    mpz_nextprime(prime_base, prime_base);
    
    
    
    // the following portion simulates the individual users who generate y_nu=(1+N)^{[x_nu]}H(t)^{p_nu} mod N^2
    time_t current_time = time(NULL);
    mpz_t temp;
    mpz_init(temp);
    for(int i = 0; i < MAX_USER_NUM; i++){
        for(int j = 0; j < MAX_DIMENSION; j++){
            for(int k = 0; k < MAX_DIMENSION; k++){
                printf("%d,%d,%d\n", i, j, k);
                encode(x_nu[i][j][k], user_gradient[i][j][k]);
                mpz_mul(temp, N, x_nu[i][j][k]);
                mpz_add_ui(y_nu[i][j][k], temp, 1);
                // y_nu[i][j][k] = (1+N)^{x_nu[i][j][k]}
                simpleHash(temp, current_time, prime_base, N); // temp = H(t)
                //mpz_powm(temp, temp, p_nu[i], N); // temp = H(t)^{p_nu[i]} mod N^2
                mpz_mul(y_nu[i][j][k], y_nu[i][j][k], temp);
                mpz_mod(y_nu[i][j][k], y_nu[i][j][k], N2); // y_nu[i][j][k] = final ciphertext to be broadcast
                
            }
        }
    }
    
    
    // the following portion simulates the aggregation at the parameter server side, who computes sum of user_gradients for each weight.
    
    clock_t begin, end;
    double time_spent;
    
    begin = clock();
    for(int j = 0; j< MAX_DIMENSION; j++){
        for(int k = 0; k < MAX_DIMENSION; k++){
            //printf("%d,%d\n", j, k);
            for(int i = 0; i < MAX_USER_NUM; i++){
                mpz_set_ui(temp, 1);
                mpz_mul(temp, temp, y_nu[i][j][k]);
            }
            mpz_mod(temp, temp, N2);
            mpz_sub_ui(temp, temp, 1);
            mpz_tdiv_q(temp, temp, N);
            decode(&(aggregate_gradient[j][k]), temp);
            
        }
    }
    end = clock();
    
    time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);

    printf("Total time: %f\n", time_spent * 1000);
    
    mpz_clear(temp);

    return 0;



    
    FILE *f;
	f = fopen("simulation.csv", "w");
	fprintf(f, "num,user1,user2,others,aggregator\n");

    //clock_t begin, end;
    //double time_spent;

	for(int i = 1; i <= 100; i++){
	
		begin = clock();

		

		end = clock();


		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);

		fprintf(f, "%d,%f,", i, time_spent * 1000);

		begin = clock();
		
        

		end = clock();

		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);

		
		fprintf(f, "%f,", time_spent * 1000);

		begin = clock();
		
        

		end = clock();
		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);

		fprintf(f, "%f,", time_spent * 1000);
		
		
		begin = clock();
		
		
        
		
		end = clock();
		
		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);
		
		fprintf(f, "%f\n", time_spent * 1000);

		fflush(f);	
	}

	fclose(f);
	

	

	gmp_randclear(state);
 
    mpz_clears(p,q,N,N2,p_inverse_mod_q,q_inverse_mod_p,NULL);
    mpz_clear(prime_base);
    
    for(int i = 0; i < MAX_USER_NUM; i++){
        for(int j = 0; j < MAX_DIMENSION; j++){
            free(user_gradient[i][j]);
        }
        free(user_gradient[i]);
    }
    free(user_gradient);
    
    
    
    clear3Darray(&x_nu);
    clear3Darray(&y_nu);

    
    return 0;
 
}


