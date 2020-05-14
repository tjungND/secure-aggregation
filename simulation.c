#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <math.h>

#define TOP 10000
#define E_FIXED_POINT 32


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



int main(){
    
    // This program uses a single machine to simulate multi-user secure aggregation

	
    int user_num[4] = {1, 2, 5, 10};
    int dimension[5] = {32, 64, 128, 256, 512};
    double user_gradient[4][512][512];

	clock_t begin, end;
	double time_spent;


	mpz_t Random[TOP][10];
    
    
    srand (time ( NULL));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 512; j++){
            for(int k = 0; k < 512; k++){
                user_gradient[i][j][k] = (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1
            }
        }
    }
    
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    
    
    mpz_t p,q,N,N2;
    mpz_inits(p,q,N,N2,NULL);
    mpz_randomb(p, state, 1024);
    mpz_randomb(q, state, 1024);
    mpz_mul(N, p, q);
    mpz_mul(N2, N, N);
    
    mpz_t s_nu_mu[10][10];
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            mpz_init(s_nu_mu[i][j]);
            if(i == j){
                mpz_set_ui(s_nu_mu[i][j], 0);
                continue;
            }
            mpz_randomb(s_nu_mu[i][j], state, 2000);
        }
    }
    
    mpz_t p_nu_mu[10][10];
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            mpz_init(p_nu_mu[i][j]);
            mpz_sub(p_nu_mu[i][j], s_nu_mu[i][j], s_nu_mu[j][i]);
        }
    }
    mpz_t p_nu[10];
    for(int i = 0; i < 10; i++){
        mpz_init(p_nu[i]);
        mpz_set_ui(p_nu[i], 0);
        for(int j = 0; j < 10; j++){
            mpz_clear(s_nu_mu[i][j]);
            mpz_add(p_nu[i], p_nu_mu[i][j]);
            mpz_clear(p_nu_mu[i][j]);
        }
    }
    
    
    
    

	
	


	int i, j, n, repeat, k, prod;



	mpz_t temp;
	mpz_init(temp);

	n = 10000;
	
	repeat = 100;
    
    FILE *f;
	f = fopen("micro_enc_dec_run_time.csv", "w");
	fprintf(f, "num,user1,user2,others,aggregator\n");



	for(i = 1; i <= 100; i++){
	
		begin = clock();

		for(k = 0; k < repeat; k++){

			mpz_set_ui(temp, 1);
	
			for(j = 0; j < n -2; j++){
				mpz_mul(temp, temp, small);
				mpz_mod(temp, temp, small);
			}
			mpz_powm(temp, small, temp, large);

		}

		end = clock();


		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC * repeat);

		fprintf(f, "%d,%f,", i, time_spent * 1000);

		begin = clock();
		for(k = 0; k < repeat; k++){


			mpz_powm(temp, small, small, large);
	
			mpz_mul(temp, temp, small);
			mpz_mod(temp, temp, large);

		}

		end = clock();

		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC * repeat);

		
		fprintf(f, "%f,", time_spent * 1000);

		begin = clock();
		for(k = 0; k < repeat + 560; k++){


			mpz_powm_ui(temp, small, i, small);
			mpz_powm(temp, temp, small, small);
			mpz_mul(temp, temp, small);
			mpz_mod(temp, temp, small);


		}

		end = clock();
		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC * (repeat + 560));

		fprintf(f, "%f,", time_spent * 1000);
		
		
		begin = clock();
		
		for(k = 0; k < repeat + 150; k++){
			mpz_set_ui(temp, 1);
			for(j = 0; j < n-1; j++){
				mpz_mul(temp, temp, large);
				mpz_mod(temp, temp, large);
			}
		}
		
		
		end = clock();
		
		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC * (repeat + 150));
		
		fprintf(f, "%f\n", time_spent * 1000);

		fflush(f);	
	}

	fclose(f);
	

	

	for(i = 0; i < TOP; i++){
		for(j = 0; j < 10; j++){
			mpz_clear(Random[i][j]);
		}
	}

	gmp_randclear(state);
 
    return 0;
}


