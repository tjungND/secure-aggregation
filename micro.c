#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <math.h>

#define TOP 10000
#define E_FIXED_POINT 32
#define 


void encoding(mpz_t x, double y){
    mpz_set_d(x, y * pow(2, E_FIXED_POINT));
}

void decoding(double *x, mpz_t y){
    
    mpf_set_default_prec(E_FIXED_POINT + 32);
    mpf_t temp;
    mpf_init(temp);
    mpf_set_z(temp, y);
    mpf_div_2exp(temp, temp, E_FIXED_POINT);
    *x = mpf_get_d(temp);
    mpf_clear(temp);
}



int main(){

	FILE *f;
    int user_num[4] = {1, 2, 5, 10};
    double user_gradient[4][DI]

	clock_t begin, end;
	double time_spent;


	mpz_t Random[TOP][10];
    
    
    srand (time ( NULL));
    (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1


	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, 123);
	


	int i, j, n, repeat, k, prod;

	mpz_t small, large;
	mpz_inits(small, large, NULL);

	

	mpz_urandomb(small, state, 1024);
	mpz_urandomb(large, state, 2048);


	mpz_t temp;
	mpz_init(temp);

	n = 10000;
	
	repeat = 100;

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


