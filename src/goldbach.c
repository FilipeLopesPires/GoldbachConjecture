#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>



int main() {
	/* contas com mpz */
	//mpz_add_ui(n,n,1); // n = n+1  // mpz_add(), mpz_sub(), mpz_sub_ui()
	//mpz_mul(n,n,n);	 // n = n*n  // mpz_mul_ui()
	//mpz_probab_prime_p(n,20) // determina se n é primo (2=sim;1=provavel;0=nao), 20 = cte (qt maior, menos 1's)
	//mpz_nextprime(x,n) // x = primo logo a seguir a n

	typedef unsigned long long u64;
	long long OFFSET=50000;

	mpz_t* dif = (mpz_t*)malloc(sizeof(mpz_t)*(10));

	//mpz_t PARAPRIMOS;  //inicio primos
	//mpz_init(PARAPRIMOS);
	//mpz_set_str(PARAPRIMOS, "", 10);

	/* variáveis iniciais */
	long long m=pow(10, 8);	//intervalo / tamanho array
	//printf("%llu \n",m);
	mpz_t n;    // início dos números pares
	int flag;   // flag de verificação de erros
	/* inicialização do n */
	mpz_init(n);
	flag = mpz_set_str(n, "4000000000000000000", 10);  //4*10^18
	assert (flag == 0); // if flag != 0, mpz_set falhou
	/* impressão do número */
	fprintf(stderr, "n = "); mpz_out_str(stderr,10,n); fprintf(stderr, "\n");
	
	/* memória reservada para os números pares da gama pretendida */
	mpz_t* pares=(mpz_t*)malloc(sizeof(mpz_t)*(m));
	for(int a=0; a<(m); a++) {mpz_init(pares[a]);}
	/* memória reservada para os números pares à medida que são encontrados */
	u64* found=(u64*)malloc(sizeof(u64)*(m));
	for(int a=0; a<(m); a++) {found[a]=0;}
	/* memória reservada para os números primos (de 0 a 50 000) 
	*/
	u64* little_primos=(u64*)malloc(sizeof(u64)*(OFFSET/5));
	for(int a=0; a<(OFFSET/5); a++) {little_primos[a]=0;}
	//mpz_t* primos=(mpz_t*)malloc(sizeof(mpz_t)*(m));
	//for(int a=0; a<(m); a++){mpz_init(primos[a]);}

	/* calculo dos números primos da gama [n-50000;n[ */
	int countprimos=0;
	mpz_t paux;
	mpz_init(paux);
	little_primos[countprimos]=2; countprimos++;
	for(int a=3; a<OFFSET; a+=2) {
		//if(a==1){little_primos[countprimos]=2; countprimos++; continue;}
		mpz_set_ui(paux, a);
		if(mpz_probab_prime_p(paux, 20)==2 || mpz_probab_prime_p(paux, 20)==1) { 
			little_primos[countprimos]=a;
			countprimos++;
		}
	}

	/* teste à conjetura */
	mpz_t tmp, soma, aux, oldp;
	mpz_init(tmp); mpz_init(soma); mpz_init(aux); mpz_init(oldp);
	char N[10][20];
	strcat(N[0], "4000000000000000000");strcpy(N[1], "4000000000200000000");strcpy(N[2], "4000000000400000000");
	strcpy(N[3], "4000000000600000000");strcpy(N[4], "4000000000800000000");strcpy(N[5], "4000000001000000000");
	strcpy(N[6], "4000000001200000000");strcpy(N[7], "4000000001400000000");strcpy(N[8], "4000000001600000000");
	strcpy(N[9], "4000000001800000000");
	for(int z=0; z<10; z++){
		mpz_set_str(n, N[z], 10);

		/* calculo de todos os pares na gama de valores [n;n+2m] e criação do array de pares */
		mpz_t aux;
		mpz_init(aux);
		for(int a=0; a<2*m; a+=2) {
			mpz_add_ui(aux, n, a);
			mpz_set(pares[a/2],aux);
		}

		mpz_init_set_ui(dif[z], 0);

		int ctrl = 1;

		for(int a=(-OFFSET); a<2*m; a+=2) {
			if(a+1<0) {
				mpz_sub_ui(tmp, n, abs(a+1));
			} else {
				mpz_add_ui(tmp, n, (a+1));
			}

			if(mpz_probab_prime_p(tmp, 20)==2 || mpz_probab_prime_p(tmp, 20)==1) {

				if(ctrl==1){
					mpz_set(oldp, tmp);
					ctrl = 0;
					//gmp_printf("%Zd\n",oldp);
				}

				mpz_out_str(stderr, 10, tmp);
				fprintf(stderr, "\r");
				for(int b=0; b<countprimos; b++) {
					mpz_add_ui(soma, tmp, little_primos[b]);
					for(int c=0; c<m; c++) {
						if(mpz_cmp(soma, pares[c])==0) {
							if(found[c]==0) {
								found[c]=little_primos[b];
							} else {
								if(little_primos[b]<found[c]){
									found[c]=little_primos[b];
								}
							}
						}
						else if(mpz_cmp(soma, pares[c])<0) {break;}
					}
				}
				mpz_sub(aux, tmp, oldp);
				//gmp_printf("%Zd - %Zd = %Zd\n",tmp,oldp, aux);

				if(mpz_cmp(aux, dif[z])>0){
					mpz_set(dif[z], aux);
				}
				mpz_set(oldp, tmp);

			}
		}
		mpz_set_ui(tmp, 0);
		mpz_set_ui(soma, 0);

	}

	/* impressão dos primos encontrados */
	int mal=0;
	for(int a=0; a<m; a++) {
		if(found[a]==0){
			mal=1;
			break;
		}
	}
	if(mal){fprintf(stderr, "Não conseguimos comprovar a conjetura!\n");}else{fprintf(stderr, "Conjetura confirmada!\n");}

	for(int index=0; index<10; index++){
		gmp_printf("%Zd\n", dif[index]);
	}

	/* limpexa dos mpz_t Handles */
	mpz_clear(n);
	mpz_clear(aux);
	mpz_clear(paux);
	mpz_clear(tmp);
	mpz_clear(soma);

	/* libertação da memória */
	free(pares);
	free(found);
	//free(primos);
	free(little_primos);
	return 1;
}
