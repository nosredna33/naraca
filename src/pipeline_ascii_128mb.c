#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <gmp.h>

/* 128 MB = 134.217.728 digitos decimais */
#define TAMANHO_DIGITOS    134217728UL
#define NUM_THREADS        8

#define CAMINHO_SEMENTE    "/cygdrive/d/Projetos/Mersenne/semente_128mb.txt"
#define CAMINHO_CHECKPOINT "/cygdrive/d/Projetos/Mersenne/checkpoint_128mb.txt"
#define CAMINHO_PRIMO_TXT  "/cygdrive/d/Projetos/Mersenne/primo_128mb.txt"

/* CONSTANTE 1: Tabela O(1) das 400 terminacoes validas nos 3 ultimos digitos (impares != 5) */
static unsigned char CONST_FINAIS_3DIG[1000];

/* CONSTANTE 2A: 25 Primoriais < 2^32 cobrindo primos de 3 a 547 (sem o 5, ja filtrado no Estagio 1) */
static const uint64_t CONST_PRIMORIAIS[25] = {
    353892219681ULL % 4294967296ULL ? 646969323ULL : 646969323ULL, /* 3*7*11*13*17*19*23*29 */
    95092081ULL,   /* 31*37*41*43*47 */
    909687937ULL,  /* 53*59*61*67*71 */
    3223462537ULL, /* 73*79*83*89*97 */
    131637181ULL,  /* 101*103*107*109 */
    214712171ULL,  /* 113*127*131*137 */
    307144969ULL,  /* 139*149*151*157 */
    468544381ULL,  /* 163*167*173*179 */
    686632379ULL,  /* 181*191*193*197 */
    1004098667ULL, /* 199*211*223*227 */
    1314380927ULL, /* 229*233*239*241 */
    1848056089ULL, /* 251*257*263*269 */
    2305556287ULL, /* 271*277*281*283 */
    3333902267ULL, /* 293*307*311*313 */
    4158599947ULL, /* 317*331*337*347 */
    50275993ULL,   /* 349*353*359 */
    53684141ULL,   /* 367*373*379 */
    62597143ULL,   /* 383*389*397 */
    71124053ULL,   /* 401*409*419 */
    79850161ULL,   /* 421*431*433 */
    89807177ULL,   /* 439*443*449 */
    99640159ULL,   /* 457*461*463 */
    111819389ULL,  /* 467*479*487 */
    123027931ULL,  /* 491*499*503 */
    139724087ULL   /* 509*521*523 */
};

/* CONSTANTE 2B: Os 97 primos individuais de 3 a 523 (sem o 5) */
static const uint64_t CONST_PRIMOS[97] = {
    3,7,11,13,17,19,23,29, 31,37,41,43,47, 53,59,61,67,71, 73,79,83,89,97,
    101,103,107,109, 113,127,131,137, 139,149,151,157, 163,167,173,179,
    181,191,193,197, 199,211,223,227, 229,233,239,241, 251,257,263,269,
    271,277,281,283, 293,307,311,313, 317,331,337,347, 349,353,359,
    367,373,379, 383,389,397, 401,409,419, 421,431,433, 439,443,449,
    457,461,463, 467,479,487, 491,499,503, 509,521,523
};

static const int CONST_MAPA_PRIMORIAL[97] = {
    0,0,0,0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3,
    4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,7,7,
    8,8,8,8, 9,9,9,9, 10,10,10,10, 11,11,11,11,
    12,12,12,12, 13,13,13,13, 14,14,14,14, 15,15,15,
    16,16,16, 17,17,17, 18,18,18, 19,19,19, 20,20,20,
    21,21,21, 22,22,22, 23,23,23, 24,24,24
};

static void inicializar_constantes(void) {
    int i;
    memset(CONST_FINAIS_3DIG, 0, sizeof(CONST_FINAIS_3DIG));
    for (i = 1; i < 1000; i += 2) {
        if ((i % 5) != 0) {
            CONST_FINAIS_3DIG[i] = 1;
        }
    }
}

/* Soma um deslocamento par (salto) diretamente na string ASCII mantendo o comprimento fixo */
static void somar_salto_ascii(char *dst, const char *src, size_t len, unsigned long salto) {
    size_t i = len;
    unsigned long carry = salto;

    memcpy(dst, src, len + 1);
    while (i > 0 && carry > 0UL) {
        i--;
        unsigned long soma = (unsigned long)(dst[i] - '0') + (carry % 10UL);
        dst[i] = (char)('0' + (soma % 10UL));
        carry = (carry / 10UL) + (soma / 10UL);
    }
}

/* Pipeline deterministico: Estagio 1 (3 digitos finais) + Estagio 2 (Primos 3..547 sem o 5) */
static int pipeline_constantes_100_primos(const char *str, size_t len) {
    uint64_t restos[25];
    uint64_t resto_541 = 0ULL;
    uint64_t resto_547 = 0ULL;
    size_t i = 0;
    int k;

    /* Estagio 1: Consulta O(1) dos 3 ultimos digitos (elimina multiplos de 2 e 5) */
    int final_3 = (str[len - 3] - '0') * 100 +
                  (str[len - 2] - '0') * 10 +
                  (str[len - 1] - '0');
    if (!CONST_FINAIS_3DIG[final_3]) return 0;

    /* Estagio 2: Reducao modular nos 128 MB contra primos de 3 a 547 (sem repetir 2 e 5) */
    memset(restos, 0, sizeof(restos));

    while (i + 9 <= len) {
        uint64_t bloco = (uint64_t)(str[i]     - '0') * 100000000ULL +
                         (uint64_t)(str[i + 1] - '0') * 10000000ULL +
                         (uint64_t)(str[i + 2] - '0') * 1000000ULL +
                         (uint64_t)(str[i + 3] - '0') * 100000ULL +
                         (uint64_t)(str[i + 4] - '0') * 10000ULL +
                         (uint64_t)(str[i + 5] - '0') * 1000ULL +
                         (uint64_t)(str[i + 6] - '0') * 100ULL +
                         (uint64_t)(str[i + 7] - '0') * 10ULL +
                         (uint64_t)(str[i + 8] - '0');

        for (k = 0; k < 25; k++) {
            restos[k] = (restos[k] * 1000000000ULL + bloco) % CONST_PRIMORIAIS[k];
        }
        resto_541 = (resto_541 * 1000000000ULL + bloco) % 541ULL;
        resto_547 = (resto_547 * 1000000000ULL + bloco) % 547ULL;
        i += 9;
    }

    while (i < len) {
        uint64_t dig = (uint64_t)(str[i] - '0');
        for (k = 0; k < 25; k++) {
            restos[k] = (restos[k] * 10ULL + dig) % CONST_PRIMORIAIS[k];
        }
        resto_541 = (resto_541 * 10ULL + dig) % 541ULL;
        resto_547 = (resto_547 * 10ULL + dig) % 547ULL;
        i++;
    }

    if (resto_541 == 0ULL || resto_547 == 0ULL) return 0;

    for (k = 0; k < 97; k++) {
        if ((restos[CONST_MAPA_PRIMORIAL[k]] % CONST_PRIMOS[k]) == 0ULL) {
            return 0;
        }
    }

    return 1;
}

static void salvar_checkpoint(unsigned long proximo_salto, unsigned long testados_gmp, unsigned long avaliados_total) {
    FILE *fc = fopen(CAMINHO_CHECKPOINT, "w");
    if (fc != NULL) {
        fprintf(fc, "%lu %lu %lu\n", proximo_salto, testados_gmp, avaliados_total);
        fclose(fc);
    }
}

int main(void) {
    FILE *fp;
    char *semente_ascii = (char *)malloc(TAMANHO_DIGITOS + 1);
    unsigned long salto_inicial = 0UL;
    unsigned long testados_gmp = 0UL;
    unsigned long avaliados_total = 0UL;
    int retomando = 0;
    volatile int encontrado = 0;
    unsigned long proximo_salto_global = 0UL;
    unsigned long maior_salto_concluido = 0UL;
    double t_ini = omp_get_wtime();

    if (semente_ascii == NULL) {
        printf("Erro: memoria insuficiente para alocar a semente de 128 MB.\n");
        return 1;
    }

    inicializar_constantes();

    /* 1. Tenta carregar savepoint existente (semente + checkpoint) */
    fp = fopen(CAMINHO_SEMENTE, "rb");
    if (fp != NULL) {
        if (fread(semente_ascii, 1, TAMANHO_DIGITOS, fp) == TAMANHO_DIGITOS) {
            semente_ascii[TAMANHO_DIGITOS] = '\0';
            retomando = 1;
        }
        fclose(fp);
    }

    if (retomando) {
        fp = fopen(CAMINHO_CHECKPOINT, "r");
        if (fp != NULL) {
            if (fscanf(fp, "%lu %lu %lu", &salto_inicial, &testados_gmp, &avaliados_total) < 1) {
                salto_inicial = 0UL;
            }
            fclose(fp);
        }
        printf("[SAVEPOINT CARREGADO] Retomando 128 MB do Salto +%lu (GMP testados: %lu | Avaliados: %lu)\n",
               salto_inicial, testados_gmp, avaliados_total);
    } else {
        /* 2. Primeiro inicio: sorteia a semente ASCII de 128 MB e grava no disco */
        FILE *frand = fopen("/dev/urandom", "rb");
        unsigned char *rnd_buf = (unsigned char *)malloc(TAMANHO_DIGITOS);
        size_t i;

        if (frand == NULL || rnd_buf == NULL) {
            printf("Erro ao abrir /dev/urandom ou alocar buffer.\n");
            free(semente_ascii);
            return 1;
        }
        if (fread(rnd_buf, 1, TAMANHO_DIGITOS, frand) != TAMANHO_DIGITOS) {
            fclose(frand);
            free(rnd_buf);
            free(semente_ascii);
            return 1;
        }
        fclose(frand);

        semente_ascii[0] = '1' + (rnd_buf[0] % 9);
        for (i = 1; i < TAMANHO_DIGITOS - 1; i++) {
            semente_ascii[i] = '0' + (rnd_buf[i] % 10);
        }
        /* Garante ultimo digito impar na base para avancar sempre de +2 em +2 */
        semente_ascii[TAMANHO_DIGITOS - 1] = '1' + ((rnd_buf[TAMANHO_DIGITOS - 1] % 5) * 2);
        semente_ascii[TAMANHO_DIGITOS] = '\0';
        free(rnd_buf);

        fp = fopen(CAMINHO_SEMENTE, "wb");
        if (fp != NULL) {
            fwrite(semente_ascii, 1, TAMANHO_DIGITOS, fp);
            fclose(fp);
        }
        salvar_checkpoint(0UL, 0UL, 0UL);
        printf("[NOVA SEMENTE 128 MB SALVA] Iniciando busca do Salto +0 com %d threads...\n", NUM_THREADS);
    }
    fflush(stdout);

    proximo_salto_global = salto_inicial;
    maior_salto_concluido = salto_inicial;

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        char *candidato_local = (char *)malloc(TAMANHO_DIGITOS + 1);
        mpz_t cand_mpz;
        unsigned long meu_salto;
        double t_teste;

        if (candidato_local == NULL) {
            #pragma omp critical
            printf("[Thread %d] Erro ao alocar buffer local de 128 MB!\n", tid);
        } else {
            mpz_init(cand_mpz);

            while (!encontrado) {
                /* Pega o proximo salto impar (+2) de forma atomica */
                #pragma omp critical
                {
                    meu_salto = proximo_salto_global;
                    proximo_salto_global += 2UL;
                    avaliados_total++;
                }

                somar_salto_ascii(candidato_local, semente_ascii, TAMANHO_DIGITOS, meu_salto);

                /* Pipeline deterministico sobre constantes (3 digitos + primos 3..547) */
                if (!pipeline_constantes_100_primos(candidato_local, TAMANHO_DIGITOS)) {
                    #pragma omp critical
                    {
                        if (!encontrado && meu_salto + 2UL > maior_salto_concluido && (avaliados_total % 100UL) == 0UL) {
                            maior_salto_concluido = meu_salto + 2UL;
                            salvar_checkpoint(maior_salto_concluido, testados_gmp, avaliados_total);
                        }
                    }
                    continue;
                }

                #pragma omp critical
                {
                    testados_gmp++;
                    printf("[Thread %d] Salto +%lu passou nas constantes -> Iniciando GMP 1 (GMP #%lu)...\n",
                           tid, meu_salto, testados_gmp);
                    fflush(stdout);
                }

                /* ESTAGIO GMP 1 */
                t_teste = omp_get_wtime();
                mpz_set_str(cand_mpz, candidato_local, 10);

                if (mpz_probab_prime_p(cand_mpz, 1) > 0) {
                    /* ESTAGIO GMP 7 */
                    if (mpz_probab_prime_p(cand_mpz, 7) > 0) {
                        #pragma omp critical
                        {
                            if (!encontrado) {
                                FILE *fout;
                                encontrado = 1;
                                printf("\n[PRIMO 128MB ENCONTRADO!] Salto vencedor: +%lu | Tempo: %.2fs\n",
                                       meu_salto, omp_get_wtime() - t_ini);
                                fout = fopen(CAMINHO_PRIMO_TXT, "w");
                                if (fout != NULL) {
                                    fprintf(fout, "%s\n", candidato_local);
                                    fclose(fout);
                                }
                                remove(CAMINHO_CHECKPOINT);
                            }
                        }
                    }
                } else {
                    #pragma omp critical
                    {
                        if (!encontrado) {
                            if (meu_salto + 2UL > maior_salto_concluido) {
                                maior_salto_concluido = meu_salto + 2UL;
                            }
                            salvar_checkpoint(maior_salto_concluido, testados_gmp, avaliados_total);
                            printf("[Thread %d] GMP 1 descartou Salto +%lu em %.2fs | Savepoint atualizado (+%lu)\n",
                                   tid, meu_salto, omp_get_wtime() - t_teste, maior_salto_concluido);
                            fflush(stdout);
                        }
                    }
                }
            }
            mpz_clear(cand_mpz);
            free(candidato_local);
        }
    }

    free(semente_ascii);
    return 0;
}
