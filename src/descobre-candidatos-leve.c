/***
 * COMPILAR NO CYGWIN: 
 * cd /cygdrive/d/Projetos/Mersenne/src
 * gcc -O3 -march=native -fopenmp ./proth_42m.c -I/usr/local/include -L/usr/local/lib -o ./proth_42m.exe -lgmp
 *
 * RODAR:
 * cd ..
 * ./proth_42m.exe
 
***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>

#define N_EXP              140000000ULL
#define NUM_THREADS        8
#define TAM_FILA           256
#define BLOCO_URANDOM      4096   /* Le 4096 uint64_t por vez do /dev/urandom */
#define RODADAS_GMP_FINAL  7      /* GMP 7+ para confirmacao final */

/* Meta estatistica media de testes no Filtro 2 para k*2^140M + 1 (apos Filtro 1 de 3 digitos) */
#define META_MEDIA_GMP1    14556090.0

#define CAMINHO_SEMENTE    "/cygdrive/d/Projetos/Mersenne/semente_k.txt"
#define CAMINHO_CHECKPOINT "/cygdrive/d/Projetos/Mersenne/checkpoint_pipeline.txt"
#define CAMINHO_PRIMO_TXT  "/cygdrive/d/Projetos/Mersenne/primo_vencedor.txt"

/* Lista constante de radicais de 3 digitos permitidos no final ASCII */
static const uint16_t LISTA_RADICAIS[] = {
    01,   3,   7,   9,  
    11,  13,  17,  19,  
	     23,       29,
    31,       37,  
    41,  43,  47,  
    53,            59,  
    61,       67,  
    71,  73,       79,  
    83,            89,  
    97,            
   101, 103, 107, 109, 
	    113,
             127, 
   131,      137, 139, 
                  149, 
   151,      157, 
        163, 167, 
        173,      179, 
   181, 
   191, 193, 197, 199, 
   211, 223, 227, 229,
        233,      239, 
   241, 
   251,      257, 
        263,      269, 
   271,      277, 
   281, 283, 
   	    293, 
   		     307, 
   311, 313, 317, 
   331,      337, 
             347, 349,
        353,      359, 
             367, 
        373,      379, 
        383,      389, 
             397, 
   401,           409,
                  419, 
   421, 
   431, 433,      439, 
        443,      449,
             457, 
   461, 463, 467, 
   471,           479, 
             487, 
   491,           499, 
        503,      509, 
   511, 
   521,      527, 
   551, 
   591, 
             607, 
             647, 
   671, 
   711, 
             727, 
   751, 
   791, 
   871, 
             887, 
			 927, 
   951,      953, 
   991
};

static unsigned char MASCARA_3DIG[1000];
static uint32_t POW2_MOD_1000;

/* Fila 1: Filtro 1 (Leve) -> Filtro 2 (GMP 1) */
static uint64_t fila_gmp1[TAM_FILA];
static int f1_ini = 0, f1_fim = 0, f1_qtd = 0;

/* Fila 2: Filtro 2 (GMP 1) -> Filtro 3 (GMP 7+) */
static uint64_t fila_gmp7[TAM_FILA];
static int f2_ini = 0, f2_fim = 0, f2_qtd = 0;

static uint32_t pow_mod_u32(uint64_t base, uint64_t exp, uint32_t mod) {
    uint64_t res = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) res = (res * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return (uint32_t)res;
}

static void inicializar_constantes(void) {
    int i;
    memset(MASCARA_3DIG, 0, sizeof(MASCARA_3DIG));
    for (i = 0; i < 98; i++) {
        MASCARA_3DIG[LISTA_RADICAIS[i]] = 1;
    }
    POW2_MOD_1000 = pow_mod_u32(2ULL, N_EXP, 1000U);
}

/* Filtro 1: olha apenas os 3 ultimos digitos decimais contra a lista constante */
static int filtro_leve_3dig(uint64_t k) {
    uint32_t ultimos_3 = (uint32_t)(((k % 1000ULL) * (uint64_t)POW2_MOD_1000 + 1ULL) % 1000ULL);
    return MASCARA_3DIG[ultimos_3];
}

static void formatar_eta(double seg_restantes, char *buf, size_t buf_len) {
    time_t agora = time(NULL);
    time_t fim = agora + (time_t)seg_restantes;
    struct tm *tm_info = localtime(&fim);
    if (tm_info != NULL) {
        strftime(buf, buf_len, "%d/%m/%Y %H:%M", tm_info);
    } else {
        snprintf(buf, buf_len, "indisponivel");
    }
}

static void salvar_checkpoint(uint64_t ultimo_k, uint64_t sorteios_f1, uint64_t sobreviventes_f1,
                              uint64_t concluidos_g1, uint64_t concluidos_g7,
                              double tempo_acum_f1, double soma_duracao_g1, double tempo_parede_total) {
    FILE *fc = fopen(CAMINHO_CHECKPOINT, "w");
    if (fc != NULL) {
        fprintf(fc, "%llu %llu %llu %llu %llu %.6f %.4f %.2f\n",
                (unsigned long long)ultimo_k,
                (unsigned long long)sorteios_f1,
                (unsigned long long)sobreviventes_f1,
                (unsigned long long)concluidos_g1,
                (unsigned long long)concluidos_g7,
                tempo_acum_f1,
                soma_duracao_g1,
                tempo_parede_total);
        fclose(fc);
    }
}

int main(void) {
    FILE *fp;
    uint64_t k_retomada = 0ULL;
    uint64_t sorteios_f1 = 0ULL;
    uint64_t sobreviventes_f1 = 0ULL;
    uint64_t concluidos_g1 = 0ULL;
    uint64_t concluidos_g7 = 0ULL;
    double tempo_previo_f1 = 0.0;
    double soma_duracao_g1 = 0.0;
    double tempo_previo_parede = 0.0;
    volatile int encontrado = 0;
    double t_ini_parede;

    inicializar_constantes();

    /* Carrega savepoint detalhado se existir */
    fp = fopen(CAMINHO_CHECKPOINT, "r");
    if (fp != NULL) {
        unsigned long long ck_k = 0, ck_s1 = 0, ck_a1 = 0, ck_g1 = 0, ck_g7 = 0;
        double ck_tf1 = 0.0, ck_sg1 = 0.0, ck_tp = 0.0;
        if (fscanf(fp, "%llu %llu %llu %llu %llu %lf %lf %lf",
                   &ck_k, &ck_s1, &ck_a1, &ck_g1, &ck_g7, &ck_tf1, &ck_sg1, &ck_tp) >= 1 && ck_k > 0) {
            k_retomada = (uint64_t)ck_k;
            sorteios_f1 = (uint64_t)ck_s1;
            sobreviventes_f1 = (uint64_t)ck_a1;
            concluidos_g1 = (uint64_t)ck_g1;
            concluidos_g7 = (uint64_t)ck_g7;
            tempo_previo_f1 = ck_tf1;
            soma_duracao_g1 = ck_sg1;
            tempo_previo_parede = ck_tp;
            printf("[SAVEPOINT] Retomando k=%llu | Filtro1: %llu/%llu | Filtro2(GMP1): %llu | Filtro3(GMP%d): %llu | Tempo: %.1fs\n",
                   (unsigned long long)k_retomada,
                   (unsigned long long)sobreviventes_f1,
                   (unsigned long long)sorteios_f1,
                   (unsigned long long)concluidos_g1,
                   RODADAS_GMP_FINAL,
                   (unsigned long long)concluidos_g7,
                   tempo_previo_parede);
        }
        fclose(fp);
    }

    if (k_retomada == 0ULL) {
        FILE *frand = fopen("/dev/urandom", "rb");
        uint64_t rnd_val = 0ULL;
        if (frand != NULL) {
            if (fread(&rnd_val, sizeof(rnd_val), 1, frand) != 1) rnd_val = 1000000001ULL;
            fclose(frand);
        } else {
            rnd_val = 1000000001ULL;
        }
        k_retomada = (rnd_val % 1000000000000ULL) | 1ULL;
        if (k_retomada < 1000000001ULL) k_retomada += 1000000001ULL;

        fp = fopen(CAMINHO_SEMENTE, "w");
        if (fp != NULL) {
            fprintf(fp, "%llu\n", (unsigned long long)k_retomada);
            fclose(fp);
        }
        salvar_checkpoint(k_retomada, 0ULL, 0ULL, 0ULL, 0ULL, 0.0, 0.0, 0.0);
        printf("[INICIO] Semente k=%llu | Filtro 1: Thread 0 | Filtro 2 (GMP 1): Threads 1..6 | Filtro 3 (GMP %d+): Thread 7\n",
               (unsigned long long)k_retomada, RODADAS_GMP_FINAL);
    }
    fflush(stdout);

    t_ini_parede = omp_get_wtime();

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();

        /* =====================================================================
         * FILTRO 1 (THREAD 0): LEVE (3 DIGITOS) + NOVO URANDOM A CADA DESCARTE
         * ===================================================================== */
        if (tid == 0) {
             uint64_t *buf_rnd = (uint64_t *)malloc(BLOCO_URANDOM * sizeof(uint64_t));
            size_t idx_buf = BLOCO_URANDOM;
            uint64_t k_atual = k_retomada;
            int primeiro_passo = 1;

            if (frand == NULL || buf_rnd == NULL) {
                #pragma omp critical(io_lock)
                printf("[Thread 0 - Filtro 1] Erro ao abrir /dev/urandom ou alocar buffer!\n");
            } else {
                while (!encontrado) {
                    double t_f1_ini = omp_get_wtime();
                    uint64_t k_candidato;

                    if (primeiro_passo) {
                        k_candidato = k_atual;
                        primeiro_passo = 0;
                    } else {
                        /* Descarte anterior gera novo sorteio no /dev/urandom */
                        if (idx_buf >= BLOCO_URANDOM) {
                            if (fread(buf_rnd, sizeof(uint64_t), BLOCO_URANDOM, frand) != BLOCO_URANDOM) {
                                break;
                            }
                            idx_buf = 0;
                        }
                        k_candidato = (buf_rnd[idx_buf++] % 1000000000000ULL) | 1ULL;
                        if (k_candidato < 1000000001ULL) k_candidato += 1000000001ULL;
                    }

                    int sobreviveu = filtro_leve_3dig(k_candidato);
                    double dur_f1 = omp_get_wtime() - t_f1_ini;

                    #pragma omp critical(stats_lock)
                    {
                        sorteios_f1++;
                        tempo_previo_f1 += dur_f1;
                    }

                    if (!sobreviveu) {
                        /* Descartado no Filtro 1 -> sorteia novo urandom imediatamente */
                        continue;
                    }

                    /* Sobreviveu ao Filtro 1 -> Repassa para a fila do Filtro 2 (GMP 1) */
                    int inserido = 0;
                    while (!inserido && !encontrado) {
                        #pragma omp critical(fila1_lock)
                        {
                            if (f1_qtd < TAM_FILA) {
                                fila_gmp1[f1_fim] = k_candidato;
                                f1_fim = (f1_fim + 1) % TAM_FILA;
                                f1_qtd++;
                                sobreviventes_f1++;
                                inserido = 1;
                            }
                        }
                        if (!inserido) {
                            struct timespec ts = {0, 1000000L}; /* 1ms */
                            nanosleep(&ts, NULL);
                        }
                    }
                }
            }
            if (frand != NULL) fclose(frand);
            free(buf_rnd);
        }
        /* =====================================================================
         * FILTRO 2 (THREADS 1..6): GMP 1 -> REPASSA SOBREVIVENTES PARA GMP 7+
         * ===================================================================== */
        else if (tid >= 1 && tid <= 6) {
            mpz_t cand;
            mpz_init(cand);

            while (!encontrado) {
                uint64_t k_cand = 0ULL;
                int pegou = 0;

                #pragma omp critical(fila1_lock)
                {
                    if (f1_qtd > 0) {
                        k_cand = fila_gmp1[f1_ini];
                        f1_ini = (f1_ini + 1) % TAM_FILA;
                        f1_qtd--;
                        pegou = 1;
                    }
                }

                if (!pegou) {
                    struct timespec ts = {0, 1000000L};
                    nanosleep(&ts, NULL);
                    continue;
                }

                double t_g1 = omp_get_wtime();
                mpz_set_ui(cand, (unsigned long)k_cand);
                mpz_mul_2exp(cand, cand, (unsigned long)N_EXP);
                mpz_add_ui(cand, cand, 1UL);

                int sobreviveu_g1 = (mpz_probab_prime_p(cand, 1) > 0);
                double dur_g1 = omp_get_wtime() - t_g1;

                if (sobreviveu_g1) {
                    int inserido_f2 = 0;
                    #pragma omp critical(io_lock)
                    {
                        printf("[Filtro 2 - T%d] k=%llu passou na GMP 1 em %.2fs -> Repassando para Filtro 3 (GMP %d+)...\n",
                               tid, (unsigned long long)k_cand, dur_g1, RODADAS_GMP_FINAL);
                        fflush(stdout);
                    }
                    while (!inserido_f2 && !encontrado) {
                        #pragma omp critical(fila2_lock)
                        {
                            if (f2_qtd < TAM_FILA) {
                                fila_gmp7[f2_fim] = k_cand;
                                f2_fim = (f2_fim + 1) % TAM_FILA;
                                f2_qtd++;
                                inserido_f2 = 1;
                            }
                        }
                        if (!inserido_f2) {
                            struct timespec ts = {0, 1000000L};
                            nanosleep(&ts, NULL);
                        }
                    }
                } else {
                    #pragma omp critical(io_lock)
                    {
                        if (!encontrado) {
                            double tempo_parede_atual = tempo_previo_parede + (omp_get_wtime() - t_ini_parede);
                            double media_f1_us = (sorteios_f1 > 0) ? (tempo_previo_f1 * 1e6 / (double)sorteios_f1) : 0.0;
                            double media_g1_individual;
                            double ritmo_parede_g1;
                            double seg_restantes;
                            char eta_str[64];

                            concluidos_g1++;
                            soma_duracao_g1 += dur_g1;

                            salvar_checkpoint(k_cand, sorteios_f1, sobreviventes_f1, concluidos_g1, concluidos_g7,
                                              tempo_previo_f1, soma_duracao_g1, tempo_parede_atual);

                            media_g1_individual = soma_duracao_g1 / (double)concluidos_g1;
                            ritmo_parede_g1 = tempo_parede_atual / (double)concluidos_g1;
                            seg_restantes = (META_MEDIA_GMP1 - (double)concluidos_g1) * ritmo_parede_g1;
                            if (seg_restantes < 0.0) seg_restantes = 0.0;
                            formatar_eta(seg_restantes, eta_str, sizeof(eta_str));

                            printf("[Filtro 2 - T%d] Descartou k=%llu | F1: %.2fus/sorteio (%llu/%llu) | GMP1: %.2fs (media %.2fs) | Total GMP1: %llu | ETA: %s\n",
                                   tid, (unsigned long long)k_cand,
                                   media_f1_us, (unsigned long long)sobreviventes_f1, (unsigned long long)sorteios_f1,
                                   dur_g1, media_g1_individual,
                                   (unsigned long long)concluidos_g1, eta_str);
                            fflush(stdout);
                        }
                    }
                }
            }
            mpz_clear(cand);
        }
        /* =====================================================================
         * FILTRO 3 (THREAD 7): CONFIRMACAO FINAL GMP 7+
         * ===================================================================== */
        else if (tid == 7) {
            mpz_t cand;
            mpz_init(cand);

            while (!encontrado) {
                uint64_t k_cand = 0ULL;
                int pegou = 0;

                #pragma omp critical(fila2_lock)
                {
                    if (f2_qtd > 0) {
                        k_cand = fila_gmp7[f2_ini];
                        f2_ini = (f2_ini + 1) % TAM_FILA;
                        f2_qtd--;
                        pegou = 1;
                    }
                }

                if (!pegou) {
                    struct timespec ts = {0, 5000000L};
                    nanosleep(&ts, NULL);
                    continue;
                }

                double t_g7 = omp_get_wtime();
                mpz_set_ui(cand, (unsigned long)k_cand);
                mpz_mul_2exp(cand, cand, (unsigned long)N_EXP);
                mpz_add_ui(cand, cand, 1UL);

                int passou_g7 = (mpz_probab_prime_p(cand, RODADAS_GMP_FINAL) > 0);
                double dur_g7 = omp_get_wtime() - t_g7;

                #pragma omp critical(io_lock)
                {
                    double tempo_parede_atual = tempo_previo_parede + (omp_get_wtime() - t_ini_parede);
                    concluidos_g7++;
                    if (passou_g7 && !encontrado) {
                        FILE *fout;
                        encontrado = 1;
                        printf("\n[FILTRO 3 - PRIMO CONFIRMADO GMP %d+!] %llu * 2^%llu + 1 | Tempo GMP%d: %.2fs | Total: %.2fs\n",
                               RODADAS_GMP_FINAL, (unsigned long long)k_cand, N_EXP,
                               RODADAS_GMP_FINAL, dur_g7, tempo_parede_atual);
                        fout = fopen(CAMINHO_PRIMO_TXT, "w");
                        if (fout != NULL) {
                            fprintf(fout, "%llu * 2^%llu + 1\n", (unsigned long long)k_cand, N_EXP);
                            fclose(fout);
                        }
                        remove(CAMINHO_CHECKPOINT);
                    } else if (!encontrado) {
                        salvar_checkpoint(k_cand, sorteios_f1, sobreviventes_f1, concluidos_g1, concluidos_g7,
                                          tempo_previo_f1, soma_duracao_g1, tempo_parede_atual);
                        printf("[Filtro 3 - T7] Descartou k=%llu na GMP %d em %.2fs (Total GMP%d: %llu)\n",
                               (unsigned long long)k_cand, RODADAS_GMP_FINAL, dur_g7,
                               RODADAS_GMP_FINAL, (unsigned long long)concluidos_g7);
                        fflush(stdout);
                    }
                }
            }
            mpz_clear(cand);
        }
    }

    return 0;
}
