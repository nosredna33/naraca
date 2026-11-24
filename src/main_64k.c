#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>
#include "main_64k.h"

/* Gera todos os primos impares de 3 ate LIMITE_CRIVO usando Crivo de Eratostenes */
static size_t gerar_primos_pequenos(unsigned long **lista_primos) {
    unsigned char *eh_composto;
    unsigned long i, j, count, idx;

    eh_composto = (unsigned char *)calloc(LIMITE_CRIVO + 1UL, 1);
    if (eh_composto == NULL) return 0;

    for (i = 3UL; i * i <= LIMITE_CRIVO; i += 2UL) {
        if (!eh_composto[i]) {
            for (j = i * i; j <= LIMITE_CRIVO; j += 2UL * i) {
                eh_composto[j] = 1;
            }
        }
    }

    count = 0;
    for (i = 3UL; i <= LIMITE_CRIVO; i += 2UL) {
        if (!eh_composto[i]) count++;
    }

    *lista_primos = (unsigned long *)malloc(count * sizeof(unsigned long));
    if (*lista_primos == NULL) {
        free(eh_composto);
        return 0;
    }

    idx = 0;
    for (i = 3UL; i <= LIMITE_CRIVO; i += 2UL) {
        if (!eh_composto[i]) {
            (*lista_primos)[idx++] = i;
        }
    }

    free(eh_composto);
    return count;
}

int main(void) {
    FILE *fp;
    unsigned char *buffer;
    unsigned char *crivo_bloco;
    unsigned long *primos;
    size_t qtd_primos, idx_p;
    unsigned long passo_base, k_inicial, k, r, p, primeiro_k;
    unsigned long bloco, sobreviventes, testados_mr, maior_k_concluido;
    int retomando, encontrado;
    unsigned long salto_vencedor;
    double t_ini, t_bloco_ini;
    mpz_t semente_mpz, base_bloco_mpz;
    char *str_primo;

    t_ini = omp_get_wtime();
    qtd_primos = gerar_primos_pequenos(&primos);
    if (qtd_primos == 0) return 1;

    buffer = (unsigned char *)malloc(TAMANHO_BYTES);
    crivo_bloco = (unsigned char *)malloc(TAMANHO_BLOCO);
    if (buffer == NULL || crivo_bloco == NULL) return 1;

    passo_base = 0UL;
    k_inicial = 0UL;
    retomando = 0;

    fp = fopen(CAMINHO_CHECKPOINT, "r");
    if (fp != NULL) {
        if (fscanf(fp, "%lu %lu", &passo_base, &k_inicial) >= 1) {
            FILE *fe = fopen(CAMINHO_ENTROPIA, "rb");
            if (fe != NULL && fread(buffer, 1, TAMANHO_BYTES, fe) == TAMANHO_BYTES) {
                retomando = 1;
                printf("Retomando 64KB do Salto +%lu (indice k=%lu | %lu primos no crivo)!\n",
                       passo_base, k_inicial, (unsigned long)qtd_primos);
            }
            if (fe != NULL) fclose(fe);
        }
        fclose(fp);
    }

    if (!retomando) {
        fp = fopen("/dev/urandom", "rb");
        if (fp == NULL || fread(buffer, 1, TAMANHO_BYTES, fp) != TAMANHO_BYTES) {
            if (fp != NULL) fclose(fp);
            return 1;
        }
        fclose(fp);

        /* Garante 524.288 bits exatos (MSB = 1) e numero impar (LSB = 1) */
        buffer[0] |= 0x80;
        buffer[TAMANHO_BYTES - 1] |= 0x01;

        fp = fopen(CAMINHO_ENTROPIA, "wb");
        if (fp != NULL) {
            fwrite(buffer, 1, TAMANHO_BYTES, fp);
            fclose(fp);
        }
        printf("Nova semente 64KB (524.288 bits) gerada! (%lu primos no crivo)\n",
               (unsigned long)qtd_primos);
    }

    mpz_init(semente_mpz);
    mpz_init(base_bloco_mpz);
    mpz_import(semente_mpz, TAMANHO_BYTES, 1, 1, 1, 0, buffer);
    encontrado = 0;
    salto_vencedor = 0UL;

    for (bloco = 0UL; bloco < MAX_BLOCOS && !encontrado; bloco++) {
        t_bloco_ini = omp_get_wtime();
        mpz_add_ui(base_bloco_mpz, semente_mpz, passo_base);
        memset(crivo_bloco, 0, TAMANHO_BLOCO);

        for (idx_p = 0; idx_p < qtd_primos; idx_p++) {
            p = primos[idx_p];
            r = mpz_fdiv_ui(base_bloco_mpz, p);
            if (r == 0UL) {
                primeiro_k = 0UL;
            } else {
                primeiro_k = ((p - r) % 2UL == 0UL) ? ((p - r) / 2UL) : ((2UL * p - r) / 2UL);
            }
            for (k = primeiro_k; k < TAMANHO_BLOCO; k += p) {
                crivo_bloco[k] = 1;
            }
        }

        sobreviventes = 0UL;
        for (k = k_inicial; k < TAMANHO_BLOCO; k++) {
            if (!crivo_bloco[k]) sobreviventes++;
        }
        printf("[Bloco +%lu | k>=%lu] Sobreviventes (p <= 50M): %lu (crivo: %.2fs)\n",
               passo_base, k_inicial, sobreviventes, omp_get_wtime() - t_bloco_ini);
        fflush(stdout);

        testados_mr = 0UL;
        maior_k_concluido = k_inicial;

        #pragma omp parallel num_threads(NUM_THREADS)
        {
            mpz_t cand_local, exp_local, base2_local, res_local;
            unsigned long k_priv, salto_atual;
            int eh_primo;

            mpz_init(cand_local); mpz_init(exp_local);
            mpz_init(base2_local); mpz_init(res_local);
            mpz_set_ui(base2_local, 2UL);

            #pragma omp for schedule(dynamic, 1)
            for (k_priv = k_inicial; k_priv < TAMANHO_BLOCO; k_priv++) {
                if (encontrado || crivo_bloco[k_priv]) continue;

                salto_atual = passo_base + (2UL * k_priv);
                mpz_add_ui(cand_local, semente_mpz, salto_atual);
                mpz_sub_ui(exp_local, cand_local, 1UL);
                mpz_powm(res_local, base2_local, exp_local, cand_local);

                eh_primo = 0;
                if (mpz_cmp_ui(res_local, 1UL) == 0) {
                    eh_primo = mpz_probab_prime_p(cand_local, 7);
                }

                #pragma omp critical
                {
                    testados_mr++;
                    if (eh_primo > 0 && !encontrado) {
                        encontrado = 1;
                        salto_vencedor = salto_atual;
                        printf("\n[SUCESSO 524.288-bit!] Salto +%lu | Tempo: %.2fs\n",
                               salto_vencedor, omp_get_wtime() - t_ini);
                        str_primo = mpz_get_str(NULL, 10, cand_local);
                        if (str_primo != NULL) {
                            fp = fopen(CAMINHO_PRIMO_TXT, "w");
                            if (fp != NULL) { fprintf(fp, "%s\n", str_primo); fclose(fp); }
                            free(str_primo);
                        }
                        remove(CAMINHO_CHECKPOINT);
                    } else if (!encontrado) {
                        if (k_priv + 1UL > maior_k_concluido) {
                            maior_k_concluido = k_priv + 1UL;
                            fp = fopen(CAMINHO_CHECKPOINT, "w");
                            if (fp != NULL) {
                                fprintf(fp, "%lu %lu\n", passo_base, maior_k_concluido);
                                fclose(fp);
                            }
                        }
                        if ((testados_mr % 10UL) == 0UL) {
                            printf("  -> %lu/%lu testados (Salto +%lu | %.1fs)...\n",
                                   testados_mr, sobreviventes, salto_atual, omp_get_wtime() - t_ini);
                            fflush(stdout);
                        }
                    }
                }
            }
            mpz_clear(cand_local); mpz_clear(exp_local);
            mpz_clear(base2_local); mpz_clear(res_local);
        }

        if (!encontrado) {
            passo_base += (2UL * TAMANHO_BLOCO);
            k_inicial = 0UL;
            fp = fopen(CAMINHO_CHECKPOINT, "w");
            if (fp != NULL) { fprintf(fp, "%lu 0\n", passo_base); fclose(fp); }
        }
    }

    mpz_clear(semente_mpz); mpz_clear(base_bloco_mpz);
    free(primos); free(buffer); free(crivo_bloco);
    return 0;
}
