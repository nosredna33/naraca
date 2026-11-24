#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mini-gmp.h"
#include "main.h"

int main(void) {
    FILE *fp_urandom;
    FILE *fp_arquivo;
    FILE *fp_chk;
    size_t num_bytes;
    size_t limite_blocos;
    size_t i;
    unsigned long tentativa;
    unsigned long rejeitados_b1;
    unsigned long rejeitados_b2;
    unsigned long long passo_acumulado;
    unsigned long long passo_salvo;
    unsigned char *buffer;
    unsigned long long palavra32;
    unsigned long long b_atual;
    unsigned long long r1, r2, r3, r4, r5, r6, r7, r8, r9, r10;
    unsigned long long r11, r12, r13, r14, r15, r16, r17, r18, r19, r20;
    unsigned long long r21, r22, r23, r24, r25;
    mpz_t semente_mpz;
    mpz_t candidato;
    mpz_t exp_fermat;
    mpz_t base_dois;
    mpz_t res_fermat;
    clock_t inicio_total;
    clock_t t_rodada_ini;
    clock_t t_rodada_fim;
    double tempo_total;
    double tempo_rodada;
    int resultado_mr;
    int retomando;
    char *str_primo;

    inicio_total = clock();
    num_bytes = (size_t)TAMANHO_BYTES;
    limite_blocos = num_bytes - (num_bytes % 4);
    rejeitados_b1 = 0UL;
    rejeitados_b2 = 0UL;
    passo_acumulado = 0ULL;
    retomando = 0;

    buffer = (unsigned char *)malloc(num_bytes);
    if (buffer == NULL) {
        return 1;
    }

    /* Verifica se existe checkpoint anterior para continuar */
    fp_chk = fopen(CAMINHO_CHECKPOINT, "r");
    if (fp_chk != NULL) {
        if (fscanf(fp_chk, "%llu", &passo_salvo) == 1) {
            fp_arquivo = fopen(CAMINHO_ENTROPIA, "rb");
            if (fp_arquivo != NULL && fread(buffer, 1, num_bytes, fp_arquivo) == num_bytes) {
                passo_acumulado = passo_salvo;
                retomando = 1;
                printf("Retomando busca de 8192 bits do Salto +%llu!\n", passo_acumulado);
            }
            if (fp_arquivo != NULL) {
                fclose(fp_arquivo);
            }
        }
        fclose(fp_chk);
    }
    if (!retomando) {
        fp_urandom = fopen("/dev/urandom", "rb");
        if (fp_urandom != NULL) {
            if (fread(buffer, 1, num_bytes, fp_urandom) != num_bytes) {
                fclose(fp_urandom);
                free(buffer);
                return 1;
            }
            fclose(fp_urandom);
            printf("Nova semente de 8192 bits lida de /dev/urandom\n");
        } else {
            if (system(COMANDO_DD) != 0) {
                free(buffer);
                return 1;
            }
            fp_arquivo = fopen(CAMINHO_ENTROPIA, "rb");
            if (fp_arquivo == NULL || fread(buffer, 1, num_bytes, fp_arquivo) != num_bytes) {
                if (fp_arquivo != NULL) fclose(fp_arquivo);
                free(buffer);
                return 1;
            }
            fclose(fp_arquivo);
            printf("Nova semente de 8192 bits lida via dd\n");
        }

        /* Garante MSB = 1 (8192 bits exatos) e LSB = 1 (impar) */
        buffer[0] |= 0x80;
        buffer[num_bytes - 1] |= 0x01;

        /* Salva a semente base em entropia.bin para permitir checkpoint */
        fp_arquivo = fopen(CAMINHO_ENTROPIA, "wb");
        if (fp_arquivo != NULL) {
            fwrite(buffer, 1, num_bytes, fp_arquivo);
            fclose(fp_arquivo);
        }
    }

    mpz_init(semente_mpz);
    mpz_init(candidato);
    mpz_init(exp_fermat);
    mpz_init(base_dois);
    mpz_init(res_fermat);
    mpz_set_ui(base_dois, 2UL);
    mpz_import(semente_mpz, num_bytes, 1, 1, 1, 0, buffer);

    /* Aplica o salto restaurado do 
     * checkpoint (se houver) 
     * SEM zerar passo_acumulado e 
     * SEM alterar semente_mpz */
    if (passo_acumulado > 0ULL) {
        mpz_add_ui(candidato, semente_mpz, (unsigned long)passo_acumulado);
        mpz_export(buffer, NULL, 1, 1, 1, 0, candidato);
    }

    r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = 0ULL;
    r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = r20 = 0ULL;
    r21 = r22 = r23 = r24 = r25 = 0ULL;
    for (i = 0; i < limite_blocos; i += 4) {
        palavra32 = ((unsigned long long)buffer[i] << 24) |
                    ((unsigned long long)buffer[i + 1] << 16) |
                    ((unsigned long long)buffer[i + 2] << 8)  |
                    ((unsigned long long)buffer[i + 3]);

        r1  = ((r1  << 32) | palavra32) % P64_1;  r2  = ((r2  << 32) | palavra32) % P64_2;
        r3  = ((r3  << 32) | palavra32) % P64_3;  r4  = ((r4  << 32) | palavra32) % P64_4;
        r5  = ((r5  << 32) | palavra32) % P64_5;  r6  = ((r6  << 32) | palavra32) % P64_6;
        r7  = ((r7  << 32) | palavra32) % P64_7;  r8  = ((r8  << 32) | palavra32) % P64_8;
        r9  = ((r9  << 32) | palavra32) % P64_9;  r10 = ((r10 << 32) | palavra32) % P64_10;
        r11 = ((r11 << 32) | palavra32) % P64_11; r12 = ((r12 << 32) | palavra32) % P64_12;
        r13 = ((r13 << 32) | palavra32) % P64_13; r14 = ((r14 << 32) | palavra32) % P64_14;
        r15 = ((r15 << 32) | palavra32) % P64_15; r16 = ((r16 << 32) | palavra32) % P64_16;
        r17 = ((r17 << 32) | palavra32) % P64_17; r18 = ((r18 << 32) | palavra32) % P64_18;
        r19 = ((r19 << 32) | palavra32) % P64_19; r20 = ((r20 << 32) | palavra32) % P64_20;
        r21 = ((r21 << 32) | palavra32) % P64_21; r22 = ((r22 << 32) | palavra32) % P64_22;
        r23 = ((r23 << 32) | palavra32) % P64_23; r24 = ((r24 << 32) | palavra32) % P64_24;
        r25 = ((r25 << 32) | palavra32) % P64_25;
    }

    for (; i < num_bytes; i++) {
        b_atual = (unsigned long long)buffer[i];
        r1  = ((r1  << 8) | b_atual) % P64_1;  r2  = ((r2  << 8) | b_atual) % P64_2;
        r3  = ((r3  << 8) | b_atual) % P64_3;  r4  = ((r4  << 8) | b_atual) % P64_4;
        r5  = ((r5  << 8) | b_atual) % P64_5;  r6  = ((r6  << 8) | b_atual) % P64_6;
        r7  = ((r7  << 8) | b_atual) % P64_7;  r8  = ((r8  << 8) | b_atual) % P64_8;
        r9  = ((r9  << 8) | b_atual) % P64_9;  r10 = ((r10 << 8) | b_atual) % P64_10;
        r11 = ((r11 << 8) | b_atual) % P64_11; r12 = ((r12 << 8) | b_atual) % P64_12;
        r13 = ((r13 << 8) | b_atual) % P64_13; r14 = ((r14 << 8) | b_atual) % P64_14;
        r15 = ((r15 << 8) | b_atual) % P64_15; r16 = ((r16 << 8) | b_atual) % P64_16;
        r17 = ((r17 << 8) | b_atual) % P64_17; r18 = ((r18 << 8) | b_atual) % P64_18;
        r19 = ((r19 << 8) | b_atual) % P64_19; r20 = ((r20 << 8) | b_atual) % P64_20;
        r21 = ((r21 << 8) | b_atual) % P64_21; r22 = ((r22 << 8) | b_atual) % P64_22;
        r23 = ((r23 << 8) | b_atual) % P64_23; r24 = ((r24 << 8) | b_atual) % P64_24;
        r25 = ((r25 << 8) | b_atual) % P64_25;
    }

    for (tentativa = 1UL; tentativa <= MAX_TENTATIVAS; tentativa++) {
        if (OK_1(r1)   && OK_2(r2)   && OK_3(r3)   && OK_4(r4)   && OK_5(r5)   &&
            OK_6(r6)   && OK_7(r7)   && OK_8(r8)   && OK_9(r9)   && OK_10(r10) &&
            OK_11(r11) && OK_12(r12) && OK_13(r13) && OK_14(r14) && OK_15(r15) &&
            OK_16(r16) && OK_17(r17) && OK_18(r18) && OK_19(r19) && OK_20(r20) &&
            OK_21(r21) && OK_22(r22) && OK_23(r23) && OK_24(r24) && OK_25(r25)) {

            t_rodada_ini = clock();
            mpz_add_ui(candidato, semente_mpz, (unsigned long)passo_acumulado);
            mpz_sub_ui(exp_fermat, candidato, 1UL);
            mpz_powm(res_fermat, base_dois, exp_fermat, candidato);

            resultado_mr = 0;
            if (mpz_cmp_ui(res_fermat, 1UL) == 0) {
                /* Passou Base 2: confirma com 7 rodadas de Miller-Rabin 
                   resultado_mr = mpz_probab_prime_p(candidato, 7);
                   1 rodada rapida descarta 99,999999999999% dos compostos
				*/
                resultado_mr = mpz_probab_prime_p(candidato, 1);
                if (resultado_mr > 0) {
                   /* Apenas o candidato vencedor passa pelas 7 rodadas extras */
                   resultado_mr = mpz_probab_prime_p(candidato, 7);
                }
            }
            t_rodada_fim = clock();
            tempo_rodada = (double)(t_rodada_fim - t_rodada_ini) / CLOCKS_PER_SEC;

            if (resultado_mr > 0) {
                tempo_total = (double)(t_rodada_fim - inicio_total) / CLOCKS_PER_SEC;
                printf("\n[SUCESSO 8192-bit] Salto +%llu (dt: %.3fs | Total: %.2fs)\n",
                       passo_acumulado, tempo_rodada, tempo_total);
                mpz_export(buffer, NULL, 1, 1, 1, 0, candidato);
                fp_arquivo = fopen(CAMINHO_ENTROPIA, "wb");
                if (fp_arquivo != NULL) {
                    fwrite(buffer, 1, num_bytes, fp_arquivo);
                    fclose(fp_arquivo);
                }
                remove(CAMINHO_CHECKPOINT);
                str_primo = mpz_get_str(NULL, 10, candidato);
                if (str_primo != NULL) {
                    fp_arquivo = fopen(CAMINHO_PRIMO_TXT, "w");
                    if (fp_arquivo != NULL) {
                        fprintf(fp_arquivo, "%s\n", str_primo);
                        fclose(fp_arquivo);
                    }
                    free(str_primo);
                }

                break;
            }
            rejeitados_b2++;
            tempo_total = (double)(t_rodada_fim - inicio_total) / CLOCKS_PER_SEC;
            printf("[Salto +%llu] MR #%lu rejeitado | dt: %.3fs | Total: %.2fs\n",
                   passo_acumulado, rejeitados_b2, tempo_rodada, tempo_total);
            fflush(stdout);

            /* Grava checkpoint a cada 10 candidatos testados no MR */
            if ((rejeitados_b2 % 10UL) == 0UL) {
                fp_chk = fopen(CAMINHO_CHECKPOINT, "w");
                if (fp_chk != NULL) {
                    fprintf(fp_chk, "%llu\n", passo_acumulado + 2ULL);
                    fclose(fp_chk);
                }
            }
        } else {
            rejeitados_b1++;
        }

        passo_acumulado += 2ULL;
        AVANCA_2(r1, P64_1);   AVANCA_2(r2, P64_2);   AVANCA_2(r3, P64_3);
        AVANCA_2(r4, P64_4);   AVANCA_2(r5, P64_5);   AVANCA_2(r6, P64_6);
        AVANCA_2(r7, P64_7);   AVANCA_2(r8, P64_8);   AVANCA_2(r9, P64_9);
        AVANCA_2(r10, P64_10); AVANCA_2(r11, P64_11); AVANCA_2(r12, P64_12);
        AVANCA_2(r13, P64_13); AVANCA_2(r14, P64_14); AVANCA_2(r15, P64_15);
        AVANCA_2(r16, P64_16); AVANCA_2(r17, P64_17); AVANCA_2(r18, P64_18);
        AVANCA_2(r19, P64_19); AVANCA_2(r20, P64_20); AVANCA_2(r21, P64_21);
        AVANCA_2(r22, P64_22); AVANCA_2(r23, P64_23); AVANCA_2(r24, P64_24);
        AVANCA_2(r25, P64_25);
    }

    mpz_clear(semente_mpz); 
    mpz_clear(candidato); 
    mpz_clear(exp_fermat);
    mpz_clear(base_dois);   
    mpz_clear(res_fermat); 
    free(buffer);
    return 0;
}
