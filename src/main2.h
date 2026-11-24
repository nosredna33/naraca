/* ============================================================================
 * main.h — Motor de Pré-Triagem de Primos Gigantes (8.192 bits)
 * ----------------------------------------------------------------------------
 * Cabeçalho central: primoriais de 64 bits, macros de filtragem, caminhos
 * absolutos e protótipos. Projetado para C ANSI (C89/C90) + extensões de
 * 64 bits, rodando no CxStudio (Android/ARM64) contra Mini-GMP.
 *
 * Autor: (auditoria aberta)
 * Licença: domínio público / use por sua conta e risco.
 *
 * Compilando no CygWin
 
 * clang -std=c89 -O2 -o ../motor_primos2 ./main2.c ./mini-gmp.c
 * chnod +x ./motor_primos
 * 
 * ==========================================================================*/

#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>   /* C99; em C89 puro, defina typedefs abaixo */
#include <time.h>

/* --- Fallback para compiladores sem <stdint.h> (C89 estrito) ------------- */
#ifndef UINT64_MAX
typedef unsigned long long uint64_t;
typedef unsigned int       uint32_t;
typedef unsigned char      uint8_t;
#endif

/* --- Mini-GMP ----------------------------------------------------------- */
#include "mini-gmp.h"

/* ============================================================================
 * 1. CONSTANTES GERAIS
 * ==========================================================================*/

#define NUM_BITS        8192
#define NUM_BYTES       (NUM_BITS / 8)          /* 1024 */
#define NUM_PALAVRAS    (NUM_BYTES / 4)         /* 256 palavras de 32 bits */
#define MR_RODADAS      7                       /* rodadas Miller-Rabin      */

#define CAMINHO_BASE    "../Mersenne/"
#define ARQ_ENTROPIA    CAMINHO_BASE "entropia.bin"
#define ARQ_CHECKPOINT  CAMINHO_BASE "checkpoint.txt"
#define ARQ_PRIMO       CAMINHO_BASE "primo_candidato.txt"

/* ============================================================================
 * 2. TABELA DE PRIMORIAIS (25 acumuladores, todos < 2^32)
 * ----------------------------------------------------------------------------
 * Cada P_k é o produto de um bloco de primos consecutivos, escolhido de modo
 * que P_k < 2^32.  A redução modular de N (lido em palavras de 32 bits) por
 * P_k preserva N mod q para todo primo q | P_k, permitindo descartar
 * candidatos sem aritmética de precisão arbitrária.
 *
 * Invariante crítica: r < P_k < 2^32  =>  (r << 32) | w  <  2^64  =>  sem
 * overflow em uint64_t durante a Regra de Horner.
 * ==========================================================================*/

#define NUM_PRIMORIAIS 25

static const uint64_t PRIMORIAIS[NUM_PRIMORIAIS] = {
    UINT64_C(3234846615), /* r1 : 3..29          */
    UINT64_C(95092081),   /* r2 : 31..47         */
    UINT64_C(909687937),  /* r3 : 53..71         */
    UINT64_C(3223462537), /* r4 : 73..97         */
    UINT64_C(131637181),  /* r5 : 101..109       */
    UINT64_C(214712171),  /* r6 : 113..137       */
    UINT64_C(307144969),  /* r7 : 139..157       */
    UINT64_C(468544381),  /* r8 : 163..179       */
    UINT64_C(686632379),  /* r9 : 181..197       */
    UINT64_C(1004098667), /* r10: 199..227       */
    UINT64_C(1314380927), /* r11: 229..241       */
    UINT64_C(1848056089), /* r12: 251..269       */
    UINT64_C(2305556287), /* r13: 271..283       */
    UINT64_C(3333902267), /* r14: 293..313       */
    UINT64_C(4158599947), /* r15: 317..347       */
    UINT64_C(50275993),   /* r16: 349..359       */
    UINT64_C(53684141),   /* r17: 367..379       */
    UINT64_C(62597143),   /* r18: 383..397       */
    UINT64_C(71124053),   /* r19: 401..419       */
    UINT64_C(79850161),   /* r20: 421..433       */
    UINT64_C(89807177),   /* r21: 439..449       */
    UINT64_C(99640159),   /* r22: 457..463       */
    UINT64_C(111819389),  /* r23: 467..487       */
    UINT64_C(123027931),  /* r24: 491..503       */
    UINT64_C(139724087)   /* r25: 509..523       */
};

/* Soma total de primos cobertos (informativo): 99 primos, 3..523 */
#define PRIMOS_COBERTOS 99

/* ============================================================================
 * 3. MACROS DE FILTRAGEM
 * ==========================================================================*/

/* Avanço incremental de +2 em todos os acumuladores, sem divisão.
 * Pré-condição: 0 <= r[k] < P[k].  Pós-condição: idem.
 * (r + 2 pode exceder P apenas uma vez, pois P >= 3.) */
#define AVANCA_2(r)                                       \
    do {                                                  \
        int _k;                                           \
        for (_k = 0; _k < NUM_PRIMORIAIS; ++_k) {         \
            (r)[_k] += 2u;                                \
            if ((r)[_k] >= PRIMORIAIS[_k])                \
                (r)[_k] -= PRIMORIAIS[_k];                \
        }                                                 \
    } while (0)

/* Retorna 1 se algum acumulador é zero (candidato divisível por algum primo
 * coberto), 0 caso contrário (candidato sobrevive à Barreira 1). */
#define REPROVADO_BARREIRA1(r)                            \
    ( ((r)[0]==0) || ((r)[1]==0) || ((r)[2]==0) ||        \
      ((r)[3]==0) || ((r)[4]==0) || ((r)[5]==0) ||        \
      ((r)[6]==0) || ((r)[7]==0) || ((r)[8]==0) ||        \
      ((r)[9]==0) || ((r)[10]==0)|| ((r)[11]==0)||        \
      ((r)[12]==0)|| ((r)[13]==0)|| ((r)[14]==0)||        \
      ((r)[15]==0)|| ((r)[16]==0)|| ((r)[17]==0)||        \
      ((r)[18]==0)|| ((r)[19]==0)|| ((r)[20]==0)||        \
      ((r)[21]==0)|| ((r)[22]==0)|| ((r)[23]==0)||        \
      ((r)[24]==0) )

/* ============================================================================
 * 4. PROTÓTIPOS
 * ==========================================================================*/

/* Entropia: preenche buf[NUM_BYTES]; retorna 0 em sucesso. */
int  ler_entropia(unsigned char *buf, size_t n);

/* Ajusta MSB e LSB do buffer (big-endian: buf[0] é o byte mais significativo). */
void alinhar_bits(unsigned char *buf);

/* Calcula os 25 restos iniciais de N (big-endian, 32 bits/palavra). */
void calcular_restos(const unsigned char *buf, uint64_t *restos);

/* Recalcula o buffer a partir da semente + salto (necessário antes de
 * instanciar mpz_t no candidato atual).  salto é par (passo +2). */
void materializar_candidato(const unsigned char *semente,
                            uint64_t salto,
                            unsigned char *destino);

/* Checkpoint: grava/lê o salto absoluto. */
int  gravar_checkpoint(uint64_t salto);
int  ler_checkpoint(uint64_t *salto);

/* Grava a representação decimal do primo em ARQ_PRIMO. */
int  gravar_primo(const char *decimal);

/* Executa a Barreira 2 (Fermat base 2 + Miller-Rabin) sobre candidato. */
int  barreira2(const unsigned char *buf);

/* Loop principal; retorna 0 em sucesso (primo encontrado e gravado). */
int  executar_busca(void);

#endif /* MAIN_H */