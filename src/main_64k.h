#ifndef MAIN_64K_H
#define MAIN_64K_H

#define CAMINHO_ENTROPIA   "/cygdrive/d/Projetos/Mersenne/entropia_64k.bin"
#define CAMINHO_CHECKPOINT "/cygdrive/d/Projetos/Mersenne/checkpoint_64k.txt"
#define CAMINHO_PRIMO_TXT  "/cygdrive/d/Projetos/Mersenne/primo_64k.txt"

#define TAMANHO_BYTES   65536       /* 64 KB = 524.288 bits (~157.826 digitos) */
#define LIMITE_CRIVO    50000000UL  /* Filtra todos os primos ate 50 milhoes na RAM */
#define TAMANHO_BLOCO   50000UL     /* 50.000 impares (salto +100.000) por bloco */
#define MAX_BLOCOS      1000UL      /* Ate salto +100.000.000 */
#define NUM_THREADS     8           /* 8 nucleos do Core i7 */

#endif /* MAIN_64K_H */
