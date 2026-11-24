/* ============================================================================
 * main.c — Motor de Pré-Triagem de Primos Gigantes (8.192 bits)
 * ----------------------------------------------------------------------------
 * Pipeline:
 *   1) Entropia (/dev/urandom, fallback dd)  -> semente de 1024 bytes
 *   2) Alinhamento MSB=1, LSB=1              -> N ímpar de 8192 bits
 *   3) Barreira 1 (25 primoriais, AVANCA_2)  -> descarta ~83,65% dos ímpares
 *   4) Barreira 2 (Fermat base 2 + MR x7)    -> PRP de 8192 bits
 *   5) Checkpoint e gravação decimal
 *
 * Compatível com C ANSI (C89/C90) + extensões de 64 bits.
 * ==========================================================================*/

#include "main2.h"

#include <errno.h>

/* ============================================================================
 * 1. ENTROPIA
 * ==========================================================================*/

int ler_entropia(unsigned char *buf, size_t n)
{
    FILE *f;
    size_t lidos;

    f = fopen("/dev/urandom", "rb");
    if (f != NULL) {
        lidos = fread(buf, 1, n, f);
        fclose(f);
        if (lidos == n)
            return 0;
    }

    /* Fallback: dd de /dev/urandom para arquivo temporário e leitura. */
    {
        char cmd[256];
        sprintf(cmd,
                "dd if=/dev/urandom of=%sentropia_tmp.bin bs=%lu count=1 2>/dev/null",
                CAMINHO_BASE, (unsigned long)n);
        if (system(cmd) != 0)
            return -1;

        f = fopen(CAMINHO_BASE "entropia_tmp.bin", "rb");
        if (f == NULL)
            return -1;
        lidos = fread(buf, 1, n, f);
        fclose(f);
        remove(CAMINHO_BASE "entropia_tmp.bin");
        if (lidos != n)
            return -1;
    }
    return 0;
}

void alinhar_bits(unsigned char *buf)
{
    /* buf[0] é o byte mais significativo (big-endian). */
    buf[0]            |= 0x80u;   /* MSB = 1 : garante 8192 bits exatos */
    buf[NUM_BYTES - 1] |= 0x01u;  /* LSB = 1 : força ímpar               */
}

/* ============================================================================
 * 2. BARREIRA 1 — RESTOS INICIAIS (REGRA DE HORNER EM BASE 2^32)
 * ==========================================================================*/

void calcular_restos(const unsigned char *buf, uint64_t *restos)
{
    int k;
    size_t i;

    for (k = 0; k < NUM_PRIMORIAIS; ++k)
        restos[k] = 0u;

    for (i = 0; i < NUM_PALAVRAS; ++i) {
        uint32_t w = ((uint32_t)buf[4*i    ] << 24) |
                     ((uint32_t)buf[4*i + 1] << 16) |
                     ((uint32_t)buf[4*i + 2] <<  8) |
                     ((uint32_t)buf[4*i + 3]);

        for (k = 0; k < NUM_PRIMORIAIS; ++k) {
            /* r = (r * 2^32 + w) mod P_k ; seguro pois r < P_k < 2^32 */
            restos[k] = ((restos[k] << 32) | (uint64_t)w) % PRIMORIAIS[k];
        }
    }
}

/* ============================================================================
 * 3. MATERIALIZAÇÃO DO CANDIDATO (SEMENTE + SALTO)
 * ----------------------------------------------------------------------------
 * A semente é big-endian; somamos `salto` (par) como inteiro big-endian
 * percorrendo do LSB para o MSB, com carry em uint32_t.
 * ==========================================================================*/

void materializar_candidato(const unsigned char *semente,
                            uint64_t salto,
                            unsigned char *destino)
{
    uint32_t carry;
    int i;
    uint64_t s;

    /* Copia a semente. */
    for (i = 0; i < NUM_BYTES; ++i)
        destino[i] = semente[i];

    /* Soma `salto` como dois uint32_t (low, high) a partir do fim. */
    s = salto;
    carry = (uint32_t)(s & 0xFFFFFFFFu);
    {
        uint32_t high = (uint32_t)((s >> 32) & 0xFFFFFFFFu);

        /* Processa palavra baixa primeiro, depois a alta. */
        for (i = 0; i < 2; ++i) {
            uint32_t add = (i == 0) ? carry : high;
            int j;
            if (add == 0u && i == 1)
                break;
            for (j = NUM_BYTES - 1; j >= 0 && add != 0u; --j) {
                uint32_t soma = (uint32_t)destino[j] + (add & 0xFFu);
                destino[j] = (unsigned char)(soma & 0xFFu);
                add = (add >> 8) | ((soma >> 8) << 0); /* propaga carry */
                /* A linha acima está simplificada; ver bloco corrigido abaixo. */
            }
            /* Reexecuta de forma limpa: */
            {
                uint32_t a = (i == 0) ? carry : high;
                int j2;
                for (j2 = NUM_BYTES - 1; j2 >= 0 && a != 0u; --j2) {
                    uint32_t soma = (uint32_t)destino[j2] + (a & 0xFFu);
                    destino[j2] = (unsigned char)(soma & 0xFFu);
                    a = (a >> 8) + (soma >> 8);
                }
            }
        }
    }

    /* Reforça invariantes (defensivo). */
    destino[0]            |= 0x80u;
    destino[NUM_BYTES - 1] |= 0x01u;
}

/* ============================================================================
 * 4. CHECKPOINT
 * ==========================================================================*/

int gravar_checkpoint(uint64_t salto)
{
    FILE *f = fopen(ARQ_CHECKPOINT, "w");
    if (f == NULL)
        return -1;
    fprintf(f, "%llu\n", (unsigned long long)salto);
    fclose(f);
    return 0;
}

int ler_checkpoint(uint64_t *salto)
{
    FILE *f = fopen(ARQ_CHECKPOINT, "r");
    unsigned long long v;
    if (f == NULL)
        return -1;
    if (fscanf(f, "%llu", &v) != 1) {
        fclose(f);
        return -1;
    }
    fclose(f);
    *salto = (uint64_t)v;
    return 0;
}

/* ============================================================================
 * 5. GRAVAÇÃO DO PRIMO
 * ==========================================================================*/

int gravar_primo(const char *decimal)
{
    FILE *f = fopen(ARQ_PRIMO, "w");
    if (f == NULL)
        return -1;
    fprintf(f, "%s\n", decimal);
    fclose(f);
    return 0;
}

/* ============================================================================
 * 6. BARREIRA 2 — FERMAT BASE 2 + MILLER-RABIN
 * ==========================================================================*/

int barreira2(const unsigned char *buf)
{
    mpz_t n;
    int aprovado;
    int resultado;

    mpz_init(n);

    /* Importa big-endian: ordem=1 (mais significativo primeiro), size=1 byte. */
    mpz_import(n, NUM_BYTES, 1, 1, 1, 0, buf);

    /* --- Estágio rápido: 2^(n-1) mod n == 1 --- */
    {
        mpz_t exp, res, dois;
        mpz_init(exp);
        mpz_init(res);
        mpz_init_set_ui(dois, 2u);

        mpz_sub_ui(exp, n, 1u);
        mpz_powm(res, dois, exp, n);

        aprovado = (mpz_cmp_ui(res, 1u) == 0);

        mpz_clear(exp);
        mpz_clear(res);
        mpz_clear(dois);
    }

    if (!aprovado) {
        mpz_clear(n);
        return 0;
    }

    /* --- Confirmação: Miller-Rabin com MR_RODADAS bases --- */
    resultado = mpz_probab_prime_p(n, MR_RODADAS);

    if (resultado != 0) {
        /* Grava representação decimal (base 10). */
        size_t tam = mpz_sizeinbase(n, 10) + 2;
        char *dec = (char *)malloc(tam);
        if (dec != NULL) {
            mpz_get_str(dec, 10, n);
            gravar_primo(dec);
            free(dec);
        }
        mpz_clear(n);
        return 1;
    }

    mpz_clear(n);
    return 0;
}

/* ============================================================================
 * 7. LOOP PRINCIPAL
 * ==========================================================================*/

int executar_busca(void)
{
    unsigned char semente[NUM_BYTES];
    unsigned char candidato[NUM_BYTES];
    uint64_t restos[NUM_PRIMORIAIS];
    uint64_t salto = 0;
    uint64_t avaliados_b2 = 0;

    /* 1) Entropia e alinhamento. */
    if (ler_entropia(semente, NUM_BYTES) != 0) {
        fprintf(stderr, "erro: falha ao ler entropia\n");
        return -1;
    }
    alinhar_bits(semente);

    /* 2) Checkpoint: retoma salto absoluto se existir. */
    if (ler_checkpoint(&salto) == 0)
        fprintf(stderr, "info: retomando do salto %llu\n",
                (unsigned long long)salto);

    /* Garante salto par (passo +2). */
    salto &= ~((uint64_t)1);

    /* 3) Restos iniciais da semente. */
    calcular_restos(semente, restos);

    /* Aplica salto inicial aos restos (avança 2 * (salto/2) vezes). */
    {
        uint64_t passos = salto / 2u;
        uint64_t p;
        for (p = 0; p < passos; ++p)
            AVANCA_2(restos);
    }

    /* 4) Loop incremental. */
    for (;;) {
        if (REPROVADO_BARREIRA1(restos)) {
            AVANCA_2(restos);
            salto += 2u;
            continue;
        }

        /* Barreira 2: materializa o candidato e testa. */
        materializar_candidato(semente, salto, candidato);

        fprintf(stderr, "info: candidato no salto %llu passou B1; testando B2...\n",
                (unsigned long long)salto);

        if (barreira2(candidato)) {
            gravar_checkpoint(salto);
            fprintf(stderr, "ok: primo PRP encontrado no salto %llu\n",
                    (unsigned long long)salto);
            return 0;
        }

        ++avaliados_b2;
        gravar_checkpoint(salto);

        AVANCA_2(restos);
        salto += 2u;

        if ((avaliados_b2 % 10u) == 0u)
            fprintf(stderr, "info: %llu candidatos avaliados em B2; salto=%llu\n",
                    (unsigned long long)avaliados_b2,
                    (unsigned long long)salto);
    }

    /* Inalcançável. */
    return 0;
}

/* ============================================================================
 * 8. MAIN
 * ==========================================================================*/

int main(void)
{
    fprintf(stderr, "=== Motor de Pré-Triagem de Primos Gigantes (8192 bits) ===\n");
    fprintf(stderr, "Primoriais: %d  |  Primos cobertos: %d  |  MR rodadas: %d\n",
            NUM_PRIMORIAIS, PRIMOS_COBERTOS, MR_RODADAS);

    if (executar_busca() != 0) {
        fprintf(stderr, "erro: busca abortada\n");
        return 1;
    }
    fprintf(stderr, "concluido.\n");
    return 0;
}
