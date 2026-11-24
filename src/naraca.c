#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <gmp.h>

/*
 * Compilação:
 *   gcc -std=c99 -O3 primo_peneira_tempo.c -o primo_peneira_tempo -lgmp
 *
 * Uso:
 *   ./primo_peneira_tempo <bits> [-f|--filtro]
 *
 *   -f, --filtro  Ativa a peneira da lista (termina em 1,3,7,9).
 *                 Sem esta flag, o programa testa todos os ímpares.
 */

#include <string.h>

#include <gmp.h>
#include <stdint.h>

#define CAMINHO_CHECKPOINT "./checkpoint_naraca.txt"

/* Tabela de consulta (0 = descarta, 1 = aceita) para todos os finais de numeros candidatos */
static uint8_t tabela_terminacoes[1000] = {0};

/* Chame esta função UMA VEZ no início do main() para preencher a tabela */
void inicializar_peneiraOld(void) {
    static const unsigned int lista[] = {
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
	   951, 
	   991
	};
    size_t tamanho = sizeof(lista) / sizeof(lista[0]);
    for (size_t i = 0; i < tamanho; i++) {
        tabela_terminacoes[lista[i]] = 1;
    }
}

/* Peneira O(1) direta no mpz_t: zero strings, zero loops */
static inline int peneira_listaOld(const mpz_t candidato) {
    /* Extrai os últimos 3 dígitos (candidato mod 1000) em nanossegundos */
    unsigned long ultimos_3 = mpz_fdiv_ui(candidato, 1000);
    return tabela_terminacoes[ultimos_3];
}

/**
A Base Matemática

A razão é simples e direta: no sistema decimal, qualquer número que termine 
em 0, 2, 4, 6 ou 8 é par (divisível por 2), e qualquer número que termine 
em 0 ou 5 é divisível por 5.

Como os números primos (com exceção do 2 e do 5) não podem ser divisíveis 
por 2 nem por 5, eles obrigatoriamente precisam terminar em um dos quatro dígitos 
restantes: 1, 3, 7 ou 9 . Isso é apenas uma consequência da definição de número 
primo e da forma como escrevemos os números na base 10.

A Teoria Mais Profunda (Teorema de Dirichlet)
Existe um resultado muito mais forte que confirma que todos esses quatro 
finais são igualmente possíveis para números primos, não importa quão grande 
eles sejam.

O Teorema de Dirichlet sobre progressões aritméticas (provado em 1837) 
garante que, para quaisquer dois números inteiros e que sejam coprimos 
(sem divisores comuns além do 1), existem infinitos números primos que deixam resto 
quando divididos por m.
**/
static inline int peneira_lista(const mpz_t candidato) {
    unsigned long u = mpz_fdiv_ui(candidato, 10);
    return (u == 1 || u == 3 || u == 7 || u == 9);
}

/* Retorna o tempo em segundos (com precisão de nanossegundos) */
double tempo_atual(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
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

void gerar_candidato(mpz_t n, size_t nbytes) {
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) { perror("urandom"); return; }

    unsigned char *buf = malloc(nbytes);
    if (!buf) { fclose(f); return; }

    size_t lidos = fread(buf, 1, nbytes, f);
    fclose(f);
    if (lidos != nbytes) { free(buf); return; }

    mpz_import(n, nbytes, 1, 1, 0, 0, buf);
    free(buf);
    mpz_setbit(n, 0); /* torna ímpar */
}

static void showHelp (const char* prog) {
	printf("Uso: %s \n");
	printf("    [-b|--bitsLen <tam em bits>] [-f|--filtro] [-p|--progress <n.m segundos>] [-s|--seed <arq-gerado.bin>]\n", prog);
	printf("    -b, --bitsLen       A quantidade de bits de uma semente aleatorória\n");
	printf("    -f, --filtro        Ativa a peneira da lista (termina em 1,3,7,9)\n");
	printf("    -p, --progress N.M  Ativa a barra de progresso a cada N.M segundos \n");
	printf("    -s, --seed arquivo gerado com dd if=/dev/urandom...  \n\n");
	printf("# Gera a semente com dd (1250 bytes = 10000 bits) \n");
	printf("dd if=/dev/urandom of=semente.bin bs=1250 count=1 \n");
}


int main(int argc, char *argv[]) {
	char *nomeSeed = NULL;
	int semente_fixada = 0;
    int usar_filtro = 0;
    unsigned long bits = 0;
	int mostrar_progresso = 0;
    int encontrado = 0;
    unsigned long tentativas = 0;
    unsigned long rejeitados_peneira = 0;
    unsigned long rejeitados_teste = 0;
    unsigned long chamadas_teste1 = 0;
    unsigned long chamadas_teste7 = 0;
    double intervalo_progresso = 5.0; /* segundos */
	mpz_t candidato;
	mpz_init(candidato);
	
    /* Parsing de opções com getopt_long */
    static struct option long_options[] = {
        {"bitsLen" , required_argument , 0, 'b'},
        {"progress", required_argument , 0, 'p'},
		{"semente" , required_argument , 0, 's'},
        {"filtro"  , no_argument       , 0, 'f'},
        {"help"    , no_argument       , 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "b:s:p:fh", long_options, NULL)) != -1) {
        switch (opt) {
			case 's':
				semente_fixada = 1;
				nomeSeed = strdup(optarg); 
				break;			
            case 'b': 
				bits = strtoul(optarg, NULL, 10);
			    break;
			case 'p': 
			    mostrar_progresso = 1;
				intervalo_progresso = atof(optarg);
			    break;
            case 'f': usar_filtro = 1; break;
            case 'h':
				showHelp(argv[0]);
                return 0;
            default:
                fprintf(stderr, "Opção inválida. Use -h para ajuda.\n");
                return 1;
        }
    }
	
	if (!semente_fixada && bits == 0) {
		fprintf(stderr, "Erro: informe -b <bits> ou -s <arquivo>.\n\n");
		showHelp(argv[0]);
		return 2;
	}
	
	if ( semente_fixada  == 0 ) { 
		
		if (bits < 64) {
			fprintf(stderr, "Erro: use pelo menos 64 bits.\n");
			showHelp(argv[0]);
			return 3;
		}
		gerar_candidato(candidato, (bits + 7) / 8);
	} else {
		FILE *fs = fopen(nomeSeed, "rb");
		if (!fs) { perror(nomeSeed); return 5; }

		/* Descobre o tamanho do arquivo */
		fseek(fs, 0, SEEK_END);
		long tam = ftell(fs);
		fseek(fs, 0, SEEK_SET);

		unsigned char *buf = malloc(tam);
		if (fread(buf, 1, tam, fs) != (size_t)tam) { /* erro */ }
		fclose(fs);

		mpz_import(candidato, tam, 1, 1, 0, 0, buf);
		free(buf);

		mpz_setbit(candidato, 0); /* torna ímpar */
		mpz_setbit(candidato, 8 * tam - 1);

		bits = mpz_sizeinbase(candidato, 2); /* pega a qtd de bits do número importado */	
	}
	
    printf("Configuração:\n");
	if ( semente_fixada  == 1 ) { 
		printf("\tSemente do arquivo: %s\n", nomeSeed);
		printf("\tSemente fixada: %lu bits\n", bits);	
	} else { 
		printf("\tBits para geração candidato: %lu\n", bits);
	}
    printf("\tPeneira da lista: %s\n", usar_filtro ? "ATIVADA" : "DESATIVADA");
    printf("\tTeste de primalidade: mpz_probab_prime_p com reps=1 e 7\n\n");

    double t_inicio = tempo_atual();
    double t_ultimo_progresso = tempo_atual();
	
	if ( usar_filtro == 1) {
		double agora = 0;		
		while (!encontrado) {
			tentativas++;
			
			// Indicador de progresso
			agora = tempo_atual();
			if (mostrar_progresso && (agora - t_ultimo_progresso >= intervalo_progresso)){
				double decorrido = agora - t_inicio;
				double taxa = tentativas / decorrido;
				printf("\rTentativas: %lu | Peneira: %lu | Teste: %lu | %.0f tent/s | %.1f s",
					   tentativas, rejeitados_peneira, rejeitados_teste, taxa, decorrido);
				fflush(stdout);
				t_ultimo_progresso = agora;
			}

			if (!peneira_lista(candidato)) {
				rejeitados_peneira++;
				// gerar_candidato(candidato, (bits + 7) / 8);
				mpz_add_ui(candidato, candidato, 2UL);
				continue;
			}

			int resultado = mpz_probab_prime_p(candidato, 1);
			chamadas_teste1++;

			if (resultado == 0) {
				rejeitados_teste++;
				// gerar_candidato(candidato, (bits + 7) / 8);
				mpz_add_ui(candidato, candidato, 2UL);
				continue;
			}

			resultado = mpz_probab_prime_p(candidato, 7);
			chamadas_teste7++;

			if (resultado != 0 ) {
				encontrado = 1;
				break;
			}

			/* Se chegou aqui, reps=7 retornou 0 (composto) ou > 0 (provavelmente primo
			   mas não definitivo). Em ambos os casos, não é o que queremos:
			   incrementa +2 e tenta o próximo ímpar. */
			rejeitados_teste++;
			mpz_add_ui(candidato, candidato, 2UL);
		}
	} else {
		double agora = 0;
		while (!encontrado) {
			tentativas++;
			// Indicador de progresso
			agora = tempo_atual();
			if (mostrar_progresso && (agora - t_ultimo_progresso >= intervalo_progresso)){
	
				double decorrido = agora - t_inicio;
				double taxa = tentativas / decorrido;
				printf("\rTentativas: %lu | Peneira: %lu | Teste: %lu | %.0f tent/s | %.1f s",
					   tentativas, rejeitados_peneira, rejeitados_teste, taxa, decorrido);
				fflush(stdout);
				t_ultimo_progresso = agora;
			}

			int resultado = mpz_probab_prime_p(candidato, 1);
			chamadas_teste1++;

			if (resultado == 0) {
				/* Composto no reps=1: incrementa +2 e tenta o próximo ímpar */
				rejeitados_teste++;
				mpz_add_ui(candidato, candidato, 2UL);
				continue;
			}

			/* Passou no reps=1: testa com reps=7 */
			resultado = mpz_probab_prime_p(candidato, 7);
			chamadas_teste7++;

			if (resultado != 0 ) {
				/* Primo definitivo (só ocorre para números pequenos) */
				encontrado = 1;
				break;
			}

			/* Se chegou aqui, reps=7 retornou 0 (composto) ou 1 (provavelmente primo
			   mas não definitivo). Em ambos os casos, não é o que queremos:
			   incrementa +2 e tenta o próximo ímpar. */
			rejeitados_teste++;
			mpz_add_ui(candidato, candidato, 2UL);
		}

		
	}
	
	/* Encontrou candidato que passou em todos os testes */
	double t_fim = tempo_atual();
	double tempo_total = t_fim - t_inicio;
	
	char *str = mpz_get_str(NULL, 10, candidato);
	
	printf("\n\n=== CANDIDATO ENCONTRADO ===\n");
	printf("Tentativas totais: %lu\n", tentativas);
	printf("Descartados pela peneira: %lu\n", rejeitados_peneira);
	printf("Rejeitado pelo teste GMP_1: %lu\n", chamadas_teste1);
	printf("Rejeitado pelo teste GMP_7: %lu\n", chamadas_teste7);
	printf("Descartados pelo teste de primalidade: %lu\n", rejeitados_teste);
	printf("Últimos 3 dígitos: %s\n", str + strlen(str) - 3);
	printf("Número completo:\n%s\n", str);
	printf("Tá na lista como é esperado: %s\n", (peneira_lista(candidato)?"Sim":"Não"));

	printf("\n--- TEMPO TOTAL ---\n");
	printf("Tempo total de execução: %.6f segundos\n", tempo_total);
	printf("Tempo médio por tentativa: %.9f segundos\n", tempo_total / (double)tentativas);

	free(str);
	if ( nomeSeed != NULL) free(nomeSeed);
	encontrado = 1;

    mpz_clear(candidato);
    // gmp_randclear(state);
    return 0;
}