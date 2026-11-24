# README.md — Naraca: Busca de Primos com Pré-Filtro de Terminações

## Visão Geral

O **Naraca**, que pronuncia-se **"Na Raça"** é uma ferramenta experimental escrita em C (ANSI C99) para busca de números PROVAVEIS números primos de tamanho arbitrário, com foco em **pré-filtragem de candidatos** baseada nas terminações decimais possíveis de primos (1, 3, 7, 9). O projeto investiga o impacto real dessa peneira no desempenho da busca, usando a biblioteca **GMP** para aritmética de precisão arbitrária e testes de primalidade.

O objetivo não é competir com projetos como o GIMPS (que busca primos de Mersenne), mas sim **medir e documentar** o ganho (ou a ausência dele) de um pré-filtro simples em buscas de primos genéricos.

---

## Fundamentação Matemática

### 1. Terminações possíveis de primos

Todo número primo maior que 5 termina, na base 10, em um dos dígitos **1, 3, 7 ou 9**. Isso é uma consequência direta da definição de primo:

- Números terminados em **0, 2, 4, 6, 8** são pares (divisíveis por 2).
- Números terminados em **0 ou 5** são divisíveis por 5.

Como primos (exceto 2 e 5) não podem ser divisíveis por 2 nem por 5, eles obrigatoriamente terminam em 1, 3, 7 ou 9.

### 2. Teorema de Dirichlet

O **Teorema de Dirichlet sobre progressões aritméticas** (1837) garante que, para quaisquer inteiros \(a\) e \(m\) coprimos, existem infinitos primos congruentes a \(a \pmod{m}\). Aplicado a \(m = 10\), isso significa que existem infinitos primos terminados em 1, em 3, em 7 e em 9 — nenhuma dessas terminações é "privilegiada" ou "rara".

### 3. O que a peneira realmente faz

No código, o candidato já é gerado **ímpar** (bit 0 forçado a 1). Portanto, a peneira de terminações 1,3,7,9 rejeita apenas os candidatos terminados em **5** — ou seja, **20% dos ímpares**. Ela **não** rejeita pares (já eliminados) nem múltiplos de 3, 7, 11, etc.

### 4. Limitações do pré-filtro

A peneira reduz o número de chamadas ao teste de primalidade em 20%, mas **não reduz o tempo na mesma proporção**, porque os candidatos terminados em 5 são rejeitados muito rapidamente pelo `mpz_probab_prime_p` (divisão por 5 é uma operação desprezível). O ganho de tempo medido foi de **~1,4%** para primos de 20.000 bits.

Para ganhos maiores, seria necessário um pré-filtro mais agressivo (divisibilidade por 3, 7, 11, 13, ...), que elimina candidatos que o teste de primalidade demora mais para rejeitar.

---

## Funcionalidades

- Geração de candidatos aleatórios de \(N\) bits via `/dev/urandom`.
- Leitura de semente de um arquivo binário (`-s`), permitindo comparações honestas (mesma semente com e sem filtro).
- Pré-filtro opcional de terminações (1, 3, 7, 9) via `-f`.
- Teste de primalidade com `mpz_probab_prime_p` (reps=1 e reps=7).
- Barra de progresso configurável (`-p N.M` segundos).
- Relatório final com contadores, tempos e terminação do primo encontrado.

---

## Instalação do GMP 6.3.0 a partir do código-fonte (agnóstica de ambiente)

Esta seção descreve como baixar e compilar o GMP 6.3.0 manualmente. Os passos funcionam no **Cygwin**, **Linux** e **macOS** (com pequenas adaptações de prefixo). O método é o mesmo; apenas o comando de download e o local de instalação podem variar.

### Pré-requisitos

- Compilador C (GCC)
- `make`
- `curl`
- `xz` (ou `lzip`, dependendo do formato do arquivo)
- Ferramentas de build (autoconf, automake, libtool podem ser necessárias em alguns sistemas)

**No Cygwin:**
```bash
apt-cyg install gcc-core make curl xz
```

**No Debian/Ubuntu:**
```bash
sudo apt install build-essential curl xz-utils
```

**No macOS (Homebrew):**
```bash
brew install gcc make curl xz
```

### Passo a passo

**1. Baixar o código-fonte**

O GMP 6.3.0 está disponível em `.tar.xz`. Use `curl` para baixar:

```bash
curl -O https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
```

*Alternativa com `wget`:*
```bash
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
```

**2. Extrair o pacote**

```bash
tar -xf gmp-6.3.0.tar.xz
```

O comando `tar` moderno detecta automaticamente a compressão xz. Se necessário, use `tar -xJf gmp-6.3.0.tar.xz`.

**3. Entrar no diretório**

```bash
cd gmp-6.3.0
```

**4. Configurar (configure)**

O `./configure` detecta automaticamente o ambiente (Cygwin, Linux, macOS) e o compilador.

```bash
./configure --prefix=/usr/local --enable-cxx
```

- `--prefix=/usr/local` instala em `/usr/local` (padrão). No Cygwin, você pode usar um caminho dentro do seu home, por exemplo: `--prefix=$HOME/gmp`.
- `--enable-cxx` habilita a interface C++ (opcional, mas útil).
- Para biblioteca compartilhada (DLL no Cygwin), adicione `--disable-static --enable-shared`.

**5. Compilar**

```bash
make
```

**6. Executar os testes (altamente recomendado)**

A documentação do GMP enfatiza que o `make check` **não deve ser pulado**, porque o risco de compilação incorreta é real.

```bash
make check
```

**7. Instalar**

```bash
make install
```

No Cygwin, se você usou `--prefix=$HOME/gmp`, não precisa de `sudo`. Em `/usr/local`, pode ser necessário `sudo make install` (dependendo das permissões).

### Resumo dos comandos (copiar e colar)

```bash
curl -O https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --prefix=/usr/local --enable-cxx
make
make check
make install
```

### Verificação

Após a instalação, verifique se a biblioteca está acessível:

```bash
gmp-config --version
```

Ou compile o Naraca com `-lgmp` e execute um teste rápido.

### Notas específicas do Cygwin

- O Cygwin fornece um ambiente POSIX sobre Windows. O `./configure` do GMP detecta o Cygwin corretamente.
- Se ocorrer erro de `make` relacionado a "colon expected", consulte a documentação do Cygwin sobre Makefiles.
- O GMP compilado no Cygwin produz bibliotecas `.dll.a` e `.a`, compatíveis com o GCC do Cygwin.

### Notas para outros ambientes

- **Linux:** os passos são idênticos. O `sudo` pode ser necessário no `make install` se o prefixo for `/usr/local`.
- **macOS:** o Homebrew já fornece GMP (`brew install gmp`), mas a compilação manual funciona com `--prefix=$HOME/gmp` para evitar conflitos com o sistema.
- **Windows nativo (MinGW/MSYS2):** os passos são semelhantes, mas o prefixo deve ser um caminho Windows (ex.: `--prefix=C:/gmp`).

---

## Compilação

### Pré-requisitos

- GCC (ou compilador C99 compatível)
- Biblioteca GMP (instalada conforme seção anterior)
- Make (opcional, para automação)

### Compilação manual (recomendada)

O comando usado nos experimentos foi:

```bash
gcc -O3 -march=native -fopenmp ./naraca.c \
    -I/usr/local/include \
    -L/usr/local/lib \
    -lgmp \
    -o ../naraca.exe
```

**Explicação das flags:**

| Flag | Descrição |
|------|-----------|
| `-O3` | Otimização máxima do GCC (loop unrolling, inlining, vetorização). |
| `-march=native` | Gera código específico para a CPU da máquina (usa instruções AVX2, etc., se disponíveis). **Não portável** entre máquinas diferentes. |
| `-fopenmp` | Habilita OpenMP (paralelismo com threads). No código atual, não há diretivas `#pragma omp`, mas a flag está presente para futuras extensões. |
| `-I/usr/local/include` | Caminho para os headers do GMP instalado manualmente. |
| `-L/usr/local/lib` | Caminho para as bibliotecas do GMP instalado manualmente. |
| `-lgmp` | Linka a biblioteca GMP. |
| `-o ../naraca.exe` | Nome do executável de saída (no diretório pai, no Cygwin). |

### Compilação via Makefile

O Makefile incluso automatiza o processo:

```bash
make          # compila o binário naraca
make clean    # remove binário e arquivos objeto
make run      # executa um exemplo com 10000 bits e filtro
```

### Notas específicas do Cygwin

- O executável é gerado com extensão `.exe` (ex.: `../naraca.exe`).
- O Cygwin resolve `/usr/local` para o caminho correspondente dentro da instalação Cygwin (ex.: `C:\cygwin64\usr\local`).
- A flag `-march=native` funciona no Cygwin/GCC, mas o binário resultante **não é portável** para outras máquinas com CPU diferente.

### Notas para outros ambientes

- **Linux:** os mesmos comandos funcionam, apenas troque `../naraca.exe` por `../naraca` (sem extensão).
- **macOS:** o Homebrew instala GMP em `/opt/homebrew` (Apple Silicon) ou `/usr/local` (Intel). Ajuste `-I` e `-L` conforme o caso.
- **Windows nativo (MinGW/MSYS2):** use caminhos Windows (ex.: `-IC:/gmp/include -LC:/gmp/lib`).

---

## Uso

```
./naraca -b <bits> [-f] [-p <segundos>]
./naraca -s <arquivo.bin> [-f] [-p <segundos>]
```

### Opções

| Opção | Descrição |
|-------|-----------|
| `-b, --bitsLen <N>` | Gera semente aleatória de N bits |
| `-s, --semente <arq>` | Lê semente de arquivo binário (gerado com `dd`) |
| `-f, --filtro` | Ativa a peneira de terminações (1,3,7,9) |
| `-p, --progress <N.M>` | Mostra progresso a cada N.M segundos |
| `-h, --help` | Exibe ajuda |

### Exemplos

**Gerar semente de 20.000 bits (2500 bytes):**
```bash
dd if=/dev/urandom of=semente.bin bs=2500 count=1
```

**Buscar primo com filtro:**
```bash
./naraca -s semente.bin -f -p 5.0
```

**Buscar primo sem filtro (mesma semente, para comparação honesta):**
```bash
./naraca -s semente.bin -p 5.0
```

**Buscar primo de 10.000 bits com semente aleatória:**
```bash
./naraca -b 10000 -f -p 5.0
```

---

## Resultados Experimentais

### Ambiente de execução

Os experimentos foram conduzidos em:

- **Sistema operacional:** Windows 11
- **Camada de compatibilidade:** Cygwin **3.6.10**
- **Hardware:** Notebook Samsung
- **Processador:** Intel Core i7
- **Memória RAM:** 32 GB
- **Compilador:** GCC (Cygwin)
- **Biblioteca:** GMP 6.3.0 (compilada a partir do código-fonte)

### Comparação honesta (mesma semente, 20.000 bits)

| Métrica | Com filtro | Sem filtro | Diferença |
|---------|-----------|------------|-----------|
| Tempo total | 137,96 s | 139,96 s | **-1,4%** |
| Tentativas totais | 1.784 | 1.784 | iguais |
| Descartados pela peneira | 357 (20,0%) | 0 | — |
| Chamadas ao teste (reps=1) | 1.427 | 1.784 | **-20,0%** |
| Tempo médio por tentativa | 0,07733 s | 0,07845 s | **-1,4%** |
| Últimos 3 dígitos | 237 | 237 | mesmo primo |

### Observações sobre o ambiente

- O **Cygwin 3.6.10** adiciona uma camada de emulação POSIX sobre o Windows, o que pode introduzir overhead em chamadas de sistema (como `fopen`/`fread` de `/dev/urandom`) e em operações de I/O. Os tempos medidos refletem esse ambiente; em Linux nativo, espera-se desempenho ligeiramente superior.
- O **Core i7** e os **32 GB de RAM** são mais do que suficientes para candidatos de 20.000 bits (2,5 KB por `mpz_t`). O gargalo é o teste de primalidade, não a memória.
- A versão do Cygwin pode ser verificada com `cygcheck -V` (neste caso, 3.6.10).
- Os resultados são **específicos desta máquina**. Para comparações entre ambientes, seria necessário rodar o mesmo experimento em Linux nativo e comparar.

### Interpretação

- A peneira **reduz as chamadas ao teste em 20%**, como esperado.
- O ganho de tempo é **modesto (~1,4%)**, porque os candidatos eliminados (terminados em 5) são rejeitados rapidamente pelo teste de primalidade.
- O número de tentativas até achar um primo é **altamente variável** (distribuição geométrica); uma única rodada não é representativa. Para conclusões estatísticas, são necessárias múltiplas sementes.

### Limitações

- O teste `mpz_probab_prime_p` é **probabilístico**, não determinístico. Para prova formal de primalidade, são necessárias ferramentas como ECPP (Primo, CM, FLINT).
- O programa não implementa checkpoint; para execuções longas (100.000+ bits), interrupções perdem o progresso.
- A peneira atual é minimalista; peneiras mais completas (múltiplos de 3, 7, 11, ...) podem trazer ganhos maiores.

---

## Estrutura do Código

```
naraca.c          # Código-fonte principal
Makefile          # Automação de compilação
README.md         # Este documento
```

### Funções principais

| Função | Descrição |
|--------|-----------|
| `gerar_candidato` | Lê `/dev/urandom` e importa como `mpz_t` ímpar |
| `peneira_lista` | Verifica se o último dígito é 1, 3, 7 ou 9 |
| `tempo_atual` | Retorna tempo monotônico em segundos |
| `showHelp` | Exibe ajuda |
| `main` | Parsing de argumentos, loop de busca, relatório |

---

## Makefile

```makefile
# Makefile — Naraca
CC = gcc
CFLAGS = -O3 -march=native -fopenmp -Wall -Wextra
CPPFLAGS = -I/usr/local/include
LDFLAGS = -L/usr/local/lib
LDLIBS = -lgmp
TARGET = naraca
SRC = naraca.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o

run: $(TARGET)
	./$(TARGET) -b 10000 -f -p 5.0

.PHONY: all clean run
```

---

## Referências

### Matemática

- **Teorema de Dirichlet**: Dirichlet, P. G. L. (1837). *Beweis des Satzes, dass jede unbegrenzte arithmetische Progression, deren erstes Glied und Differenz ganze Zahlen ohne gemeinschaftlichen Factor sind, unendlich viele Primzahlen enthält*. Abhandlungen der Königlichen Preußischen Akademie der Wissenschaften zu Berlin.
- **Teoria elementar dos primos**: Hardy, G. H., & Wright, E. M. (2008). *An Introduction to the Theory of Numbers* (6th ed.). Oxford University Press.
- **GMP — GNU Multiple Precision Arithmetic Library**: Documentação oficial em https://gmplib.org/manual/
- **Teste de Baillie-PSW**: Baillie, R., & Wagstaff, S. S. (1980). *Lucas pseudoprimes*. Mathematics of Computation, 35(152), 1391–1417.

### Ferramentas

- **GMP**: https://gmplib.org/
- **Primo (ECPP)**: http://www.ellipsa.eu/public/primo/primo.html
- **CM (ECPP)**: https://www.multiprecision.org/cm/
- **FLINT**: https://flintlib.org/
- **PARI/GP**: https://pari.math.u-bordeaux.fr/
- **GIMPS**: https://www.mersenne.org/

### Artigos e discussões

- Documentação do GMP sobre `mpz_probab_prime_p`: https://gmplib.org/manual/Number-Theoretic-Functions
- Discussões sobre custo de testes de primalidade para números gigantes: lista de e-mails do GMP (arquivos públicos).

---

## Licença

Este projeto é distribuído sob a **MIT License**. Veja o arquivo `LICENSE` para detalhes.

---

## Agradecimentos

- À comunidade GMP, pelo desenvolvimento da biblioteca que torna este experimento possível.
- A Anderson Silva, pela curiosidade e pela condução dos experimentos que motivaram este documento.

---

## Como Contribuir

1. Faça um fork do repositório.
2. Crie uma branch para sua feature (`git checkout -b feature/nova-peneira`).
3. Commit suas mudanças (`git commit -am 'Adiciona peneira de múltiplos de 3'`).
4. Push para a branch (`git push origin feature/nova-peneira`).
5. Abra um Pull Request.

### Ideias para contribuições

- Implementar peneira mais completa (múltiplos de 3, 7, 11, 13, ...).
- Adicionar checkpoint (candidato + contadores) para execuções longas.
- Implementar teste de primalidade determinístico (ECPP) para prova formal.
- Adicionar suporte a múltiplas threads para busca paralela.
- Medir tempo de geração vs teste separadamente.

---

### Dados do Projeto
* *Última atualização: 2026-11-24*
* *Versão: 1.0.0*
* *Autor: Anderson Silva*
* *Co-autores:*
	* IA *Genini* - mais atrapalhou do que ajudou, grande dificuldade de entender o objetvo do projeto
	* IA *DeepSeek* - Ajudou na codificação conforme as minha esoecificações e apresentou sugestões de melhorai de desepenho

	
