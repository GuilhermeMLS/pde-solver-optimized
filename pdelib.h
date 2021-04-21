/**
 * \file pdelib.h
 * \brief Definições das funções da biblioteca <code>pdelib</code>.
 * \author Guilherme M. Lopes
 */
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/time.h>

#define pi 3.14159265359
#define MAXIMUM_ERROR 6.0E-5
#define calculate_index(i, j, n) (((i) * (n)) + (j))

typedef double t_float;

/**
 * \struct t_LS5Diag
 * \brief Estrutura de dados que representa um sistema linear pentadiagonal.
 */
typedef struct t_LS5Diag {
    double main_diagonal, /**< Constante que representa a diagonal principal */
           bottom_diagonal, /**< Constante que representa a diagonal inferior */
           upper_diagonal, /**< Constante que representa a diagonal superior */
           away_bottom_diagonal, /**< Constante que representa a diagonal inferior afastada */
           away_upper_diagonal; /**< Constante que representa a diagonal superior afastada */
    double *b; /**< Vetor B (termos independetes) */
    unsigned int n; /**< Variável utilizada para calcular o tamanho total da malha */
    unsigned int nx; /**< Número de pontos para discretização na direção das abscissas */
    unsigned int ny; /**< Número de pontos para discretização na direção das ordenadas */
} t_LS5Diag;

void generateOutputFile(
        unsigned int iterations,
        char* output_file,
        double gauss_seidel_total_time,
        double residue_total_time,
        t_float *residues,
        t_LS5Diag *SL,
        t_float *u
);
int get_options(int argc,
    char **argv,
    int *nx,
    int *ny,
    int *max_iterations,
    char **output_file
);

void allocate_and_start_linear_system(t_LS5Diag **SL, int nx, int ny);

void allocate_and_start_solution( t_float **u, t_LS5Diag *SL );

t_float calculate_residues(t_LS5Diag *SL, t_float *u);

void gaussSeidel( t_LS5Diag **SL, t_float **u );


void show_help(char *name);

double timestamp(void);
