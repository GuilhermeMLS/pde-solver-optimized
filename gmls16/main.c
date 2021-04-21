/**
 * \file main.c
 * \author Guilherme M. Lopes
 * \brief Função principal do programa, utiliza a biblioteca <code>pdelib</code> para resolver uma PDE.
 *
 * \mainpage A program to find an approximate solution to a Partial Differential Equation using the Gauss-Seidel iterative method.
 *
 * \section introSec About the program
 * Given a Partial Differential Equation, this program discretizes the domain using the Central Finite Differences method.\n
 * Then, it uses Gauss-Seidel method to find a solution of the generated linear system.
 *
 * **Matéria/Course:** \n
 * Introdução à Computação Científica <i>(Introduction to Scientific Computing)</i> \n
 * **Professor:** \n
 * Daniel Weingaertner \n
 * **Aluno/Student:** \n
 *  Guilherme Morais Lopes dos Santos - GRR20163043
 *
 * \date 24 out 2019
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdelib.h"

// Autovectorization
#pragma GCC optimize("O3","unroll-loops","omit-frame-pointer","inline") //Optimization flags
#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX
#pragma GCC target("avx")  //Enable AVX
#include <x86intrin.h> //AVX/SSE Extensions

//#define DEBUG

// This block enables to compile the code with and without the likwid header in place
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

int main(int argc, char *argv[]) {
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    int nx; // M
    int ny; // N
    int max_iterations;
    char *output_file;

    if (argc < 2) {
        show_help(argv[0]);
    }

    get_options(argc, argv, &nx, &ny, &max_iterations, &output_file);


    t_LS5Diag *SL; // Sistema Linear Pentadiagonal
    allocate_and_start_linear_system(&SL, nx, ny);

    t_float *u; // Solução
    allocate_and_start_solution(&u, SL);

    // Residuos
    t_float *residues = malloc(sizeof(t_float)*max_iterations);

    double initial_time;
    double gauss_seidel_total_time = 0;
    double residue_total_time = 0;

    unsigned int k = 0;
    do {
        k++;

        initial_time = timestamp();
        LIKWID_MARKER_START("gaussSeidel");
        gaussSeidel(&SL, &u);
        LIKWID_MARKER_STOP("gaussSeidel");
        gauss_seidel_total_time += (timestamp() - initial_time);

        // Calcula Residuo da iteração
        initial_time = timestamp();
	LIKWID_MARKER_START("calculate_residues");
        residues[k-1] = calculate_residues(SL, u);
	LIKWID_MARKER_STOP("calculate_residues");
        residue_total_time += (timestamp() - initial_time);

    } while(((k <= max_iterations) && (residues[k-1] >= MAXIMUM_ERROR)));


    generateOutputFile(
        k,
        output_file,
        gauss_seidel_total_time,
        residue_total_time,
        residues,
        SL,
        u
    );
    LIKWID_MARKER_CLOSE;
    return 0;
}
