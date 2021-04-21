# optimized Partial Differential Equation Solver
An optimized version of my [pdeSolver](https://github.com/GuilhermeMLS/pdeSolver) with several performance gains provided by optimization techniques such as _loop unrolling_ and C programming best practices.

## What is pdeSolver
[pdeSolver](https://github.com/GuilhermeMLS/pdeSolver) is a program to find an approximate solution to a Partial Differential Equation (PDE) using the Gauss-Seidel iterative method to solve the pentadiagonal matrix that represents the Linear System found with the aproximation of the partial derivatives of the equation.

## Some of the performance results

 - The graphics were plotted with [Gnuplot](http://www.gnuplot.info)
 - The indicators were mesured using [Likwid](https://github.com/RRZE-HPC/likwid)
 - Details about the computer were the data were collected can be found at `reports_v2/Artigo.pdf/`
 - Other performance results can be found at `/reports_v2/gaussSeidel` and `/reports_v2/residues`

![MFLOP/s](/reports_v2/gaussSeidel/gauss_seidel_group_flopsdp.png)

MFLOP/s after/before optimization

![CACHEL2/s](/reports_v2/gaussSeidel/gauss_seidel_group_l2cache.png)

Cache L2 usage after/before optimization




