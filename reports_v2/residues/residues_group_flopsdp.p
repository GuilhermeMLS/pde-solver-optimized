set terminal png size 640,480;
set output 'residues_group_flopsdp.png'
set ylabel 'Milhões de operações em ponto flutuante [MFLOP/s]'
set xlabel 'Número de pontos em cada dimensão (nx = ny)'
set title 'Resultados do grupo FLOPS-DP do LiKWID para a função calculate-residues()'
plot 'residues_group_flopsdp' title 'Depois da otimização' with lines, 'trab1_residues_group_flopsdp' title 'Antes da otimização' with lines
