set terminal png size 640,480;
set output 'residues_group_l2cache.png'
set ylabel 'Taxa de miss (%)'
set xlabel 'Número de pontos em cada dimensão (nx = ny)'
set title 'Resultados do grupo L2CACHE do LiKWID para a função calculate-residues()'
plot 'residues_group_l2cache' title 'Depois da otimização' with lines, 'trab1_residues_group_l2cache' title 'Antes da otimização' with lines
