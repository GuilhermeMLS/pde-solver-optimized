set terminal png size 640,480;
set output 'gauss_seidel_group_l3.png'
set ylabel 'Memory bandwidth [MBytes/s]'
set xlabel 'Número de pontos em cada dimensão (nx = ny)'
set title 'Resultados do grupo L3 do LiKWID para a função gaussSeidel()'
plot 'gauss_seidel_group_l3' title 'Depois da otimização' with lines, 'trab1_gauss_seidel_group_l3' title 'Antes da otimização' with lines
