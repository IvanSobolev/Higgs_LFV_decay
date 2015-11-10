set grid
unset key
set logscale y
set format y "10^{%L}"
set xlabel "ms"
set ylabel "Br"
set terminal jpeg enhanced size 900,500
set output 'Y_ms_8F_neg.jpeg'
plot "LFV_output" u 1:2 w p pt 1 ps 1,1.51e-02,0.078e-02
