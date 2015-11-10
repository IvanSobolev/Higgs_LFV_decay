set grid
unset key
set logscale x
set logscale y
set format y "10^{%L}"
set format x "10^{%L}"
set xlabel "Br(tau->mu+gm)"
set ylabel "Br(h-> tau+mu)"
set terminal jpeg enhanced size 900,500
set output 'br_ms_8F_neg.jpeg'
plot "LFV_output" u 11:2 w p pt 1 ps 1,1.51e-02,0.078e-02