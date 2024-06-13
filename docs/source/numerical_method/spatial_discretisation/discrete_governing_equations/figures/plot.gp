# physical coordinate
{
  reset
  xmin = - 0.5
  xmax = + 4.0
  ymin = - 0.5
  ymax = + 4.0
  lx = xmax-xmin
  ly = ymax-ymin
  set terminal epslatex standalone color size lx,ly font ',12'
  set output 'car.tex'
  unset border
  set lmargin 0
  set rmargin 0
  set bmargin 0
  set tmargin 0
  unset xlabel
  unset ylabel
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  unset xtics
  unset ytics
  set size ratio -1
  set samples 10000
  set style line 1 lc rgb '#000000' lw 10
  num = 8
  xmin = 0.
  xmax = 3.5
  ymin = 0.
  ymax = 3.5
  deltax = (xmax - xmin) / num
  deltay = (ymax - ymin) / num
  do for [i = 0 : num : 1] {
    if (0 == i) {
      x = xmin
    } else if (num == i) {
      x = xmax
    } else {
      x = i - num / 2
      x = xmin + (xmax - xmin) * 1. / (1. + exp(- 1. * x))
    }
    set arrow from first x, first ymin to first x, first ymax nohead ls 1
    y = ymin + i * deltay
    set arrow from first xmin, first y to first xmax, first y nohead ls 1
  }
  plot \
    NaN notitle
}

# computational coordinate
{
  reset
  xmin = - 0.5
  xmax = + 4.0
  ymin = - 0.5
  ymax = + 4.0
  lx = xmax-xmin
  ly = ymax-ymin
  set terminal epslatex standalone color size lx,ly font ',12'
  set output 'comp.tex'
  unset border
  set lmargin 0
  set rmargin 0
  set bmargin 0
  set tmargin 0
  unset xlabel
  unset ylabel
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  unset xtics
  unset ytics
  set size ratio -1
  set samples 10000
  set style line 1 lc rgb '#000000' lw 10
  num = 8
  xmin = 0.
  xmax = 3.5
  ymin = 0.
  ymax = 3.5
  deltax = (xmax - xmin) / num
  deltay = (ymax - ymin) / num
  do for [i = 0 : num : 1] {
    x = xmin + i * deltax
    set arrow from first x, first ymin to first x, first ymax nohead ls 1
    y = ymin + i * deltay
    set arrow from first xmin, first y to first xmax, first y nohead ls 1
  }
  plot \
    NaN notitle
}

# merge above elements
{
  reset
  xmin = 0.0
  xmax = 9.0
  ymin = 0.0
  ymax = 5.25
  lx = xmax - xmin
  ly = ymax - ymin
  set terminal epslatex standalone color size lx,ly font ',20'
  set output 'result.tex'
  unset border
  set lmargin 0
  set rmargin 0
  set bmargin 0
  set tmargin 0
  unset xlabel
  unset ylabel
  set xrange [xmin:xmax]
  set yrange [ymin:ymax]
  unset xtics
  unset ytics
  set size ratio -1
  ref = 4.5
  set label 'Physical      coordinate' center at first 0.5 * ref, first 1. * ref + 0.35 textcolor rgb '#000000'
  set label 'Computational coordinate' center at first 1.5 * ref, first 1. * ref + 0.35 textcolor rgb '#000000'
  set label '$\left( x, y \right)$'         center at first 0.5 * ref, first 1. * ref - 0.05 textcolor rgb '#000000'
  set label '$\left( \xi^1, \xi^2 \right)$' center at first 1.5 * ref, first 1. * ref - 0.05 textcolor rgb '#000000'
  set label '\includegraphics[width=4.500in, height=4.500in]{car.pdf}'  center at first 0.5 * ref, first 0.5 * ref
  set label '\includegraphics[width=4.500in, height=4.500in]{comp.pdf}' center at first 1.5 * ref, first 0.5 * ref
  plot \
    NaN notitle
}

