
set terminal svg background "#ffffff" 
set grid lw 2
set border lw 2
set xlabel "X" font "Arial,16"
set ylabel "Y" font "Arial,16"
set key font "Arial,16"

$well_defined_example_data << EOD
-1.0 1.0   
0.0 0.0
1.0 1.0
EOD

$first_and_overdetermined_example_data << EOD
0.0 1.0   
1.0 2.1
2.0 2.9
3.0 4.2
EOD

$under_determined_example_data << EOD
0.0 1.0   
1.0 2.0
EOD

$rank_deficient_example_data << EOD
0.0 1.0
0.0 1.25   
1.0 2.0
EOD

fit_example(x) = 1.0399999999999996 + 0.890000000000001 * x + 0.04999999999999968 * x**2
well_defined_fit(x) = x**2
overdetermined_fit(x) = 0.99 + 1.04 * x
under_determined_fit(x) = 0.9999999999999999 + 0.33333333333333326 * x + 0.3333333333333332 * x**2 + 0.3333333333333335 * x**3
under_determined_candidate_1(x) = 1.0 + x + (1 - x) * x**2
under_determined_candidate_2(x) = 1.0 + x**3
rank_deficient_fit(x) = 1.1249999999999998 + 0.43750000000000017 * x + 0.43750000000000006 * x**2

set output "./polynomial_fit/original_scatter.svg"
set key right bottom
plot $first_and_overdetermined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     #fit_example(x) w l lw 4 lc 14 title "Fit by polynomial degree 2", \
     #overdetermined_fit(x) w l lw 4 lc 15 title "Fit by polynomial degree 1"

set output "./polynomial_fit/degree1_vs_degree2.svg"
set key right bottom
plot $first_and_overdetermined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     fit_example(x) w l lw 4 lc -1 lt -1 title "Fit by polynomial degree 2", \
     overdetermined_fit(x) w l lw 4 lc -1 dt (5,3) title "Fit by polynomial degree 1"

set output "./polynomial_fit/fit_example.svg"
set key right bottom
plot $first_and_overdetermined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     fit_example(x) w l lw 4 lc -1 title "Fit by polynomial degree 2"

set output "./polynomial_fit/well_defined.svg"
set key center top
plot $well_defined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
      well_defined_fit(x) w l lw 4 lc -1 title "Fit: p(x) =  x^2"

set output "./polynomial_fit/over_determined.svg"
set key right bottom
plot $first_and_overdetermined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     overdetermined_fit(x) w l lw 4 lc -1 title "Fit: p(x) =  0.99 + 1.04x"

set output "./polynomial_fit/under_determined.svg"
set key right bottom  Left width -12
plot $under_determined_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     under_determined_candidate_1(x) w l lw 2 lc -1 dt (5,3) title "Candidate 1: p(x) =  1 + x + (1-x)x^2",\
     under_determined_candidate_2(x) w l lw 2 lc -1 dt (3,5)title "Candidate 2: p(x) =  1 + x^3",\
     under_determined_fit(x) w l lw 4 lc -1 title "Fit: p(x) =  1 + (x+x^2+x^3)/3"

set output "./polynomial_fit/rank_deficient.svg"
set key right bottom
plot $rank_deficient_example_data using 1:2 w p pt 7 ps 1.5 lc -1 title "Data",\
     rank_deficient_fit(x) w l lw 4 lc -1 title "Fit: p(x) =  1.125 + 0.4375x + 0.4375x^2"

reset