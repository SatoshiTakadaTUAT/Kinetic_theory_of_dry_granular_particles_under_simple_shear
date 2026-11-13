set terminal postscript eps enhanced color 32
#set colorsequence classic
set border lw 1
set key maxrows 4
set key right top

set key font "Arial, 28"
set xlabel ''
set ylabel ''
set format x 

##############################
g0CS(phi) = (1.-phi/2.)/(1.-phi)**3
g0T(phi)  = g0CS(0.49)*(0.639-0.49)/(0.639-phi)
g0(phi) = (phi < 0.49)? g0CS(phi) : g0T(phi)
H(phi, e0) = phi*(1. + 2.*(1.+e0)*phi*g0(phi))
etak(phi, e0) = (1. - 0.4*(1.+e0)*(1.-3.*e0)*phi*g0(phi)) / ((1. - 0.25*(1.-e0)**2 - (5./24.)*(1.-e0**2))*g0(phi))
etac(phi, e0) = 0.8*(1.+e0)*phi*g0(phi)*etak(phi, e0)
etab(phi, e0) = (384./(25.*pi))*(1.+e0)*phi**2*g0(phi)
eta(phi, e0) = etak(phi, e0) + etac(phi, e0) + etab(phi, e0)
J(phi, e0) = (5.*sqrt(pi)/96.)*eta(phi, e0)
K(phi, e0) = (12./sqrt(pi))*phi**2*g0(phi)*(1.-e0**2)
etaGD(phi, e0) = sqrt(J(phi, e0)**3/K(phi, e0))*(6./pi)
thetaGD(phi, e0) = J(phi, e0)/K(phi, e0)

dens(x)=(6./pi)*x
##############################

##############################
# eta
##############################
set xrange[0:0.5]
set xtics 0.1
set logscale y
set yrange[1e-1:1e3]
set format y '10^{%L}'
set ytics 10
set key at graph 0.95, 0.90
set output 'eta_phi.eps'
plot \
'e_0.95.txt' u 1:3 w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:3 w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:3 w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:3 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):(-($5)/dens($1)):(-($6)/dens($1))  w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(-($5)/dens($1)):(-($6)/dens($1))  w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):(-($5)/dens($1)):(-($6)/dens($1))  w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(-($5)/dens($1)):(-($6)/dens($1))  w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[1e-1:1e5]
set key at graph 0.6, 0.90
set output 'eta_e.eps'
plot \
'phi_0.45.txt' u 1:3 w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:3 w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:3 w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:3 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):(-($5)/dens($1)):(-($6)/dens($1)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(-($5)/dens($1)):(-($6)/dens($1)) w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(-($5)/dens($1)):(-($6)/dens($1)) w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):(-($5)/dens($1)):(-($6)/dens($1)) w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# theta
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[1e-1:1e4]
set key at graph 0.95, 0.90
set output 'theta_phi.eps'
plot \
'e_0.95.txt' u 1:2 w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:2 w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:2 w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:2 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):3:4 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9? $1 : 1/0):3:4  w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8? $1 : 1/0):3:4  w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7? $1 : 1/0):3:4  w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[1e-2:1e6]
set ytics 100
set key at graph 0.6, 0.90
set output 'theta_e.eps'
plot \
'phi_0.45.txt' u 1:2 w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:2 w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:2 w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:2 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):3:4 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):3:4 w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):3:4 w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):3:4 w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# eta ratio
##############################
set nologscale y
set format y 

set xrange[0:0.5]
set xtics 0.1
set yrange[0:2.5]
set ytics 0.5
set key at graph 0.6, 0.90

set output 'eta_ratio_phi.eps'
plot \
'e_0.95.txt' u 1:(($3)*dens($1)/etaGD($1, 0.95)) w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:(($3)*dens($1)/etaGD($1, 0.90)) w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:(($3)*dens($1)/etaGD($1, 0.80)) w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:(($3)*dens($1)/etaGD($1, 0.70)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):(-($5)/etaGD($1, 0.95)):(-($6)/etaGD($1, 0.95)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(-($5)/etaGD($1, 0.90)):(-($6)/etaGD($1, 0.90)) w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):(-($5)/etaGD($1, 0.80)):(-($6)/etaGD($1, 0.80)) w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(-($5)/etaGD($1, 0.70)):(-($6)/etaGD($1, 0.70)) w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set output 'eta_ratio_e.eps'
plot \
'phi_0.45.txt' u 1:(($3)*dens(0.45)/etaGD(0.45, $1)) w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:(($3)*dens(0.30)/etaGD(0.30, $1)) w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:(($3)*dens(0.10)/etaGD(0.10, $1)) w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:(($3)*dens(0.01)/etaGD(0.01, $1)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):(-($5)/etaGD(0.45, $2)):(-($6)/etaGD(0.45, $2)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(-($5)/etaGD(0.30, $2)):(-($6)/etaGD(0.30, $2)) w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(-($5)/etaGD(0.10, $2)):(-($6)/etaGD(0.10, $2)) w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):(-($5)/etaGD(0.01, $2)):(-($6)/etaGD(0.01, $2)) w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# theta ratio
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[0:2.5]
set ytics 0.5
set output 'theta_ratio_phi.eps'
plot \
'e_0.95.txt' u 1:(($2)/thetaGD($1, 0.95)) w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:(($2)/thetaGD($1, 0.90)) w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:(($2)/thetaGD($1, 0.80)) w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:(($2)/thetaGD($1, 0.70)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):(($3)/thetaGD($1, 0.95)):(($4)/thetaGD($1, 0.95)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(($3)/thetaGD($1, 0.90)):(($4)/thetaGD($1, 0.90)) w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):(($3)/thetaGD($1, 0.80)):(($4)/thetaGD($1, 0.80)) w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(($3)/thetaGD($1, 0.70)):(($4)/thetaGD($1, 0.70)) w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[0:2]
set output 'theta_ratio_e.eps'
plot \
'phi_0.45.txt' u 1:(($2)/thetaGD(0.45, $1)) w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:(($2)/thetaGD(0.30, $1)) w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:(($2)/thetaGD(0.10, $1)) w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:(($2)/thetaGD(0.01, $1)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):(($3)/thetaGD(0.45, $2)):(($4)/thetaGD(0.45, $2)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(($3)/thetaGD(0.30, $2)):(($4)/thetaGD(0.30, $2)) w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(($3)/thetaGD(0.10, $2)):(($4)/thetaGD(0.10, $2)) w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):(($3)/thetaGD(0.01, $2)):(($4)/thetaGD(0.01, $2)) w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# N1
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[-0.2:1]
set ytics 0.2
set key at graph 0.95, 0.90
set output 'N1_phi.eps'
plot \
'e_0.95.txt' u 1:4 w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:4 w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:4 w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:4 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):9:10 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):9:10 w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):9:10 w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):9:10 w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[-0.5:2.5]
set ytics 0.5
set output 'N1_e.eps'
plot \
'phi_0.45.txt' u 1:4 w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:4 w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:4 w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:4 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):9:10 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):9:10 w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):9:10 w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):9:10 w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# N2
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[-0.1:0.3]
set ytics 0.1
set key at graph 0.60, 0.90
set output 'N2_phi.eps'
plot \
'e_0.95.txt' u 1:5 w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:5 w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:5 w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:5 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):11:12 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):11:12 w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):11:12 w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):11:12 w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[-0.2:0.6]
set ytics 0.2
set key at graph 0.95, 0.90
set output 'N2_e.eps'
plot \
'phi_0.45.txt' u 1:5 w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:5 w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:5 w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:5 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):11:12 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):11:12 w yerrorbars pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):11:12 w yerrorbars pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):11:12 w yerrorbars pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# mu
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[0:1]
set ytics 0.2
set key at graph 0.95, 0.90
set output 'mu_phi.eps'
plot \
'e_0.95.txt' u 1:6 w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:6 w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:6 w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:6 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):(-$5/$7) pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(-$5/$7) pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):(-$5/$7) pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(-$5/$7) pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[0:1.2]
set output 'mu_e.eps'
plot \
'phi_0.45.txt' u 1:6 w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:6 w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:6 w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:6 w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):(-$5/$7) pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(-$5/$7) pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(-$5/$7) pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):(-$5/$7) pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'

##############################
# angle
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[-0.6:0.2]
set ytics 0.2
set key at graph 0.95, 0.3
set output 'angle_phi.eps'
plot \
'e_0.95.txt' u 1:($7-(pi/4.)) w l lw 4 lt 1 dt 1 title '\',\
'e_0.90.txt' u 1:($7-(pi/4.)) w l lw 4 lt 6 dt 2 title '\',\
'e_0.80.txt' u 1:($7-(pi/4.)) w l lw 4 lt 7 dt 3 title '\',\
'e_0.70.txt' u 1:($7-(pi/4.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.95? $1 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($2)==0.8?  $1 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.95',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.90',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.80',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.70'

set xrange[0:1]
set xtics 0.2
set yrange[-1.2:0.2]
set output 'angle_e.eps'
plot \
'phi_0.45.txt' u 1:($7-(pi/4.)) w l lw 4 lt 1 dt 1 title '\',\
'phi_0.30.txt' u 1:($7-(pi/4.)) w l lw 4 lt 6 dt 2 title '\',\
'phi_0.10.txt' u 1:($7-(pi/4.)) w l lw 4 lt 7 dt 3 title '\',\
'phi_0.01.txt' u 1:($7-(pi/4.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.45? $2 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 4 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 8 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u (($1)==0.01? $2 : 1/0):(-0.5*atan2(2*($5),($9))-(pi/4.)) pt 2 ps 4 lw 4 lt 8 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '0.45',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 6 title '0.30',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 7 title '0.10',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 8 title '0.01'