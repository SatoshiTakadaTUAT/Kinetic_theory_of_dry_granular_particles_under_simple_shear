set terminal postscript eps enhanced color 32
#set colorsequence classic
set border lw 1
set key maxrows 4
set key right top

set key font "Arial, 28"
set xlabel ''
set ylabel ''
set format x 
set format y 

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
set nologscale y
set yrange[0:1.5]
set ytics 0.5
set key at graph 0.95, 0.30
set output 'eta_phi_SA.eps'
plot \
'e_0.90.txt' u 1:(($3)*dens($1)/etaGD($1, 0.90))        w l lw 4 lt 1 dt 1 title '\',\
'e_0.90_SA.txt' u 1:(-($5+$9)*dens($1)/etaGD($1, 0.90)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(-($5)/etaGD($1, 0.90)):(-($6)/etaGD($1, 0.90)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                                                             pt 6 ps 4 lw 4 lt 1 title '\',\
'e_0.70.txt' u 1:(($3)*dens($1)/etaGD($1, 0.70))       w l lw 4 lt 7 dt 3 title '\',\
'e_0.70_SA.txt' u 1:(-($5+$9)*dens($1)/etaGD($1, 0.70)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(-($5)/etaGD($1, 0.70)):(-($6)/etaGD($1, 0.70)) w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                                                             pt 4 ps 4 lw 4 lt 7 title '\'

set xrange[0:1]
set xtics 0.2
set key at graph 0.85, 0.9
set output 'eta_e_SA.eps'
plot \
'phi_0.40.txt' u 1:(($3)*dens(0.40)/etaGD(0.40, $1))        w l lw 4 lt 1 dt 1 title '\',\
'phi_0.40_SA.txt' u 1:(-($5+$9)*dens(0.40)/etaGD(0.40, $1)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($1)==0.4?  $2 : 1/0):(-($5)/etaGD(0.40, $2)):(-($6)/etaGD(0.40, $2)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                                                             pt 6 ps 4 lw 4 lt 1 title '\',\
'phi_0.10.txt' u 1:(($3)*dens(0.10)/etaGD(0.10, $1))        w l lw 4 lt 7 dt 3 title '\',\
'phi_0.10_SA.txt' u 1:(-($5+$9)*dens(0.10)/etaGD(0.10, $1)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(-($5)/etaGD(0.10, $2)):(-($6)/etaGD(0.10, $2)) w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                                                             pt 4 ps 4 lw 4 lt 7 title '\'

##############################
# theta
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[0:1.5]
set ytics 0.5
set key at graph 0.9, 0.30
set output 'theta_phi_SA.eps'
plot \
'e_0.90.txt'    u 1:(($2)/thetaGD($1, 0.90)) w l lw 4 lt 1 dt 1 title '\',\
'e_0.90_SA.txt' u 1:(($2)/thetaGD($1, 0.90)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):(($3)/thetaGD($1, 0.90)):(($4)/thetaGD($1, 0.90)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                                                               pt 6 ps 4 lw 4 lt 1 title '\',\
'e_0.70.txt'    u 1:(($2)/thetaGD($1, 0.70))       w l lw 4 lt 7 dt 3 title '\',\
'e_0.70_SA.txt' u 1:(($2)/thetaGD($1, 0.70)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):(($3)/thetaGD($1, 0.70)):(($4)/thetaGD($1, 0.70)) w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                                                               pt 4 ps 4 lw 4 lt 7 title '\'

set xrange[0:1]
set xtics 0.2
set yrange[0:2]
set ytics 0.5
set key at graph 0.9, 0.85
set output 'theta_e_SA.eps'
plot \
'phi_0.40.txt'    u 1:(($2)/thetaGD(0.40, $1))        w l lw 4 lt 1 dt 1 title '\',\
'phi_0.40_SA.txt' u 1:(($2)/thetaGD(0.40, $1)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($1)==0.4?  $2 : 1/0):(($3)/thetaGD(0.40, $2)):(($4)/thetaGD(0.40, $2)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                                                               pt 6 ps 4 lw 4 lt 1 title '\',\
'phi_0.10.txt'    u 1:(($2)/thetaGD(0.10, $1))        w l lw 4 lt 7 dt 3 title '\',\
'phi_0.10_SA.txt' u 1:(($2)/thetaGD(0.10, $1)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):(($3)/thetaGD(0.10, $2)):(($4)/thetaGD(0.10, $2)) w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                                                               pt 4 ps 4 lw 4 lt 7 title '\'

##############################
# N1
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[-0.2:1]
set ytics 0.2
set key at graph 0.95, 0.85
set output 'N1_phi_SA.eps'
plot \
'e_0.90.txt' u 1:4                          w l lw 4 lt 1 dt 1 title '\',\
'e_0.90_SA.txt' u 1:(($3+$7-$8)/($2+$6/3.)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):9:10 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                  pt 6 ps 4 lw 4 lt 1 title '\',\
'e_0.70.txt' u 1:4                          w l lw 4 lt 7 dt 3 title '\',\
'e_0.70_SA.txt' u 1:(($3+$7-$8)/($2+$6/3.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):9:10 w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                  pt 4 ps 4 lw 4 lt 7 title '\'

set xrange[0:1]
set xtics 0.2
set yrange[-0.2:2]
set ytics 0.5
set key at graph 0.95, 0.85
set output 'N1_e_SA.eps'
plot \
'phi_0.40.txt' u 1:4                          w l lw 4 lt 1 dt 1 title '\',\
'phi_0.40_SA.txt' u 1:(($3+$7-$8)/($2+$6/3.)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($1)==0.4?  $2 : 1/0):9:10 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                  pt 6 ps 4 lw 4 lt 1 title '\',\
'phi_0.10.txt' u 1:4                          w l lw 4 lt 7 dt 3 title '\',\
'phi_0.10_SA.txt' u 1:(($3+$7-$8)/($2+$6/3.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.10?  $2 : 1/0):9:10 w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                   pt 4 ps 4 lw 4 lt 7 title '\'

##############################
# N2
##############################
set xrange[0:0.5]
set xtics 0.1
set yrange[-0.1:0.3]
set ytics 0.1
set key at graph 0.75, 0.90
set output 'N2_phi_SA.eps'
plot \
'e_0.90.txt' u 1:5                                   w l lw 4 lt 1 dt 1 title '\',\
'e_0.90_SA.txt' u 1:(($4-$3+$7+2.*$8-$6)/($2+$6/3.)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($2)==0.9?  $1 : 1/0):11:12 w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                   pt 6 ps 4 lw 4 lt 1 title '\',\
'e_0.70.txt' u 1:5                                   w l lw 4 lt 7 dt 3 title '\',\
'e_0.70_SA.txt' u 1:(($4-$3+$7+2.*$8-$6)/($2+$6/3.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($2)==0.7?  $1 : 1/0):11:12 w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                   pt 4 ps 4 lw 4 lt 7 title '\'

set xrange[0:1]
set xtics 0.2
set yrange[-0.1:0.5]
set ytics 0.1
set key at graph 1, 0.90
set output 'N2_e_SA.eps'
plot \
'phi_0.40.txt' u 1:5                                   w l lw 4 lt 1 dt 1 title '\',\
'phi_0.40_SA.txt' u 1:(($4-$3+$7+2.*$8-$6)/($2+$6/3.)) w l lw 4 lt 6 dt 2 title '\',\
'coeff.txt' u (($1)==0.4?  $2 : 1/0):11:12   w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0)                                     pt 6 ps 4 lw 4 lt 1 title '\',\
'phi_0.10.txt' u 1:5                                   w l lw 4 lt 7 dt 3 title '\',\
'phi_0.10_SA.txt' u 1:(($4-$3+$7+2.*$8-$6)/($2+$6/3.)) w l lw 4 lt 8 dt 4 title '\',\
'coeff.txt' u (($1)==0.1?  $2 : 1/0):11:12   w yerrorbars pt 4 ps 4 lw 4 lt 7 title '',\
'coeff.txt' u 1:(1/0)                                     pt 4 ps 4 lw 4 lt 7 title '\'