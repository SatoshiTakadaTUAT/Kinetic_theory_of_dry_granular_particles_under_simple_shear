set terminal postscript eps enhanced color 32
#set colorsequence classic
set border lw 1

set key font"Arial, 28"
set key maxrows 4
set yrange[0:2.5]
set ytics 0.5
set key at graph 0.7, 0.90

set nologscale x
set nologscale y
set key spacing 1.3
set xrange[0:1]
set yrange[0:1.5]
set xlabel ''
set ylabel ''
set format x 
set format y 
set xtics 0.2
set ytics 0.5


dens(x)=(6./pi)*x
set key maxrows 4

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

set output 'eta_theta_e_N.eps'
plot \
'phi_0.30.txt' u 1:(($3)*dens(0.30)/etaGD(0.30, $1)) w l lw 4 lt 1 dt 1 title '\',\
'coeff_4000.txt' u (($1)==0.3?  $2 : 1/0):(-($5)/etaGD(0.30, $2)):(-($6)/etaGD(0.30, $2)) w yerrorbars pt 6 ps 4 lw 4 lt 1 title '',\
'coeff.txt'      u (($1)==0.3?  $2 : 1/0):(-($5)/etaGD(0.30, $2)):(-($6)/etaGD(0.30, $2)) w yerrorbars pt 4 ps 4 lw 4 lt 1 title '',\
'coeff.txt' u 1:(1/0) pt 6 ps 4 lw 4 lt 1 title '4000',\
'coeff.txt' u 1:(1/0) pt 4 ps 4 lw 4 lt 1 title '1000',\
'phi_0.30.txt' u 1:(($2)/thetaGD(0.30, $1)) w l lw 4 lt 6 dt 2 title '\',\
'coeff_4000.txt' u (($1)==0.3?  $2 : 1/0):(($3)/thetaGD(0.30, $2)):(($4)/thetaGD(0.30, $2)) w yerrorbars pt 8 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u (($1)==0.3?  $2 : 1/0):(($3)/thetaGD(0.30, $2)):(($4)/thetaGD(0.30, $2)) w yerrorbars pt 2 ps 4 lw 4 lt 6 title '',\
'coeff.txt' u 1:(1/0) pt 8 ps 4 lw 4 lt 6 title '\',\
'coeff.txt' u 1:(1/0) pt 2 ps 4 lw 4 lt 6 title '\'
