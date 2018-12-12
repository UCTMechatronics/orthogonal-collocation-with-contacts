sets n "colocation points" /node1 * node2000/
     p "point" /p1/
     i "coordinates" /th1,th2/
     l "links" /l1*l4/
     sgn "sign" /ps,ng/

     seed /s1*s1/

Scalars
g /9.81/
pi /3.14159265359/
;

*Quadrature parameters----------------------------------------------------------
parameters
int(p) "time interval of point within element";

int("p1") = 1;

Table  omega(p,p)  Radau matrix
         p1
p1       1
;

Parameters
m(l) "mass" /l1 1, l2 1/
len(l) "length" /l1 1, l2 1/
clen(l) "length to COM"
In(l) "intertia";

clen(l) = 0.5*len(l);
In(l) = (m(l)*len(l)**3)/12;

variables
q(n,p,i)
dq(n,p,i)
ddq(n,p,i)

tau_c(n)
tau_j(n,p);

*TIME---------------------------------------------------------------------------
scalar
TT
h_global;

positive variables
h(n)
;

Equations
time_min(n)
time_max(n)
;

time_min(n).. h(n) =g= 0.8*h_global;
time_max(n).. h(n) =l= 1.2*h_global;

variables
time(n,p);

Equations
get_true_time1(n);

get_true_time1(n).. time(n,"p1") =e= h(n)*int("p1");

*RELATIVE ANGLES----------------------------------------------------------------
variables qr(n,p,l);

Equations def_qr(n,p,l);
def_qr(n,p,l)$(ord(n) gt 1 or ord(p) eq card(p)).. qr(n,p,l) =e= 0
+ (q(n,p,"th1"))$sameas(l,"l1")
+ (q(n,p,"th2") - q(n,p,"th1"))$sameas(l,"l2");

*DYNAMICS-----------------------------------------------------------------------
$include "dynamics.txt";

alias(p,pp)

Equations
interp_q(n,p,i)
interp_dq(n,p,i)
;

interp_q(n,p,i)$(ord(n) gt 1).. q(n,p,i) =e= sum(pp$(ord(pp) eq card(pp)), q(n-1,pp,i)) + h(n)*sum(pp, omega(p,pp)*dq(n,pp,i));
interp_dq(n,p,i)$(ord(n) gt 1).. dq(n,p,i) =e= sum(pp$(ord(pp) eq card(pp)), dq(n-1,pp,i)) + h(n)*sum(pp, omega(p,pp)*ddq(n,pp,i));

*JOINT STOPS--------------------------------------------------------------------
scalar jointlim;
jointlim = pi/8;

positive variables
bound_up(n,p,sgn)
bound_lo(n,p,sgn)
tau_j_sgn(n,p,sgn)

bound_up_penalty
bound_lo_penalty;

variables
comp_penalty
;

equations
def_bound_up(n,p)
def_bound_lo(n,p)
def_tau_j(n,p)

def_bound_up_penalty
def_bound_lo_penalty
def_comp_penalty
;

def_bound_up(n,p).. bound_up(n,p,"ps") - bound_up(n,p,"ng") =e= jointlim - qr(n,p,"l2");
def_bound_lo(n,p).. bound_lo(n,p,"ps") - bound_lo(n,p,"ng") =e= -jointlim - qr(n,p,"l2");
def_tau_j(n,p).. tau_j_sgn(n,p,"ps") - tau_j_sgn(n,p,"ng") =e= tau_j(n,p);

def_bound_up_penalty.. bound_up_penalty =e= sum((n,p), bound_up(n,p,"ng"));
def_bound_lo_penalty.. bound_lo_penalty =e= sum((n,p), bound_lo(n,p,"ps"));
def_comp_penalty.. comp_penalty =e= sum((n,p), bound_up(n,p,"ps")*tau_j_sgn(n,p,"ng") + bound_lo(n,p,"ng")*tau_j_sgn(n,p,"ps"));

bound_up_penalty.fx = 0;
bound_lo_penalty.fx = 0;

*STARTING POINT-----------------------------------------------------------------
q.fx(n,p,"th1")$(ord(n) lt 2) = pi/2;
q.fx(n,p,"th2")$(ord(n) lt 2) = pi/2;
dq.fx(n,p,i)$(ord(n) lt 1) = 0; 
tau_c.fx(n) = 0;

*FINAL POINT--------------------------------------------------------------------
$ontext
positive variables
final_position(i,sgn);

equation
def_final_position(n,p,i);

def_final_position(n,p,i)$(ord(n) eq card(n) and ord(p) eq card(p)).. final_position(i,"ps") - final_position(i,"ng") =e= q(n,p,i) - pi;

final_position.fx(i,sgn) = 0;
dq.fx(n,p,i)$(ord(n) eq card(n) and ord(p) eq card(p)) = 0;
$offtext

*SOLVE--------------------------------------------------------------------------
scalar
tnext
tprev
solvetime;

variable
Jcost;

equation
cost;

cost.. Jcost =e= 0
+comp_penalty
+sum(n, h(n)*tau_c(n)*tau_c(n));
;

model pendulum2 /all/;

*dq.scale(n,p,i) = 10;
*ddq.scale(n,p,i) = 100;
*h.scale(n) = 0.1;

pendulum2.reslim = 20000;
pendulum2.workfactor = 10;
pendulum2.threads = 4;
pendulum2.limrow = 20;
pendulum2.optfile = 1;
pendulum2.scaleopt = 1;

option nlp = conopt;

file results /euler_ref_90_2000.csv/ ;
results.nd = 8; results.pw = 10000;

*RANDOMIZE======================================================================
execseed = 1 + gmillisec(jnow);

tprev = 0;
loop (seed,
loop(n,

TT = 2;
h_global = TT/card(n);
h.l(n) = h_global;

*q.l(n,p,i) = uniform(-pi, pi);
dq.l(n,p,i) = 0.01;
ddq.l(n,p,i) = 0.01;

tau_c.l(n) = 0.01;
tau_j.l(n,p) = 0.01;
tau_j_sgn.l(n,p,sgn) = 0.01;

Jcost.l = 0.01;
);

*SOLVE==========================================================================
pendulum2.optfile = 2;
solve pendulum2 using nlp minimizing comp_penalty;

pendulum2.optfile = 1;
comp_penalty.fx = 0;
solve pendulum2 using nlp minimizing Jcost;

tnext = TimeElapsed;
solvetime = tnext - tprev;
tprev = tnext;

*===============================================================================
put results;
*put_utilities 'ren' / 'results\Euler_'seed.tl:0'.csv';

put pendulum2.solvestat";" Jcost.l";" comp_penalty.l";" solvetime /;

loop(p$(ord(p) eq card(p)),
put
"0;"
qr.l("node1",p,"l1")";" qr.l("node1",p,"l2")";" tau_c.l("node1")";" tau_j.l("node1",p)";"
/;
);

loop (n$(ord(n) gt 1),
loop(p,
put
time.l(n,p)";"
qr.l(n,p,"l1")";" qr.l(n,p,"l2")";" tau_c.l(n)";" tau_j.l(n,p)";"
/;
);
);
putclose;
);

