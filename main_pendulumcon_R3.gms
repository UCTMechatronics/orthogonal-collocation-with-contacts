sets n "colocation points" /node1 * node100/
     p "point" /p1*p3/
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

int("p1") = 0.2123;
int("p2") = 0.5905;
int("p3") = 0.9114;

Table  omega(p,p)  Radau matrix
         p1                      p2                      p3
p1       0.19681547722366        -0.06553542585020       0.02377097434822
p2       0.39442431473909       0.29207341166523         -0.04154875212600
p3       0.37640306270047       0.51248582618842         0.11111111111111
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
get_true_time1(n)
get_true_time2(n,p);

get_true_time1(n).. time(n,"p1") =e= h(n)*int("p1");
get_true_time2(n,p)$(ord(p) gt 1).. time(n,p) =e= h(n)*(int(p)-int(p-1));

*RELATIVE ANGLES----------------------------------------------------------------
variables qr(n,p,l);

Equations def_qr(n,p,l);
def_qr(n,p,l)$(ord(n) gt 1 or ord(p) eq card(p)).. qr(n,p,l) =e= 0
+ (q(n,p,"th1"))$sameas(l,"l1")
+ (q(n,p,"th2") - q(n,p,"th1"))$sameas(l,"l2");

*DYNAMICS-----------------------------------------------------------------------
$include "dynamics_RadauEuler.txt";

alias(p,pp)

Equations
interp_q(n,p,i)
interp_dq(n,p,i)
;

interp_q(n,p,i)$(ord(n) gt 1).. q(n,p,i) =e= sum(pp$(ord(pp) eq card(pp)), q(n-1,pp,i)) + h(n)*sum(pp, omega(p,pp)*dq(n,pp,i));
interp_dq(n,p,i)$(ord(n) gt 1).. dq(n,p,i) =e= sum(pp$(ord(pp) eq card(pp)), dq(n-1,pp,i)) + h(n)*sum(pp, omega(p,pp)*ddq(n,"p3",i));

*JOINT STOPS--------------------------------------------------------------------
scalar jointlim;
jointlim = pi/4;

positive variables
bound_up(n,p,sgn)
bound_lo(n,p,sgn)
tau_j_sgn(n,p,sgn)

bound_up_penalty
bound_lo_penalty;

variables
compA_up(n)
compB_up(n)
compA_lo(n)
compB_lo(n)
comp_penalty
;

equations
def_bound_up(n,p)
def_bound_lo(n,p)
def_tau_j(n,p)

def_compA_up(n)
def_compB_up(n)
def_compA_lo(n)
def_compB_lo(n)
def_bound_up_penalty
def_bound_lo_penalty
def_comp_penalty
;

def_bound_up(n,p).. bound_up(n,p,"ps") - bound_up(n,p,"ng") =e= jointlim - qr(n,p,"l2");
def_bound_lo(n,p).. bound_lo(n,p,"ps") - bound_lo(n,p,"ng") =e= -jointlim - qr(n,p,"l2");
def_tau_j(n,p).. tau_j_sgn(n,p,"ps") - tau_j_sgn(n,p,"ng") =e= tau_j(n,p);

def_compA_up(n)$(ord(n) lt card(n)).. compA_up(n) =e= sum(p, bound_up(n+1,p,"ps"));
def_compB_up(n).. compB_up(n) =e= sum(p, tau_j_sgn(n,p,"ng"));
def_compA_lo(n).. compA_lo(n)$(ord(n) lt card(n))  =e= sum(p, bound_lo(n+1,p,"ng"));
def_compB_lo(n).. compB_lo(n) =e= sum(p, tau_j_sgn(n,p,"ps"));
def_bound_up_penalty.. bound_up_penalty =e= sum((n,p), bound_up(n,p,"ng"));
def_bound_lo_penalty.. bound_lo_penalty =e= sum((n,p), bound_lo(n,p,"ps"));
def_comp_penalty.. comp_penalty =e= sum(n, compA_up(n)*compB_up(n) + compA_lo(n)*compB_lo(n));

compA_up.fx(n)$(ord(n) eq card(n)) = 0;
compB_up.fx(n)$(ord(n) eq card(n)) = 0;
compA_lo.fx(n)$(ord(n) eq card(n)) = 0;
compB_lo.fx(n)$(ord(n) eq card(n)) = 0;
bound_up_penalty.fx = 0;
bound_lo_penalty.fx = 0;

*STARTING POINT-----------------------------------------------------------------
q.fx(n,p,"th1")$(ord(n) le 2) = 0;
q.fx(n,p,"th2")$(ord(n) le 2) = 0;

*FINAL POINT--------------------------------------------------------------------
positive variables
final_position(i,sgn);

equation
def_final_position(n,p,i);

def_final_position(n,p,i)$(ord(n) eq card(n) and ord(p) eq card(p)).. final_position(i,"ps") - final_position(i,"ng") =e= q(n,p,i) - pi;

final_position.fx(i,sgn) = 0;
dq.fx(n,p,i)$(ord(n) eq card(n) and ord(p) eq card(p)) = 0;

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

pendulum2.reslim = 600;
pendulum2.workfactor = 10;
pendulum2.threads = 4;
pendulum2.limrow = 20;
pendulum2.optfile = 1;
pendulum2.scaleopt = 1;

option nlp = conopt;

file results /results.csv/ ;
results.nd = 8; results.pw = 10000;

*RANDOMIZE======================================================================
execseed = 1 + gmillisec(jnow);

tprev = 0;
loop (seed,
loop(n,

TT = 2;
h_global = TT/card(n);
h.l(n) = h_global;

q.l(n,p,i) = uniform(-pi, pi);
dq.l(n,p,i) = 0.01;
ddq.l(n,p,i) = 0.01;

tau_c.l(n) = 0.01;
tau_j.l(n,p) = 0.01;
tau_j_sgn.l(n,p,sgn) = 0.01;

Jcost.l = 0.01;
);

tau_c.up(n) = 200;
tau_c.lo(n) = -200;

*SOLVE==========================================================================
comp_penalty.up = inf;
comp_penalty.lo = -inf;
pendulum2.optfile = 2;
solve pendulum2 using nlp minimizing comp_penalty;

pendulum2.optfile = 1;
solve pendulum2 using nlp minimizing Jcost;

tnext = TimeElapsed;
solvetime = tnext - tprev;
tprev = tnext;
Display solvetime;

*===============================================================================
put results;
*put_utilities 'ren' / 'results\Radau3_200_'seed.tl:0'.csv';

put pendulum2.solvestat";" Jcost.l";" comp_penalty.l";" solvetime /;

loop(p$(ord(p) eq card(p)),
put
"0;"
qr.l("node1",p,"l1")";" qr.l("node1",p,"l2")";" tau_c.l("node1")";" tau_j.l("node1",p)";"
dq.l("node1",p,"th1")";" dq.l("node1",p,"th2")";"
/;
);

loop (n$(ord(n) gt 1),
loop(p,
put
time.l(n,p)";"
qr.l(n,p,"l1")";" qr.l(n,p,"l2")";" tau_c.l(n)";" tau_j.l(n,p)";"
dq.l(n,p,"th1")";" dq.l(n,p,"th2")";"
/;
);
);
putclose;
);

