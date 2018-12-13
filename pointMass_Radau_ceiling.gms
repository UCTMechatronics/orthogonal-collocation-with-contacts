$Title Point Mass hitting the ceiling!
file test / test.txt/;
test.pw=1000;
test.nd=8;
put test;

Sets    i     number of finite elements             /1*50/
        j     number of internal collocation points /1*3/;

Alias (j,k);

*$ontext
Table  a(j,j)  First order derivatives collocation matrix

                1                  2                  3
      1  0.19681547722366   0.39442431473909   0.37640306270047
      2 -0.06553542585020   0.29207341166523   0.51248582618842
      3  0.02377097434822  -0.04154875212600   0.11111111111111;
*$offtext

$ontext
parameter a(j,j);
a('1','1') =0.416666125000187;
a('1','2') =0.749999625000187;
a('2','1') =-0.083333125000187;
a('2','2') =0.250000374999812;
$offtext

$ontext
parameter a(j,j);
a('1','1') =0.0729985470983077;
a('1','2') =0.153774957379098;
a('1','3') =0.140062829563678;
a('1','4') =0.144894106053214;
a('1','5') =0.143713354932705;

a('2','1') = -0.0267351488651998;
a('2','2') = 0.146215310176912;
a('2','3') = 0.298967347667733;
a('2','4') = 0.276500178980935;
a('2','5') = 0.281356148864136;

a('3','1') = 0.0186767991321094;
a('3','2') = -0.0364448253412411;
a('3','3') = 0.167584366541297;
a('3','4') = 0.325797505377646;
a('3','5') = 0.311826067970553;

a('4','1') = -0.0128789879185484;
a('4','2') = 0.021233179173723;
a('4','3') = -0.0339686756088285;
a('4','4') =  0.128757291308154;
a('4','5') =  0.223104543357278;

a('5','1') = 0.00504279055333108;
a('5','2') = -0.0079356213884913;
a('5','3') = 0.0109441318361206;
a('5','4') =  -0.0157090817199498;
a('5','5') =  0.0399998848753271;
$offtext

*parameter a(j,j);
*        a('1','1') = 1.0;

parameter p(j);

p('1') = 0.155051;
p('2') = 0.644949;
p('3') = 1;
*p('1') = 1;


*p('1') = 0.333333;
*p('2') = 1;

$ontext
p('1') = 0.057104;
p('2') = 0.276843;
p('3') = 0.583590;
p('4') =  0.860240;
p('5') = 1;
$offtext

Parameter
         time /1/
         m /10/
         g /9.81/
         eps1 /1E3/
         nfe
         ncp
        ;

nfe = card(i);
ncp = card(j);

Variables y(i,j)       y
          eta(i)
          dy(i,j)      y velocity
          ddy(i,j)     y acceleration


          tt(i,j)        time
          y_0(i)

          dy_0(i)

          tt_0(i)
          phi1              objective function 1
          h(i)
;

positive variables
         penTot
         pen
         eta(i)
*         eps1(i)
         a_contact_1(i,j)
         a_contact_1_0(i)
         a_contact_1_pr(i)

         b_contact_1(i,j)
         b_contact_1_0(i)
         b_contact_1_pr(i)

         yPos_0(i)
         yPos(i,j)
         yNeg_0(i)
         yNeg(i,j)

         yPosSlack(i)

         lambda1(i)
         lambda1_0
;

 Equations  fobj1         feasible objective definition
            FECOL_y(i,j)
            FECOL_dy(i,j)
            FECOL_tt(i,j)
            CON_y(i)
            CON_dy(i)
            CON_tt(i)
            ;
* OBJECTIVE FUNCTION
fobj1..  phi1 =e= 1e3*pen + 1e1*sum(i,a_contact_1_pr(i)*lambda1(i)*sqr(h(i)-1e-5*time/card(i) ) );

* COLLOCATION CONSTRAINTS
FECOL_y(i,j)$(ord(i) le nfe).. y(i,j) =e= y_0(i) + h(i)*sum(k,a(k,j)*dy(i,k));
FECOL_dy(i,j)$(ord(i) le nfe).. dy(i,j) =e= dy_0(i) + h(i)*sum(k,a(k,j)*ddy(i,k));
FECOL_tt(i,j)$(ord(i) le nfe).. tt(i,j) =e= tt_0(i) + h(i)*sum(k,a(k,j)) ;

* CONTINUITY CONSTRAINTS
CON_y(i)$(ord(i) gt 1 and ord(i) le nfe).. y_0(i) =e= y(i-1,'3');
CON_dy(i)$(ord(i) gt 1 and ord(i) le nfe).. dy_0(i) =e= dy(i-1,'3');
CON_tt(i)$(ord(i) gt 1 and ord(i) le nfe).. tt_0(i) =e= tt_0(i-1) + h(i-1)*sum(j, a(j,'3'));

* Acceleration Equations
Equations
eq00001_dynamics(i,j)
;

eq00001_dynamics(i,j)$(ord(i) le nfe).. m*ddy(i,j) =e= m*g + lambda1(i);
*eq00001_dynamics(i)$(ord(i) le nfe).. m*ddy(i) =e= -m*g + lambda1(i);



* Time Equations
Equations
hEqnLo(i)
hEqnUp(i)
hTotal
;

hEqnLo(i)..      h(i) =g= 1e-5*time/card(i);
hEqnUp(i)..      h(i) =l= 1.2*time/card(i);
hTotal..         time =e= sum(i,h(i));

* CONTACT COMPLEMENTARITY
Equations
eq00001_friction(i,j)
eq00002_friction(i,j)
*eq00003_friction(i)

contact1_Adummy(i)
contact1_Bdummy(i)

compEqn1(i)
penTotEqn
;

eq00001_friction(i,j)$(ord(i) ge 1 ).. a_contact_1(i,j) - y(i,j) =e= 0;
eq00002_friction(i,j)$(ord(i) ge 1 ).. b_contact_1(i,j) - lambda1(i,j) =e= 0;
*eq00003_friction(i)$(ord(i) ge 1 ).. a_contact_1_0(i) - y_0(i) =e= 0;

contact1_Adummy(i)..      a_contact_1_pr(i) =e= sum(j,a_contact_1(i,j))+a_contact_1_0(i);
contact1_Bdummy(i)..      b_contact_1_pr(i) =e= sum(j,b_contact_1(i,j));


compEqn1(i)$(ord(i) lt card(i))..  a_contact_1_pr(i+1)*b_contact_1_pr(i) =l= eps1;

penTotEqn..              pen =e=  sum(i$(ord(i) lt card(i)), a_contact_1_pr(i+1)*lambda1(i) );


model pointMass_noSum /all /;





* INITIAL CONDITIONS
tt_0.fx(i)$(ord(i)=1) = 0;
y_0.fx(i)$(ord(i)=1) = 1;
*th1_0.fx(i)$(ord(i)=1) = 0;
*lambda1.fx(i,j)$(ord(i)=1) = 0;
*dy_0.fx(i)$(ord(i)=1) = 0.5;
dy_0.fx(i)$(ord(i)=1) =-4.5;
*h.l(i) = time/card(i);
*h.lo(i) = 1e-3*time/card(i);
*h.up(i) = 2*time/card(i);
option nlp = conopt;

*pointMass_noSum.iterlim = 0;
* RELAXATION
Solve pointMass_noSum minimizing phi1 using nlp;
*h.lo(i) = 1e-4*time/card(i);
*Solve pointMass_noSum minimizing phi1 using nlp;
$ontext
eps1 = eps1/10;
*pointMass_noSum.optfile = 3;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
$offtext

*Solve pointMass_noSum minimizing phi1 using nlp;
*$ontext
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
eps1 = eps1/10;
Solve pointMass_noSum minimizing phi1 using nlp;
*eps1 = eps1/10;
*Solve pointMass_noSum minimizing phi1 using nlp;
*$offtext
         put_utilities 'ren' / 'C:\Users\Amir\Documents\MATLAB\wheeledLeg\newWheeledLeg\collocationTest\pointMass_3ptFull.csv';
         put "t;y;dy;ddy;lambda1"/;
         put tt_0.l('1')";"y_0.l('1')";"dy_0.l('1')";"ddy.l('1','1')";"lambda1.l('1')";"/
         loop(i,
                 loop(j,
                         put tt.l(i,j)";"y.l(i,j)";"dy.l(i,j)";"ddy.l(i,j)";"lambda1.l(i)";"/);

                 );
*         );
         putclose;
*$offtext
