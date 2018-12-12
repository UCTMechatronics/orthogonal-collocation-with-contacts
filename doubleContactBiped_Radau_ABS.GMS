$Title Optimal Control of a Double Contact Biped  (ABSOLUTE ANGLES)
file test / test.txt/;
test.pw=1000;
test.nd=6;
put test;

set univ the universe /1*1000/;
parameter col /1/;

sets j(univ)
    lastJ(univ);

j(univ)$(ord(univ) <= col) = yes;
lastJ(univ)$(ord(univ)=col)=yes;

Sets    i     number of finite elements             /1*300/
*        j     number of internal collocation points /1*1/;

Alias (j,k);

Scalar  tTime                                /3/
        ;

parameter a(univ,univ);
        a('1','1') = 1.0;

Parameter
        m0/31.72/
        m1/10.34/
        m2/3.16/
        m3/1/
        mu/0.8/
        c/0.4/
        g/9.81/
*        L0/0/
        L0/0.52/
        L1/0.42/
        L2/0.44/
        L3/0.25/
*         d1/0.07/
        d2/0.1655/
        d3/0.0922/
        alpha1/0.4366/
        alpha2/0.8622/
*        L3/0.25/
        I0/0.7148/
        I1/0.1520/
        I2/0.0510/
        I3/0.0184/
        T1Max/400/
        T2Max/250/
        T3Max/150/
        w1SS/12/
        w2SS/12/
        w3SS/10/
*    gear/1/
        solvercounter/1/
        row0/1E0/
        iter/0/
        converge/0/
        eps1/1E3/
        gearVal/1/
        armR/18E-3/
        kPhi/0.1/
        staticP/10/
        nfe
        ncp
        nBW/5/
        velScale/1/
        ;
*h.l(i) = 1/card(i) ;
nfe = card(i);
ncp = card(j);

Variables
*    JConst
    phi1
    phi2
*    P(i,univ)
*    energy_0(i)
*    energy(i,univ)
    x(i,univ)
    y(i,univ)
    th0(i,univ)
    th1(i,univ)
    th2(i,univ)
    th3(i,univ)
    th4(i,univ)
    th5(i,univ)
    th6(i,univ)
    dx(i,univ)
    dy(i,univ)
    dth0(i,univ)
    dth1(i,univ)
    dth2(i,univ)
    dth3(i,univ)
    dth4(i,univ)
    dth5(i,univ)
    dth6(i,univ)
    ddx(i,univ)
    ddy(i,univ)
    ddth0(i,univ)
    ddth1(i,univ)
    ddth2(i,univ)
    ddth3(i,univ)
    ddth4(i,univ)
    ddth5(i,univ)
    ddth6(i,univ)
    t1(i)
    t2(i)
    t3(i)
    t4(i)
    t5(i)
    t6(i)

    x_0(i)
    y_0(i)
    th0_0(i)
    th1_0(i)
    th2_0(i)
    th3_0(i)
    th4_0(i)
    th5_0(i)
    th6_0(i)
    dx_0(i)
    dy_0(i)
    dth0_0(i)
    dth1_0(i)
    dth2_0(i)
    dth3_0(i)
    dth4_0(i)
    dth5_0(i)
    dth6_0(i)

    dummyPit(i)
;

Positive Variables
         eCost
         taxPen
         h(i)
*         tCost
         tt_0(i)
         tt(i,univ)


         lambda_y1(i,univ)
         lambda_x1_p(i,univ)
         lambda_x1_n(i,univ)
         a_contact_1(i,univ)
         a_contact_1_pr(i)
         b_contact_1(i,univ)
         b_contact_1_pr(i)
         a_friction_1(i,univ)
         a_friction_1_pr(i)
         b_friction_1(i,univ)
         b_friction_1_pr(i)
         a_slip_p_1(i,univ)
         a_slip_p_1_pr(i)
         b_slip_p_1(i,univ)
         b_slip_p_1_pr(i)
         a_slip_n_1(i,univ)
         a_slip_n_1_pr(i)
         b_slip_n_1(i,univ)
         b_slip_n_1_pr(i)
         gamma_1(i,univ)

         lambda_y2(i,univ)
         lambda_x2_p(i,univ)
         lambda_x2_n(i,univ)
         a_contact_2(i,univ)
         a_contact_2_pr(i)
         b_contact_2(i,univ)
         b_contact_2_pr(i)
         a_friction_2(i,univ)
         a_friction_2_pr(i)
         b_friction_2(i,univ)
         b_friction_2_pr(i)
         a_slip_p_2(i,univ)
         a_slip_p_2_pr(i)
         b_slip_p_2(i,univ)
         b_slip_p_2_pr(i)
         a_slip_n_2(i,univ)
         a_slip_n_2_pr(i)
         b_slip_n_2(i,univ)
         b_slip_n_2_pr(i)
         gamma_2(i,univ)

         lambda_y3(i,univ)
         lambda_x3_p(i,univ)
         lambda_x3_n(i,univ)
         a_contact_3(i,univ)
         a_contact_3_pr(i)
         b_contact_3(i,univ)
         b_contact_3_pr(i)
         a_friction_3(i,univ)
         a_friction_3_pr(i)
         b_friction_3(i,univ)
         b_friction_3_pr(i)
         a_slip_p_3(i,univ)
         a_slip_p_3_pr(i)
         b_slip_p_3(i,univ)
         b_slip_p_3_pr(i)
         a_slip_n_3(i,univ)
         a_slip_n_3_pr(i)
         b_slip_n_3(i,univ)
         b_slip_n_3_pr(i)
         gamma_3(i,univ)

         lambda_y4(i,univ)
         lambda_x4_p(i,univ)
         lambda_x4_n(i,univ)
         a_contact_4(i,univ)
         a_contact_4_pr(i)
         b_contact_4(i,univ)
         b_contact_4_pr(i)
         a_friction_4(i,univ)
         a_friction_4_pr(i)
         b_friction_4(i,univ)
         b_friction_4_pr(i)
         a_slip_p_4(i,univ)
         a_slip_p_4_pr(i)
         b_slip_p_4(i,univ)
         b_slip_p_4_pr(i)
         a_slip_n_4(i,univ)
         a_slip_n_4_pr(i)
         b_slip_n_4(i,univ)
         b_slip_n_4_pr(i)
         gamma_4(i,univ)

;
*         gamma_2_0(i);


Equation fobj1
*         fobjConst
;
$ontext

$offtext
fobj1..  phi1 =e= 1;
*fobjConst..  JConst =e= 1;
*( sqr(t1(i)) + sqr(t2(i))+ sqr(t3(i)) + sqr(t4(i)) + sqr(t5(i)) + sqr(t6(i)))*h(i) )

* TOTAL TIME
Equation
*totalTime_eqn1
stepLim_eqn1(i)
stepLim_eqn2(i);
*totalTime_eqn2;

*totalTime_eqn1.. tTime =e= sum(i,tTime*h(i));
stepLim_eqn1(i)..  h(i) =g= 0.5*tTime/card(i);
stepLim_eqn2(i)..  h(i) =l= 1.5*tTime/card(i);

* JOINT LIMITS
Equation
th1UpEqn(i)
th1LoEqn(i)
th2UpEqn(i)
th2LoEqn(i)
th3UpEqn(i)
th3LoEqn(i)
th4UpEqn(i)
th4LoEqn(i)
th5UpEqn(i)
th5LoEqn(i)
th6UpEqn(i)
th6LoEqn(i)
;
th1UpEqn(i).. th1_0(i)-th0_0(i) =l= pi/2;
th1LoEqn(i).. th1_0(i)-th0_0(i) =g= -pi/2;
th2UpEqn(i).. th2_0(i)-th1_0(i) =l= pi-(15*pi/180);
th2LoEqn(i).. th2_0(i)-th1_0(i) =g= 0;
th3UpEqn(i).. th3_0(i)-th2_0(i) =l= 135*pi/180;
th3LoEqn(i).. th3_0(i)-th2_0(i) =g= 60*pi/180;

th4UpEqn(i).. th4_0(i)-th0_0(i) =l= pi/2;
th4LoEqn(i).. th4_0(i)-th0_0(i) =g= -pi/2;
th5UpEqn(i).. th5_0(i)-th4_0(i) =l= pi-(15*pi/180);
th5LoEqn(i).. th5_0(i)-th4_0(i) =g= 0;
th6UpEqn(i).. th6_0(i)-th5_0(i) =l= 135*pi/180;
th6LoEqn(i).. th6_0(i)-th5_0(i) =g= 60*pi/180;

*$ontext
* JOINT ANGULAR VELOCITY LIMITS
Equation
dth1UpEqn(i)
dth1LoEqn(i)
dth2UpEqn(i)
dth2LoEqn(i)
dth3UpEqn(i)
dth3LoEqn(i)
dth4UpEqn(i)
dth4LoEqn(i)
dth5UpEqn(i)
dth5LoEqn(i)
dth6UpEqn(i)
dth6LoEqn(i)
;
dth1UpEqn(i).. dth1_0(i)-dth0_0(i) =l= velScale*w1SS;
dth1LoEqn(i).. dth1_0(i)-dth0_0(i) =g= -velScale*w1SS;
dth2UpEqn(i).. dth2_0(i)-dth1_0(i) =l= velScale*w2SS;
dth2LoEqn(i).. dth2_0(i)-dth1_0(i) =g= -velScale*w2SS;
dth3UpEqn(i).. dth3_0(i)-dth2_0(i) =l= velScale*w3SS;
dth3LoEqn(i).. dth3_0(i)-dth2_0(i) =g= -velScale*w3SS;

dth4UpEqn(i).. dth4_0(i)-dth0_0(i) =l= velScale*w1SS;
dth4LoEqn(i).. dth4_0(i)-dth0_0(i) =g= -velScale*w1SS;
dth5UpEqn(i).. dth5_0(i)-dth4_0(i) =l= velScale*w2SS;
dth5LoEqn(i).. dth5_0(i)-dth4_0(i) =g= -velScale*w2SS;
dth6UpEqn(i).. dth6_0(i)-dth5_0(i) =l= velScale*w3SS;
dth6LoEqn(i).. dth6_0(i)-dth5_0(i) =g= -velScale*w3SS;
*$offtext


* TAX DAY CONSTRAINTS
Equation
taxPen_x_eqn(i)
*$ontext
taxPen_y_eqn(i)
taxPen_th0_eqn(i)
taxPen_th1_eqn(i)
taxPen_th2_eqn(i)
*taxPen_th3_eqn(i)
taxPen_th4_eqn(i)
taxPen_th5_eqn(i)
*taxPen_th6_eqn(i)
*taxPen_dth0_eqn(i)
taxPen_dx_eqn(i)
*taxPen_dy_eqn(i)
*$offtext
;

taxPen_x_eqn(i)$(ord(i) eq card(i)).. x_0(i) =e= 10;

*$ontext
taxPen_y_eqn(i)$(ord(i) eq card(i)).. y_0(i) =e= 1.12;

taxPen_th0_eqn(i)$(ord(i) eq card(i)).. th0_0(i) =e= 0;
taxPen_th1_eqn(i)$(ord(i) eq card(i)).. th1_0(i) =e= 0;
taxPen_th2_eqn(i)$(ord(i) eq card(i)).. th2_0(i) =e= 0;
*taxPen_th3_eqn(i)$(ord(i) eq card(i)).. th3_0(i) =e= pi/2;
taxPen_th4_eqn(i)$(ord(i) eq card(i)).. th4_0(i) =e= 0;
taxPen_th5_eqn(i)$(ord(i) eq card(i)).. th5_0(i) =e= 0;
*taxPen_th6_eqn(i)$(ord(i) eq card(i)).. th6_0(i) =e= pi/2;
*taxPen_dth0_eqn(i)$(ord(i) eq card(i)).. dth0_0(i) =e= 0;
taxPen_dx_eqn(i)$(ord(i) eq card(i)).. dx_0(i) =e= 0;
*taxPen_dy_eqn(i)$(ord(i) eq card(i)).. dy_0(i) =e= 0;
*$offtext
*$offtext

Equations
pitEqn(i)
foot1Pit(i)
foot2Pit(i)
foot3Pit(i)
foot4Pit(i)
;

pitEqn(i)..        dummyPit(i) =e= -0.4*x_0(i)*x_0(i)+ 4*x_0(i) -9.9;
*pitEqn(i)..        dummyPit(i) =e= -0.1778*x_0(i)*x_0(i)+ 3.5556*x_0(i) -17.6778;
foot1Pit(i)..      a_contact_1_pr(i) =g= dummyPit(i);
foot2Pit(i)..      a_contact_2_pr(i) =g= dummyPit(i);
foot3Pit(i)..      a_contact_3_pr(i) =g= dummyPit(i);
foot4Pit(i)..      a_contact_4_pr(i) =g= dummyPit(i);
*$offtext

$ontext
Equations
knee1Eqn(i)
knee2Eqn(i);

knee1Eqn(i).. y_0(i) - L1*cos(th1_0(i)) - (L0*cos(th0_0(i)))/2 =g= 0;
knee2Eqn(i).. y_0(i) - L1*cos(th4_0(i)) - (L0*cos(th0_0(i)))/2 =g= 0;
$offtext

Equations
FECOL_x(i,univ)
FECOL_y(i,univ)
FECOL_th0(i,univ)
FECOL_th1(i,univ)
FECOL_th2(i,univ)
FECOL_th3(i,univ)
FECOL_th4(i,univ)
FECOL_th5(i,univ)
FECOL_th6(i,univ)
FECOL_dx(i,univ)
FECOL_dy(i,univ)
FECOL_dth0(i,univ)
FECOL_dth1(i,univ)
FECOL_dth2(i,univ)
FECOL_dth3(i,univ)
FECOL_dth4(i,univ)
FECOL_dth5(i,univ)
FECOL_dth6(i,univ)


FECOL_tt(i,univ)

*$ontext
CON_x(i,univ)
CON_y(i,univ)
CON_th0(i,univ)
CON_th1(i,univ)
CON_th2(i,univ)
CON_th3(i,univ)
CON_th4(i,univ)
CON_th5(i,univ)
CON_th6(i,univ)
CON_dx(i,univ)
CON_dy(i,univ)
CON_dth0(i,univ)
CON_dth1(i,univ)
CON_dth2(i,univ)
CON_dth3(i,univ)
CON_dth4(i,univ)
CON_dth5(i,univ)
CON_dth6(i,univ)
CON_tt(i,univ)
*$offtext

$ontext
CON_x(i)
CON_y(i)
CON_th0(i)
CON_th1(i)
CON_th2(i)
CON_th3(i)
CON_th4(i)
CON_th5(i)
CON_th6(i)
CON_dx(i)
CON_dy(i)
CON_dth0(i)
CON_dth1(i)
CON_dth2(i)
CON_dth3(i)
CON_dth4(i)
CON_dth5(i)
CON_dth6(i)
$offtext

*CON_tt(i)


*term_x(i)
*term_y(i)
*term_th0(i)
*term_th1(i)
*term_th2(i)
*term_th4(i)
*term_th5(i)
*term_dx(i)
*term_dth0(i)
;

FECOL_x(i,j)$(ord(i) le nfe).. x(i,j)=e=x_0(i) + h(i)*sum(k,a(k,j)*dx(i,k));
FECOL_y(i,j)$(ord(i) le nfe).. y(i,j)=e=y_0(i) + h(i)*sum(k,a(k,j)*dy(i,k));
FECOL_th0(i,j)$(ord(i) le nfe).. th0(i,j)=e=th0_0(i) + h(i)*sum(k,a(k,j)*dth0(i,k));
FECOL_th1(i,j)$(ord(i) le nfe).. th1(i,j)=e=th1_0(i) + h(i)*sum(k,a(k,j)*dth1(i,k));
FECOL_th2(i,j)$(ord(i) le nfe).. th2(i,j)=e=th2_0(i) + h(i)*sum(k,a(k,j)*dth2(i,k));
FECOL_th3(i,j)$(ord(i) le nfe).. th3(i,j)=e=th3_0(i) + h(i)*sum(k,a(k,j)*dth3(i,k));
FECOL_th4(i,j)$(ord(i) le nfe).. th4(i,j)=e=th4_0(i) + h(i)*sum(k,a(k,j)*dth4(i,k));
FECOL_th5(i,j)$(ord(i) le nfe).. th5(i,j)=e=th5_0(i) + h(i)*sum(k,a(k,j)*dth5(i,k));
FECOL_th6(i,j)$(ord(i) le nfe).. th6(i,j)=e=th6_0(i) + h(i)*sum(k,a(k,j)*dth6(i,k));
FECOL_dx(i,j)$(ord(i) le nfe).. dx(i,j)=e=dx_0(i) + h(i)*sum(k,a(k,j)*ddx(i,k));
FECOL_dy(i,j)$(ord(i) le nfe).. dy(i,j)=e=dy_0(i) + h(i)*sum(k,a(k,j)*ddy(i,k));
FECOL_dth0(i,j)$(ord(i) le nfe).. dth0(i,j)=e=dth0_0(i) + h(i)*sum(k,a(k,j)*ddth0(i,k));
FECOL_dth1(i,j)$(ord(i) le nfe).. dth1(i,j)=e=dth1_0(i) + h(i)*sum(k,a(k,j)*ddth1(i,k));
FECOL_dth2(i,j)$(ord(i) le nfe).. dth2(i,j)=e=dth2_0(i) + h(i)*sum(k,a(k,j)*ddth2(i,k));
FECOL_dth3(i,j)$(ord(i) le nfe).. dth3(i,j)=e=dth3_0(i) + h(i)*sum(k,a(k,j)*ddth3(i,k));
FECOL_dth4(i,j)$(ord(i) le nfe).. dth4(i,j)=e=dth4_0(i) + h(i)*sum(k,a(k,j)*ddth4(i,k));
FECOL_dth5(i,j)$(ord(i) le nfe).. dth5(i,j)=e=dth5_0(i) + h(i)*sum(k,a(k,j)*ddth5(i,k));
FECOL_dth6(i,j)$(ord(i) le nfe).. dth6(i,j)=e=dth6_0(i) + h(i)*sum(k,a(k,j)*ddth6(i,k));

FECOL_tt(i,j)$(ord(i) le nfe).. tt(i,j) =e= tt_0(i) +h(i)*sum(k,a(k,j));
*FECOL_energy(i,j)$(ord(i) le nfe).. energy(i,j)=e=energy_0(i) + h(i)*sum(k,a(k,j)*P(i,k));


*CON_x(i,lastJ)$(ord(i) NE 1).. x_0(i) =e= x(i-1,lastJ);

$ontext
CON_x(i,lastJ)$(ord(i) NE 1).. x_0(i) =e= x_0(i-1) + h(i-1)*sum(j, dx(i-1,j)*a(j,lastJ));
CON_y(i,lastJ)$((ord(i) NE 1)).. y_0(i) =e= y_0(i-1) + h(i-1)*sum(j, dy(i-1,j)*a(j,lastJ));
CON_th0(i,lastJ)$(ord(i) NE 1).. th0_0(i) =e= th0_0(i-1) + h(i-1)*sum(j, dth0(i-1,j)*a(j,lastJ));
CON_th1(i,lastJ)$(ord(i) NE 1).. th1_0(i) =e= th1_0(i-1) + h(i-1)*sum(j, dth1(i-1,j)*a(j,lastJ));
CON_th2(i,lastJ)$(ord(i) NE 1).. th2_0(i) =e= th2_0(i-1) + h(i-1)*sum(j, dth2(i-1,j)*a(j,lastJ));
CON_th3(i,lastJ)$(ord(i) NE 1).. th3_0(i) =e= th3_0(i-1) + h(i-1)*sum(j, dth3(i-1,j)*a(j,lastJ));
CON_th4(i,lastJ)$(ord(i) NE 1).. th4_0(i) =e= th4_0(i-1) + h(i-1)*sum(j, dth4(i-1,j)*a(j,lastJ));
CON_th5(i,lastJ)$(ord(i) NE 1).. th5_0(i) =e= th5_0(i-1) + h(i-1)*sum(j, dth5(i-1,j)*a(j,lastJ));
CON_th6(i,lastJ)$(ord(i) NE 1).. th6_0(i) =e= th6_0(i-1) + h(i-1)*sum(j, dth6(i-1,j)*a(j,lastJ));
CON_dx(i,lastJ)$(ord(i) NE 1).. dx_0(i) =e= dx_0(i-1) + h(i-1)*sum(j, ddx(i-1,j)*a(j,lastJ));
CON_dy(i,lastJ)$(ord(i) NE 1).. dy_0(i) =e= dy_0(i-1) + h(i-1)*sum(j, ddy(i-1,j)*a(j,lastJ));
CON_dth0(i,lastJ)$(ord(i) NE 1).. dth0_0(i) =e= dth0_0(i-1) + h(i-1)*sum(j, ddth0(i-1,j)*a(j,lastJ));
CON_dth1(i,lastJ)$(ord(i) NE 1).. dth1_0(i) =e= dth1_0(i-1) + h(i-1)*sum(j, ddth1(i-1,j)*a(j,lastJ));
CON_dth2(i,lastJ)$(ord(i) NE 1).. dth2_0(i) =e= dth2_0(i-1) + h(i-1)*sum(j, ddth2(i-1,j)*a(j,lastJ));
CON_dth3(i,lastJ)$(ord(i) NE 1).. dth3_0(i) =e= dth3_0(i-1) + h(i-1)*sum(j, ddth3(i-1,j)*a(j,lastJ));
CON_dth4(i,lastJ)$(ord(i) NE 1).. dth4_0(i) =e= dth4_0(i-1) + h(i-1)*sum(j, ddth4(i-1,j)*a(j,lastJ));
CON_dth5(i,lastJ)$(ord(i) NE 1).. dth5_0(i) =e= dth5_0(i-1) + h(i-1)*sum(j, ddth5(i-1,j)*a(j,lastJ));
CON_dth6(i,lastJ)$(ord(i) NE 1).. dth6_0(i) =e= dth6_0(i-1) + h(i-1)*sum(j, ddth6(i-1,j)*a(j,lastJ));
$offtext

*$ontext
CON_x(i,lastJ)$(ord(i) NE 1).. x_0(i) =e= x(i-1,lastJ);
CON_y(i,lastJ)$(ord(i) NE 1).. y_0(i) =e= y(i-1,lastJ);
CON_th0(i,lastJ)$(ord(i) NE 1).. th0_0(i) =e= th0(i-1,lastJ);
CON_th1(i,lastJ)$(ord(i) NE 1).. th1_0(i) =e= th1(i-1,lastJ);
CON_th2(i,lastJ)$(ord(i) NE 1).. th2_0(i) =e= th2(i-1,lastJ);
CON_th3(i,lastJ)$(ord(i) NE 1).. th3_0(i) =e= th3(i-1,lastJ);
CON_th4(i,lastJ)$(ord(i) NE 1).. th4_0(i) =e= th4(i-1,lastJ);
CON_th5(i,lastJ)$(ord(i) NE 1).. th5_0(i) =e= th5(i-1,lastJ);
CON_th6(i,lastJ)$(ord(i) NE 1).. th6_0(i) =e= th6(i-1,lastJ);

CON_dx(i,lastJ)$(ord(i) NE 1).. dx_0(i) =e= dx(i-1,lastJ);
CON_dy(i,lastJ)$(ord(i) NE 1).. dy_0(i) =e= dy(i-1,lastJ);
CON_dth0(i,lastJ)$(ord(i) NE 1).. dth0_0(i) =e= dth0(i-1,lastJ);
CON_dth1(i,lastJ)$(ord(i) NE 1).. dth1_0(i) =e= dth1(i-1,lastJ);
CON_dth2(i,lastJ)$(ord(i) NE 1).. dth2_0(i) =e= dth2(i-1,lastJ);
CON_dth3(i,lastJ)$(ord(i) NE 1).. dth3_0(i) =e= dth3(i-1,lastJ);
CON_dth4(i,lastJ)$(ord(i) NE 1).. dth4_0(i) =e= dth4(i-1,lastJ);
CON_dth5(i,lastJ)$(ord(i) NE 1).. dth5_0(i) =e= dth5(i-1,lastJ);
CON_dth6(i,lastJ)$(ord(i) NE 1).. dth6_0(i) =e= dth6(i-1,lastJ);
CON_tt(i,lastJ)$(ord(i) NE 1)..   tt_0(i) =e= tt(i-1,lastJ);
*$offtext

*CON_energy(i,lastJ)$(ord(i) NE 1).. energy_0(i) =e= energy_0(i-1) + h(i-1)*sum(j, P(i-1,j)*a(j,lastJ));


*CON_tt(i,lastJ)$(ord(i) NE 1).. tt_0(i) =e= tt_0(i-1) + h(i-1)*sum(j,a(j,lastJ));


*powerEqn(i,j)$(ord(i) le nfe).. P(i,j) =e= (t1(i)*t1(i) + t2(i)*t2(i) + t3(i)*t3(i) + t4(i)*t4(i) + t5(i)*t5(i) + t6(i)*t6(i))*1.8;

* Acceleration Equations
*$include "C:\Users\User\Documents\doubleContactAbs\_dynamics_OLD.txt"                                                                                                                                                                                                                                                                            LD.txt"

Equations
eq00001_dynamics(i,univ)
eq00002_dynamics(i,univ)
eq00003_dynamics(i,univ)
eq00004_dynamics(i,univ)
eq00005_dynamics(i,univ)
eq00006_dynamics(i,univ)
eq00007_dynamics(i,univ)
eq00008_dynamics(i,univ)
eq00009_dynamics(i,univ)
;

eq00001_dynamics(i,j).. lambda_x1_n(i,j) + lambda_x2_n(i,j) - lambda_x1_p(i,j) + lambda_x3_n(i,j) - lambda_x2_p(i,j) + lambda_x4_n(i,j) - lambda_x3_p(i,j) - lambda_x4_p(i,j) + m0*ddx(i,j) + 2*m1*ddx(i,j) + 2*m2*ddx(i,j) + 2*m3*ddx(i,j) + c*power(dx(i,j),2) - L0*m1*cos(th0(i,j))*ddth0(i,j) - L0*m2*cos(th0(i,j))*ddth0(i,j) - L0*m3*cos(th0(i,j))*ddth0(i,j) - (L1*m1*cos(th1(i,j))*ddth1(i,j))/2 - L1*m2*cos(th1(i,j))*ddth1(i,j) - L1*m3*cos(th1(i,j))*ddth1(i,j) - (L2*m2*cos(th2(i,j))*ddth2(i,j))/2 - L2*m3*cos(th2(i,j))*ddth2(i,j) - (L1*m1*cos(th4(i,j))*ddth4(i,j))/2 - L1*m2*cos(th4(i,j))*ddth4(i,j) - L1*m3*cos(th4(i,j))*ddth4(i,j) - (L2*m2*cos(th5(i,j))*ddth5(i,j))/2 - L2*m3*cos(th5(i,j))*ddth5(i,j) + L0*m1*sin(th0(i,j))*power(dth0(i,j),2) + L0*m2*sin(th0(i,j))*power(dth0(i,j),2) + L0*m3*sin(th0(i,j))*power(dth0(i,j),2) + (L1*m1*sin(th1(i,j))*power(dth1(i,j),2))/2 + L1*m2*sin(th1(i,j))*power(dth1(i,j),2) + L1*m3*sin(th1(i,j))*power(dth1(i,j),2) + (L2*m2*sin(th2(i,j))*power(dth2(i,j),2))/2 + L2*m3*sin(th2(i,j))*power(dth2(i,j),2) + (L1*m1*sin(th4(i,j))*power(dth4(i,j),2))/2 + L1*m2*sin(th4(i,j))*power(dth4(i,j),2) + L1*m3*sin(th4(i,j))*power(dth4(i,j),2) + (L2*m2*sin(th5(i,j))*power(dth5(i,j),2))/2 + L2*m3*sin(th5(i,j))*power(dth5(i,j),2) =e= 0;
eq00002_dynamics(i,j).. g*m0 - lambda_y2(i,j) - lambda_y3(i,j) - lambda_y4(i,j) - lambda_y1(i,j) + 2*g*m1 + 2*g*m2 + 2*g*m3 + m0*ddy(i,j) + 2*m1*ddy(i,j) + 2*m2*ddy(i,j) + 2*m3*ddy(i,j) + L0*m1*sin(th0(i,j))*ddth0(i,j) + L0*m2*sin(th0(i,j))*ddth0(i,j) + L0*m3*sin(th0(i,j))*ddth0(i,j) + (L1*m1*sin(th1(i,j))*ddth1(i,j))/2 + L1*m2*sin(th1(i,j))*ddth1(i,j) + L1*m3*sin(th1(i,j))*ddth1(i,j) + (L2*m2*sin(th2(i,j))*ddth2(i,j))/2 + L2*m3*sin(th2(i,j))*ddth2(i,j) + (L1*m1*sin(th4(i,j))*ddth4(i,j))/2 + L1*m2*sin(th4(i,j))*ddth4(i,j) + L1*m3*sin(th4(i,j))*ddth4(i,j) + (L2*m2*sin(th5(i,j))*ddth5(i,j))/2 + L2*m3*sin(th5(i,j))*ddth5(i,j) + L0*m1*cos(th0(i,j))*power(dth0(i,j),2) + L0*m2*cos(th0(i,j))*power(dth0(i,j),2) + L0*m3*cos(th0(i,j))*power(dth0(i,j),2) + (L1*m1*cos(th1(i,j))*power(dth1(i,j),2))/2 + L1*m2*cos(th1(i,j))*power(dth1(i,j),2) + L1*m3*cos(th1(i,j))*power(dth1(i,j),2) + (L2*m2*cos(th2(i,j))*power(dth2(i,j),2))/2 + L2*m3*cos(th2(i,j))*power(dth2(i,j),2) + (L1*m1*cos(th4(i,j))*power(dth4(i,j),2))/2 + L1*m2*cos(th4(i,j))*power(dth4(i,j),2) + L1*m3*cos(th4(i,j))*power(dth4(i,j),2) + (L2*m2*cos(th5(i,j))*power(dth5(i,j),2))/2 + L2*m3*cos(th5(i,j))*power(dth5(i,j),2) =e= 0;
eq00003_dynamics(i,j).. I0*ddth0(i,j) + (power(L0,2)*m1*ddth0(i,j))/2 + (power(L0,2)*m2*ddth0(i,j))/2 + (power(L0,2)*m3*ddth0(i,j))/2 - (L0*cos(th0(i,j))*lambda_x1_n(i,j))/2 - (L0*cos(th0(i,j))*lambda_x2_n(i,j))/2 + (L0*cos(th0(i,j))*lambda_x1_p(i,j))/2 - (L0*cos(th0(i,j))*lambda_x3_n(i,j))/2 + (L0*cos(th0(i,j))*lambda_x2_p(i,j))/2 - (L0*cos(th0(i,j))*lambda_x4_n(i,j))/2 + (L0*cos(th0(i,j))*lambda_x3_p(i,j))/2 + (L0*cos(th0(i,j))*lambda_x4_p(i,j))/2 - (L0*sin(th0(i,j))*lambda_y1(i,j))/2 - (L0*sin(th0(i,j))*lambda_y2(i,j))/2 - (L0*sin(th0(i,j))*lambda_y3(i,j))/2 - (L0*sin(th0(i,j))*lambda_y4(i,j))/2 - L0*m1*cos(th0(i,j))*ddx(i,j) - L0*m2*cos(th0(i,j))*ddx(i,j) - L0*m3*cos(th0(i,j))*ddx(i,j) + L0*m1*sin(th0(i,j))*ddy(i,j) + L0*m2*sin(th0(i,j))*ddy(i,j) + L0*m3*sin(th0(i,j))*ddy(i,j) + L0*g*m1*sin(th0(i,j)) + L0*g*m2*sin(th0(i,j)) + L0*g*m3*sin(th0(i,j)) - (L0*L1*m1*cos(th0(i,j))*sin(th1(i,j))*power(dth1(i,j),2))/4 + (L0*L1*m1*cos(th1(i,j))*sin(th0(i,j))*power(dth1(i,j),2))/4 - (L0*L1*m2*cos(th0(i,j))*sin(th1(i,j))*power(dth1(i,j),2))/2 + (L0*L1*m2*cos(th1(i,j))*sin(th0(i,j))*power(dth1(i,j),2))/2 - (L0*L1*m3*cos(th0(i,j))*sin(th1(i,j))*power(dth1(i,j),2))/2 + (L0*L1*m3*cos(th1(i,j))*sin(th0(i,j))*power(dth1(i,j),2))/2 - (L0*L2*m2*cos(th0(i,j))*sin(th2(i,j))*power(dth2(i,j),2))/4 + (L0*L2*m2*cos(th2(i,j))*sin(th0(i,j))*power(dth2(i,j),2))/4 - (L0*L2*m3*cos(th0(i,j))*sin(th2(i,j))*power(dth2(i,j),2))/2 + (L0*L2*m3*cos(th2(i,j))*sin(th0(i,j))*power(dth2(i,j),2))/2 - (L0*L1*m1*cos(th0(i,j))*sin(th4(i,j))*power(dth4(i,j),2))/4 + (L0*L1*m1*cos(th4(i,j))*sin(th0(i,j))*power(dth4(i,j),2))/4 - (L0*L1*m2*cos(th0(i,j))*sin(th4(i,j))*power(dth4(i,j),2))/2 + (L0*L1*m2*cos(th4(i,j))*sin(th0(i,j))*power(dth4(i,j),2))/2 - (L0*L1*m3*cos(th0(i,j))*sin(th4(i,j))*power(dth4(i,j),2))/2 + (L0*L1*m3*cos(th4(i,j))*sin(th0(i,j))*power(dth4(i,j),2))/2 - (L0*L2*m2*cos(th0(i,j))*sin(th5(i,j))*power(dth5(i,j),2))/4 + (L0*L2*m2*cos(th5(i,j))*sin(th0(i,j))*power(dth5(i,j),2))/4 - (L0*L2*m3*cos(th0(i,j))*sin(th5(i,j))*power(dth5(i,j),2))/2 + (L0*L2*m3*cos(th5(i,j))*sin(th0(i,j))*power(dth5(i,j),2))/2 + (L0*L1*m1*cos(th0(i,j))*cos(th1(i,j))*ddth1(i,j))/4 + (L0*L1*m2*cos(th0(i,j))*cos(th1(i,j))*ddth1(i,j))/2 + (L0*L1*m3*cos(th0(i,j))*cos(th1(i,j))*ddth1(i,j))/2 + (L0*L2*m2*cos(th0(i,j))*cos(th2(i,j))*ddth2(i,j))/4 + (L0*L2*m3*cos(th0(i,j))*cos(th2(i,j))*ddth2(i,j))/2 + (L0*L1*m1*cos(th0(i,j))*cos(th4(i,j))*ddth4(i,j))/4 + (L0*L1*m2*cos(th0(i,j))*cos(th4(i,j))*ddth4(i,j))/2 + (L0*L1*m3*cos(th0(i,j))*cos(th4(i,j))*ddth4(i,j))/2 + (L0*L2*m2*cos(th0(i,j))*cos(th5(i,j))*ddth5(i,j))/4 + (L0*L2*m3*cos(th0(i,j))*cos(th5(i,j))*ddth5(i,j))/2 + (L0*L1*m1*sin(th0(i,j))*sin(th1(i,j))*ddth1(i,j))/4 + (L0*L1*m2*sin(th0(i,j))*sin(th1(i,j))*ddth1(i,j))/2 + (L0*L1*m3*sin(th0(i,j))*sin(th1(i,j))*ddth1(i,j))/2 + (L0*L2*m2*sin(th0(i,j))*sin(th2(i,j))*ddth2(i,j))/4 + (L0*L2*m3*sin(th0(i,j))*sin(th2(i,j))*ddth2(i,j))/2 + (L0*L1*m1*sin(th0(i,j))*sin(th4(i,j))*ddth4(i,j))/4 + (L0*L1*m2*sin(th0(i,j))*sin(th4(i,j))*ddth4(i,j))/2 + (L0*L1*m3*sin(th0(i,j))*sin(th4(i,j))*ddth4(i,j))/2 + (L0*L2*m2*sin(th0(i,j))*sin(th5(i,j))*ddth5(i,j))/4 + (L0*L2*m3*sin(th0(i,j))*sin(th5(i,j))*ddth5(i,j))/2 =e= 0;
eq00004_dynamics(i,j).. I1*ddth1(i,j) - t1(i) + (power(L1,2)*m1*ddth1(i,j))/4 + power(L1,2)*m2*ddth1(i,j) + power(L1,2)*m3*ddth1(i,j) - L1*cos(th1(i,j))*lambda_x1_n(i,j) - L1*cos(th1(i,j))*lambda_x2_n(i,j) + L1*cos(th1(i,j))*lambda_x1_p(i,j) + L1*cos(th1(i,j))*lambda_x2_p(i,j) - L1*sin(th1(i,j))*lambda_y1(i,j) - L1*sin(th1(i,j))*lambda_y2(i,j) - (L1*m1*cos(th1(i,j))*ddx(i,j))/2 - L1*m2*cos(th1(i,j))*ddx(i,j) - L1*m3*cos(th1(i,j))*ddx(i,j) + (L1*m1*sin(th1(i,j))*ddy(i,j))/2 + L1*m2*sin(th1(i,j))*ddy(i,j) + L1*m3*sin(th1(i,j))*ddy(i,j) + (L1*g*m1*sin(th1(i,j)))/2 + L1*g*m2*sin(th1(i,j)) + L1*g*m3*sin(th1(i,j)) + (L0*L1*m1*cos(th0(i,j))*sin(th1(i,j))*power(dth0(i,j),2))/4 - (L0*L1*m1*cos(th1(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/4 + (L0*L1*m2*cos(th0(i,j))*sin(th1(i,j))*power(dth0(i,j),2))/2 - (L0*L1*m2*cos(th1(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 + (L0*L1*m3*cos(th0(i,j))*sin(th1(i,j))*power(dth0(i,j),2))/2 - (L0*L1*m3*cos(th1(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 - (L1*L2*m2*cos(th1(i,j))*sin(th2(i,j))*power(dth2(i,j),2))/2 + (L1*L2*m2*cos(th2(i,j))*sin(th1(i,j))*power(dth2(i,j),2))/2 - L1*L2*m3*cos(th1(i,j))*sin(th2(i,j))*power(dth2(i,j),2) + L1*L2*m3*cos(th2(i,j))*sin(th1(i,j))*power(dth2(i,j),2) + (L0*L1*m1*cos(th0(i,j))*cos(th1(i,j))*ddth0(i,j))/4 + (L0*L1*m2*cos(th0(i,j))*cos(th1(i,j))*ddth0(i,j))/2 + (L0*L1*m3*cos(th0(i,j))*cos(th1(i,j))*ddth0(i,j))/2 + (L1*L2*m2*cos(th1(i,j))*cos(th2(i,j))*ddth2(i,j))/2 + L1*L2*m3*cos(th1(i,j))*cos(th2(i,j))*ddth2(i,j) + (L0*L1*m1*sin(th0(i,j))*sin(th1(i,j))*ddth0(i,j))/4 + (L0*L1*m2*sin(th0(i,j))*sin(th1(i,j))*ddth0(i,j))/2 + (L0*L1*m3*sin(th0(i,j))*sin(th1(i,j))*ddth0(i,j))/2 + (L1*L2*m2*sin(th1(i,j))*sin(th2(i,j))*ddth2(i,j))/2 + L1*L2*m3*sin(th1(i,j))*sin(th2(i,j))*ddth2(i,j) =e= 0;
eq00005_dynamics(i,j).. I2*ddth2(i,j) - t2(i) + (power(L2,2)*m2*ddth2(i,j))/4 + power(L2,2)*m3*ddth2(i,j) - L2*cos(th2(i,j))*lambda_x1_n(i,j) - L2*cos(th2(i,j))*lambda_x2_n(i,j) + L2*cos(th2(i,j))*lambda_x1_p(i,j) + L2*cos(th2(i,j))*lambda_x2_p(i,j) - L2*sin(th2(i,j))*lambda_y1(i,j) - L2*sin(th2(i,j))*lambda_y2(i,j) - (L2*m2*cos(th2(i,j))*ddx(i,j))/2 - L2*m3*cos(th2(i,j))*ddx(i,j) + (L2*m2*sin(th2(i,j))*ddy(i,j))/2 + L2*m3*sin(th2(i,j))*ddy(i,j) + (L2*g*m2*sin(th2(i,j)))/2 + L2*g*m3*sin(th2(i,j)) + (L0*L2*m2*cos(th0(i,j))*sin(th2(i,j))*power(dth0(i,j),2))/4 - (L0*L2*m2*cos(th2(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/4 + (L0*L2*m3*cos(th0(i,j))*sin(th2(i,j))*power(dth0(i,j),2))/2 - (L0*L2*m3*cos(th2(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 + (L1*L2*m2*cos(th1(i,j))*sin(th2(i,j))*power(dth1(i,j),2))/2 - (L1*L2*m2*cos(th2(i,j))*sin(th1(i,j))*power(dth1(i,j),2))/2 + L1*L2*m3*cos(th1(i,j))*sin(th2(i,j))*power(dth1(i,j),2) - L1*L2*m3*cos(th2(i,j))*sin(th1(i,j))*power(dth1(i,j),2) + (L0*L2*m2*cos(th0(i,j))*cos(th2(i,j))*ddth0(i,j))/4 + (L0*L2*m3*cos(th0(i,j))*cos(th2(i,j))*ddth0(i,j))/2 + (L1*L2*m2*cos(th1(i,j))*cos(th2(i,j))*ddth1(i,j))/2 + L1*L2*m3*cos(th1(i,j))*cos(th2(i,j))*ddth1(i,j) + (L0*L2*m2*sin(th0(i,j))*sin(th2(i,j))*ddth0(i,j))/4 + (L0*L2*m3*sin(th0(i,j))*sin(th2(i,j))*ddth0(i,j))/2 + (L1*L2*m2*sin(th1(i,j))*sin(th2(i,j))*ddth1(i,j))/2 + L1*L2*m3*sin(th1(i,j))*sin(th2(i,j))*ddth1(i,j) =e= 0;
eq00006_dynamics(i,j).. I3*ddth3(i,j) - t3(i) + (2*L3*cos(th3(i,j))*lambda_x1_n(i,j))/3 - (L3*cos(th3(i,j))*lambda_x2_n(i,j))/3 - (2*L3*cos(th3(i,j))*lambda_x1_p(i,j))/3 + (L3*cos(th3(i,j))*lambda_x2_p(i,j))/3 + (2*L3*sin(th3(i,j))*lambda_y1(i,j))/3 - (L3*sin(th3(i,j))*lambda_y2(i,j))/3 =e= 0;
eq00007_dynamics(i,j).. I1*ddth4(i,j) - t4(i) + (power(L1,2)*m1*ddth4(i,j))/4 + power(L1,2)*m2*ddth4(i,j) + power(L1,2)*m3*ddth4(i,j) - L1*cos(th4(i,j))*lambda_x3_n(i,j) - L1*cos(th4(i,j))*lambda_x4_n(i,j) + L1*cos(th4(i,j))*lambda_x3_p(i,j) + L1*cos(th4(i,j))*lambda_x4_p(i,j) - L1*sin(th4(i,j))*lambda_y3(i,j) - L1*sin(th4(i,j))*lambda_y4(i,j) - (L1*m1*cos(th4(i,j))*ddx(i,j))/2 - L1*m2*cos(th4(i,j))*ddx(i,j) - L1*m3*cos(th4(i,j))*ddx(i,j) + (L1*m1*sin(th4(i,j))*ddy(i,j))/2 + L1*m2*sin(th4(i,j))*ddy(i,j) + L1*m3*sin(th4(i,j))*ddy(i,j) + (L1*g*m1*sin(th4(i,j)))/2 + L1*g*m2*sin(th4(i,j)) + L1*g*m3*sin(th4(i,j)) + (L0*L1*m1*cos(th0(i,j))*sin(th4(i,j))*power(dth0(i,j),2))/4 - (L0*L1*m1*cos(th4(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/4 + (L0*L1*m2*cos(th0(i,j))*sin(th4(i,j))*power(dth0(i,j),2))/2 - (L0*L1*m2*cos(th4(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 + (L0*L1*m3*cos(th0(i,j))*sin(th4(i,j))*power(dth0(i,j),2))/2 - (L0*L1*m3*cos(th4(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 - (L1*L2*m2*cos(th4(i,j))*sin(th5(i,j))*power(dth5(i,j),2))/2 + (L1*L2*m2*cos(th5(i,j))*sin(th4(i,j))*power(dth5(i,j),2))/2 - L1*L2*m3*cos(th4(i,j))*sin(th5(i,j))*power(dth5(i,j),2) + L1*L2*m3*cos(th5(i,j))*sin(th4(i,j))*power(dth5(i,j),2) + (L0*L1*m1*cos(th0(i,j))*cos(th4(i,j))*ddth0(i,j))/4 + (L0*L1*m2*cos(th0(i,j))*cos(th4(i,j))*ddth0(i,j))/2 + (L0*L1*m3*cos(th0(i,j))*cos(th4(i,j))*ddth0(i,j))/2 + (L1*L2*m2*cos(th4(i,j))*cos(th5(i,j))*ddth5(i,j))/2 + L1*L2*m3*cos(th4(i,j))*cos(th5(i,j))*ddth5(i,j) + (L0*L1*m1*sin(th0(i,j))*sin(th4(i,j))*ddth0(i,j))/4 + (L0*L1*m2*sin(th0(i,j))*sin(th4(i,j))*ddth0(i,j))/2 + (L0*L1*m3*sin(th0(i,j))*sin(th4(i,j))*ddth0(i,j))/2 + (L1*L2*m2*sin(th4(i,j))*sin(th5(i,j))*ddth5(i,j))/2 + L1*L2*m3*sin(th4(i,j))*sin(th5(i,j))*ddth5(i,j) =e= 0;
eq00008_dynamics(i,j).. I2*ddth5(i,j) - t5(i) + (power(L2,2)*m2*ddth5(i,j))/4 + power(L2,2)*                                                             m3*ddth5(i,j) - L2*cos(th5(i,j))*lambda_x3_n(i,j) - L2*cos(th5(i,j))*lambda_x4_n(i,j) + L2*cos(th5(i,j))*lambda_x3_p(i,j) + L2*cos(th5(i,j))*lambda_x4_p(i,j) - L2*sin(th5(i,j))*lambda_y3(i,j) - L2*sin(th5(i,j))*lambda_y4(i,j) - (L2*m2*cos(th5(i,j))*ddx(i,j))/2 - L2*m3*cos(th5(i,j))*ddx(i,j) + (L2*m2*sin(th5(i,j))*ddy(i,j))/2 + L2*m3*sin(th5(i,j))*ddy(i,j) + (L2*g*m2*sin(th5(i,j)))/2 + L2*g*m3*sin(th5(i,j)) + (L0*L2*m2*cos(th0(i,j))*sin(th5(i,j))*power(dth0(i,j),2))/4 - (L0*L2*m2*cos(th5(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/4 + (L0*L2*m3*cos(th0(i,j))*sin(th5(i,j))*power(dth0(i,j),2))/2 - (L0*L2*m3*cos(th5(i,j))*sin(th0(i,j))*power(dth0(i,j),2))/2 + (L1*L2*m2*cos(th4(i,j))*sin(th5(i,j))*power(dth4(i,j),2))/2 - (L1*L2*m2*cos(th5(i,j))*sin(th4(i,j))*power(dth4(i,j),2))/2 + L1*L2*m3*cos(th4(i,j))*sin(th5(i,j))*power(dth4(i,j),2) - L1*L2*m3*cos(th5(i,j))*sin(th4(i,j))*power(dth4(i,j),2) + (L0*L2*m2*cos(th0(i,j))*cos(th5(i,j))*ddth0(i,j))/4 + (L0*L2*m3*cos(th0(i,j))*cos(th5(i,j))*ddth0(i,j))/2 + (L1*L2*m2*cos(th4(i,j))*cos(th5(i,j))*ddth4(i,j))/2 + L1*L2*m3*cos(th4(i,j))*cos(th5(i,j))*ddth4(i,j) + (L0*L2*m2*sin(th0(i,j))*sin(th5(i,j))*ddth0(i,j))/4 + (L0*L2*m3*sin(th0(i,j))*sin(th5(i,j))*ddth0(i,j))/2 + (L1*L2*m2*sin(th4(i,j))*sin(th5(i,j))*ddth4(i,j))/2 + L1*L2*m3*sin(th4(i,j))*sin(th5(i,j))*ddth4(i,j) =e= 0;
eq00009_dynamics(i,j).. I3*ddth6(i,j) - t6(i) + (2*L3*cos(th6(i,j))*lambda_x3_n(i,j))/3 - (L3*cos(th6(i,j))*lambda_x4_n(i,j))/3 - (2*L3*cos(th6(i,j))*lambda_x3_p(i,j))/3 + (L3*cos(th6(i,j))*lambda_x4_p(i,j))/3 + (2*L3*sin(th6(i,j))*lambda_y3(i,j))/3 - (L3*sin(th6(i,j))*lambda_y4(i,j))/3 =e= 0;



model justDynamics /all/;

*$ontext
* MOTOR MODEL
Equations
eq00001_motor(i)
eq00002_motor(i)
eq00003_motor(i)
eq00004_motor(i)
eq00005_motor(i)
eq00006_motor(i)
eq00007_motor(i)
eq00008_motor(i)
eq00009_motor(i)
eq000010_motor(i)
eq000011_motor(i)
eq000012_motor(i);

eq00001_motor(i)$(ord(i) le nfe).. t1(i) =l= T1Max - T1Max/w1SS*( dth1_0(i)-dth0_0(i) );
eq00002_motor(i)$(ord(i) le nfe).. t1(i) =g= -T1Max - T1Max/w1SS*( dth1_0(i)-dth0_0(i) );
eq00003_motor(i)$(ord(i) le nfe).. t2(i) =l= T2Max - T2Max/w2SS*( dth2_0(i)-dth1_0(i) );
eq00004_motor(i)$(ord(i) le nfe).. t2(i) =g= -T2Max - T2Max/w2SS*( dth2_0(i)-dth1_0(i) );
eq00005_motor(i)$(ord(i) le nfe).. t3(i) =l= T3Max*gearVal - T3Max*gearVal*gearVal/w3SS*( dth3_0(i)-dth2_0(i) );
eq00006_motor(i)$(ord(i) le nfe).. t3(i) =g= -T3Max*gearVal - T3Max*gearVal*gearVal/w3SS*( dth3_0(i)-dth2_0(i) );


eq00007_motor(i)$(ord(i) le nfe).. t4(i) =l= T1Max - T1Max/w1SS*( dth4_0(i)-dth0_0(i) );
eq00008_motor(i)$(ord(i) le nfe).. t4(i) =g= -T1Max - T1Max/w1SS*( dth4_0(i)-dth0_0(i) );
eq00009_motor(i)$(ord(i) le nfe).. t5(i) =l= T2Max - T2Max/w2SS*(dth5_0(i)-dth4_0(i));
eq000010_motor(i)$(ord(i) le nfe).. t5(i) =g= -T2Max - T2Max/w2SS*(dth5_0(i)-dth4_0(i));
eq000011_motor(i)$(ord(i) le nfe).. t6(i) =l= T3Max*gearVal - T3Max*gearVal*gearVal/w3SS*( dth6_0(i)-dth5_0(i) );
eq000012_motor(i)$(ord(i) le nfe).. t6(i) =g= -T3Max*gearVal - T3Max*gearVal*gearVal/w3SS*( dth6_0(i)-dth5_0(i) );

Model justDynamics2 /all/;


* Slack variables
$include "C:\Users\User\Documents\doubleContactAbs\_friction_OLD.txt"
*$offtext


*$ontext
* Slack Equations
Equations
contact1_Adummy(i)
fric1_Adummy(i)
slipP1_Adummy(i)
slipN1_Adummy(i)
contact1_Bdummy(i)
fric1_Bdummy(i)
slipP1_Bdummy(i)
slipN1_Bdummy(i)

contact2_Adummy(i)
fric2_Adummy(i)
slipP2_Adummy(i)
slipN2_Adummy(i)
contact2_Bdummy(i)
fric2_Bdummy(i)
slipP2_Bdummy(i)
slipN2_Bdummy(i)

contact3_Adummy(i)
fric3_Adummy(i)
slipP3_Adummy(i)
slipN3_Adummy(i)
contact3_Bdummy(i)
fric3_Bdummy(i)
slipP3_Bdummy(i)
slipN3_Bdummy(i)

contact4_Adummy(i)
fric4_Adummy(i)
slipP4_Adummy(i)
slipN4_Adummy(i)
contact4_Bdummy(i)
fric4_Bdummy(i)
slipP4_Bdummy(i)
slipN4_Bdummy(i)

;

* LEG 1 - Toe
contact1_Adummy(i)..      a_contact_1_pr(i) =e= sum(j,a_contact_1(i,j));
fric1_Adummy(i)..         a_friction_1_pr(i) =e= sum(j, a_friction_1(i,j));
slipP1_Adummy(i)..        a_slip_p_1_pr(i) =e= sum(j, a_slip_p_1(i,j));
slipN1_Adummy(i)..        a_slip_n_1_pr(i) =e= sum(j, a_slip_n_1(i,j));

contact1_Bdummy(i)..      b_contact_1_pr(i) =e= sum(j,b_contact_1(i,j));
fric1_Bdummy(i)..         b_friction_1_pr(i) =e= sum(j, b_friction_1(i,j));
slipP1_Bdummy(i)..        b_slip_p_1_pr(i) =e= sum(j, b_slip_p_1(i,j));
slipN1_Bdummy(i)..        b_slip_n_1_pr(i) =e= sum(j, b_slip_n_1(i,j));

* LEG 1 - Heel
contact2_Adummy(i)..      a_contact_2_pr(i) =e= sum(j,a_contact_2(i,j));
fric2_Adummy(i)..         a_friction_2_pr(i) =e= sum(j, a_friction_2(i,j));
slipP2_Adummy(i)..        a_slip_p_2_pr(i) =e= sum(j, a_slip_p_2(i,j));
slipN2_Adummy(i)..        a_slip_n_2_pr(i) =e= sum(j, a_slip_n_2(i,j));

contact2_Bdummy(i)..      b_contact_2_pr(i) =e= sum(j,b_contact_2(i,j));
fric2_Bdummy(i)..         b_friction_2_pr(i) =e= sum(j, b_friction_2(i,j));
slipP2_Bdummy(i)..        b_slip_p_2_pr(i) =e= sum(j, b_slip_p_2(i,j));
slipN2_Bdummy(i)..        b_slip_n_2_pr(i) =e= sum(j, b_slip_n_2(i,j));


* LEG 2 - Toe
contact3_Adummy(i)..      a_contact_3_pr(i) =e= sum(j,a_contact_3(i,j));
fric3_Adummy(i)..         a_friction_3_pr(i) =e= sum(j, a_friction_3(i,j));
slipP3_Adummy(i)..        a_slip_p_3_pr(i) =e= sum(j, a_slip_p_3(i,j));
slipN3_Adummy(i)..        a_slip_n_3_pr(i) =e= sum(j, a_slip_n_3(i,j));

contact3_Bdummy(i)..      b_contact_3_pr(i) =e= sum(j,b_contact_3(i,j));
fric3_Bdummy(i)..         b_friction_3_pr(i) =e= sum(j, b_friction_3(i,j));
slipP3_Bdummy(i)..        b_slip_p_3_pr(i) =e= sum(j, b_slip_p_3(i,j));
slipN3_Bdummy(i)..        b_slip_n_3_pr(i) =e= sum(j, b_slip_n_3(i,j));



* LEG 2 - Heel
contact4_Adummy(i)..      a_contact_4_pr(i) =e= sum(j,a_contact_4(i,j));
fric4_Adummy(i)..         a_friction_4_pr(i) =e= sum(j, a_friction_4(i,j));
slipP4_Adummy(i)..        a_slip_p_4_pr(i) =e= sum(j, a_slip_p_4(i,j));
slipN4_Adummy(i)..        a_slip_n_4_pr(i) =e= sum(j, a_slip_n_4(i,j));

contact4_Bdummy(i)..      b_contact_4_pr(i) =e= sum(j,b_contact_4(i,j));
fric4_Bdummy(i)..         b_friction_4_pr(i) =e= sum(j, b_friction_4(i,j));
slipP4_Bdummy(i)..        b_slip_p_4_pr(i) =e= sum(j, b_slip_p_4(i,j));
slipN4_Bdummy(i)..        b_slip_n_4_pr(i) =e= sum(j, b_slip_n_4(i,j));



Model justDynamics3 /all/;

* COMPLEMENTARITY EQNS - Contacts
Equations
leg1_contact_eqn(i)
leg2_contact_eqn(i)
leg3_contact_eqn(i)
leg4_contact_eqn(i)

;

leg1_contact_eqn(i)$(ord(i) lt card(i))..    a_contact_1_pr(i+1)*b_contact_1_pr(i) =l= eps1;
leg2_contact_eqn(i)$(ord(i) lt card(i))..    a_contact_2_pr(i+1)*b_contact_2_pr(i) =l= eps1;
leg3_contact_eqn(i)$(ord(i) lt card(i))..    a_contact_3_pr(i+1)*b_contact_3_pr(i) =l= eps1;
leg4_contact_eqn(i)$(ord(i) lt card(i))..    a_contact_4_pr(i+1)*b_contact_4_pr(i) =l= eps1;

Model justDynamics4 /all/;

* COMPLEMENTARITY EQNS  - Friction Cone
Equations
leg1_friction_eqn(i)
leg2_friction_eqn(i)
leg3_friction_eqn(i)
leg4_friction_eqn(i)
;

leg1_friction_eqn(i)..   a_friction_1_pr(i)*b_friction_1_pr(i) =l= eps1;
leg2_friction_eqn(i)..   a_friction_2_pr(i)*b_friction_2_pr(i) =l= eps1;
leg3_friction_eqn(i)..   a_friction_3_pr(i)*b_friction_3_pr(i) =l= eps1;
leg4_friction_eqn(i)..   a_friction_4_pr(i)*b_friction_4_pr(i) =l= eps1;

Model justDynamics5 /all/;

* COMPLEMENTARITY EQNS - Slipping
Equations
leg1_slipP_eqn(i)
leg1_slipN_eqn(i)

leg2_slipP_eqn(i)
leg2_slipN_eqn(i)

leg3_slipP_eqn(i)
leg3_slipN_eqn(i)

leg4_slipP_eqn(i)
leg4_slipN_eqn(i)
;

leg1_slipP_eqn(i)..      a_slip_p_1_pr(i)*b_slip_p_1_pr(i) =l= eps1;
leg1_slipN_eqn(i)..      a_slip_n_1_pr(i)*b_slip_n_1_pr(i) =l= eps1;

leg2_slipP_eqn(i)..      a_slip_p_2_pr(i)*b_slip_p_2_pr(i) =l= eps1;
leg2_slipN_eqn(i)..      a_slip_n_2_pr(i)*b_slip_n_2_pr(i) =l= eps1;

leg3_slipP_eqn(i)..      a_slip_p_3_pr(i)*b_slip_p_3_pr(i) =l= eps1;
leg3_slipN_eqn(i)..      a_slip_n_3_pr(i)*b_slip_n_3_pr(i) =l= eps1;

leg4_slipP_eqn(i)..      a_slip_p_4_pr(i)*b_slip_p_4_pr(i) =l= eps1;
leg4_slipN_eqn(i)..      a_slip_n_4_pr(i)*b_slip_n_4_pr(i) =l= eps1;

Model justDynamics6 /all/;

Equation fobj2
*energyEqn
;
*energyEqn..      eCost =e=  sum(i, ( sqr(t1(i)) + sqr(t2(i))+ sqr(t3(i)) + sqr(t4(i)) + sqr(t5(i)) + sqr(t6(i)) )*1e-4 );
fobj2..          phi2 =e=  sum(i, ( sqr(t1(i)) + sqr(t2(i))+ sqr(t3(i)) + sqr(t4(i)) + sqr(t5(i)) + sqr(t6(i)) )*h(i) );

*sum(i, sqr( h(i)*card(i)/tTime ));

* + 1e3*sum(i, s1(i) + s2(i) +s3(i) +s4(i) +s5(i) +s6(i) +s7(i) +s8(i)+ s9(i) + s10(i) +s11(i) +s12(i) +s13(i) +s14(i) +s15(i) +s16(i));

*sum(i, sqr( h(i)*card(i)/tTime ))

*(sum(i,( t1(i)*t1(i) +t2(i)*t2(i) +t3(i)*t3(i) +t4(i)*t4(i) +t5(i)*t5(i) +t6(i)*t6(i) )*h(i)*armR/(kPhi*kPhi)))/((m0+2*m1+2*m2+2*m3)*(g*10)) + staticP*sum(i,h(i));
*sum(i, ( sqr(t1(i)) + sqr(t2(i))+ sqr(t3(i)) + sqr(t4(i)) + sqr(t5(i)) + sqr(t6(i)))*h(i) );



*$offtext
*TTT.. time =e= transtime ;
option reslim= 12E3;
Model doubleContactRadau /all/;
*****************************************************************************************************************
*$ontext
* TRY THIRD ORDER COLLOCATION
col=1;
ncp = card(j);
j(univ)$(ord(univ) <= col) = yes;
lastJ(univ)$(ord(univ)<>col)=no;
lastJ(univ)$(ord(univ)=col)=yes;


t1.lo(i) = -T1Max; t1.up(i) = T1Max;
 t2.lo(i) = -T2Max; t2.up(i) = T2Max;
 t3.lo(i) = -T3Max; t3.up(i) = T3Max;

 t4.lo(i) = -T1Max; t4.up(i) = T1Max;
 t5.lo(i) = -T2Max; t5.up(i) = T2Max;
 t6.lo(i) = -T3Max; t6.up(i) = T3Max;

*$ontext
 x_0.lo(i) = -1; x_0.up(i) = 20;
 y_0.lo(i) = 0.7; y_0.up(i) = 2;
 th0_0.lo(i) = -pi/4; th0_0.up(i) = pi/4;
*$offtext

*$ontext
lambda_y1.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x1_p.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x1_n.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;

lambda_y2.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x2_p.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x2_n.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;

lambda_y3.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x3_p.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x3_n.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;

lambda_y4.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x4_p.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
lambda_x4_n.up(i,j) = nBW*(m0+2*m1+2*m2+2*m3)*g;
*$offtext

tt_0.fx(i)$(ord(i)=1) = 0;
*energy_0.fx(i)$(ord(i) =1) = 0;
x_0.fx(i)$(ord(i)=1) = 0;
y_0.fx(i)$(ord(i)=1) = L0/2+L1*cos(0*pi/180)+L2*cos(0*pi/180)+0.01;
th0_0.fx(i)$(ord(i)=1) = 0;
th1_0.fx(i)$(ord(i)=1) = 0;
th2_0.fx(i)$(ord(i)=1) = 0;
th3_0.fx(i)$(ord(i)=1) = pi/2;
th4_0.fx(i)$(ord(i)=1) = 0;
th5_0.fx(i)$(ord(i)=1) = 0;
th6_0.fx(i)$(ord(i)=1) = pi/2;
dx_0.fx(i)$(ord(i)=1) = 0;
dy_0.fx(i)$(ord(i)=1) = 0;
dth0_0.fx(i)$(ord(i)=1) = 0;
dth1_0.fx(i)$(ord(i)=1) = 0;
dth2_0.fx(i)$(ord(i)=1) = 0;
dth3_0.fx(i)$(ord(i)=1) = 0;
dth4_0.fx(i)$(ord(i)=1) = 0;
dth5_0.fx(i)$(ord(i)=1) = 0;
dth6_0.fx(i)$(ord(i)=1) = 0;

ddy.fx(i,'1')$(ord(i)=1) = -g;


x.l(i,'1')=(ord(i)-1)*10/(card(i)-1);
*x.l(i,j)=(ord(i)-1)*10/card(i);

y.l(i,j)=uniform(0.8*1.12,1.2*1.12);
th0.l(i,j)=uniform(th0_0.lo(i),th0_0.up(i));
th1.l(i,j)=uniform(-0.1,0.1);
         th2.l(i,j)=uniform(-0.1,0.1);
          th3.l(i,j)=uniform(-0.1,0.1);
          th4.l(i,j)=uniform(-0.1,0.1);
          th5.l(i,j)=uniform(-0.1,0.1);
          th6.l(i,j)=uniform(-0.1,0.1);
*th0.l(i,j)=uniform(th0_0.lo(i),th0_0.up(i));
dx.l(i,j)=0.01;
          dy.l(i,j)=0.01;
          dth0.l(i,j)=0.01;
          dth1.l(i,j)=0.01;
          dth2.l(i,j)=0.01;
          dth3.l(i,j)=0.01;
          ddx.l(i,j)=0.01;
          ddy.l(i,j)=0.01;
          ddth0.l(i,j)=0.01;
          ddth1.l(i,j)=0.01;
          ddth2.l(i,j)=0.01;
          ddth3.l(i,j)=0.01;
          t1.l(i)=0.01;
          t2.l(i)=0.01;
          t3.l(i)=0.01;
          t4.l(i)=0.01;
          t5.l(i)=0.01;
          t6.l(i)=0.01;
lambda_y1.l(i,j)=0.01;
          lambda_x1_p.l(i,j)=0.01;
          lambda_x1_n.l(i,j)=0.01;
          lambda_y2.l(i,j)=0.01;
          lambda_x2_p.l(i,j)=0.01;
          lambda_x2_n.l(i,j)=0.01;
         lambda_y3.l(i,j)=0.01;
          lambda_x3_p.l(i,j)=0.01;
          lambda_x3_n.l(i,j)=0.01;
          lambda_y4.l(i,j)=0.01;
          lambda_x4_p.l(i,j)=0.01;
          lambda_x4_n.l(i,j)=0.01;

$ontext
s1.l(i) = 1e3;
s2.l(i) = 1e3;
s3.l(i) = 1e3;
s4.l(i) = 1e3;
s5.l(i) = 1e3;
s6.l(i) = 1e3;
s7.l(i) = 1e3;
s8.l(i) = 1e3;
s9.l(i) = 1e3;
s10.l(i) = 1e3;
s11.l(i) = 1e3;
s12.l(i) = 1e3;
s13.l(i) = 1e3;
s14.l(i) = 1e3;
s15.l(i) = 1e3;
s16.l(i) = 1e3;
$offtext

$ontext




th1.l(i,j)=uniform(th1_0.lo(i),th1_0.up(i));

         th2.l(i,j)=uniform(th2_0.lo(i),th2_0.up(i));
          th3.l(i,j)=uniform(th3_0.lo(i),th3_0.up(i));
          th4.l(i,j)=uniform(th1_0.lo(i),th1_0.up(i));
          th5.l(i,j)=uniform(th2_0.lo(i),th2_0.up(i));
          th6.l(i,j)=uniform(th3_0.lo(i),th3_0.up(i));


$offtext

*$ontext
* 5 POINT RADAU MATRIX

$ontext
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


$ontext
* 3 POINT RADAU MATRIX
a('1','1') =0.19681547722366;
a('1','2') =0.39442431473909;
a('1','3') =0.37640306270047;
a('2','1') = -0.06553542585020;
a('2','2') =0.29207341166523;
a('2','3') =0.51248582618842;
a('3','1') = 0.02377097434822;
a('3','2') =-0.04154875212600;
a('3','3') =0.11111111111111;
$offtext

$ontext
a('1','1') =0.416666125000187;
a('1','2') =0.749999625000187;
a('2','1') =-0.083333125000187;
a('2','2') =0.250000374999812;
$offtext

a('1','1') = 1.0;

$ontext
justDynamics.scaleOpt = 1;
justDynamics2.scaleOpt = 1;
justDynamics3.scaleOpt = 1;
doubleContactRadau.scaleOpt = 1;
*$include "C:\Users\User\Documents\doubleContactAbs\doubleContactScaling.txt"
$offtext

converge=0;
iter =0;
*dummyPit.l(i,j)=dummyPit.l(i,'1');

option nlp = ipopth;
*h.fx(i) = h.l(i);
eps1 = 1e3;
*$ontext
*h.fx(i) =  tTime/card(i);
justDynamics.optfile = 2;
justDynamics2.optfile = 2;
justDynamics3.optfile = 2;
justDynamics6.optfile = 2;

*$offtext
h.fx(i)= 1*tTime/card(i);
Solve justDynamics minimizing phi1 using nlp;
Solve justDynamics2 minimizing phi1 using nlp;
Solve justDynamics3 minimizing phi1 using nlp;
Solve justDynamics6 minimizing phi1 using nlp;

*monoRadau.optfile=4;
doubleContactRadau.optfile = 4;
doubleContactRadau.threads = 16;

h.up(i) =  inf;
h.lo(i) =  0;


*solve doubleContactRadau using nlp minimizing phi1;
*eps1 = eps1/10;


while (converge = 0 and iter lt 8,

*Solve Full Optimal Solution at higher tolerance, while decreasing eps

               solve doubleContactRadau using nlp minimizing phi2;
*               doubleContactRadau.optfile=3;

               if (doubleContactRadau.Modelstat lt 3,
                  eps1 = eps1/10;

               else
                  abort$(doubleContactRadau.solvestat <> 1) "numerical errors in recipe simulation";
*                  converge = 0;


               );

            iter = iter+1;

);

         put_utilities 'ren' / 'C:\Users\User\Documents\doubleContactAbs\doubleContact_10m_withGap_1pt_300N.csv';
         put "t;h; x;y;th0;th1;th2;th3;th4;th5;th6;t1;t2;t3;t4;t5;t6;";
         put "lambda_y1;lambda_x1_p;lambda_x1_n;lambda_y2;lambda_x2_p;lambda_x2_n; lambda_y3;lambda_x3_p;lambda_x3_n;lambda_y4;lambda_x4_p;lambda_x4_n;";
         put "dx; dy; dth0; dth1; dth2; dth3;dth4; dth5; dth6; ddx; ddy; ddth0; ddth1; ddth2; ddth3;ddth4;ddth5; ddth6; "/;

         loop(i,
*                 put tt.l(i,j)";"x.l(i,j)";"y.l(i,j)";"th0.l(i,j)";"th1.l(i,j)";"th2.l(i,j)";"t1.l(i,j)";"dx.l(i,j)";"dy.l(i,j)";"dth0.l(i,j)";"dth1.l(i,j)";"dth2.l(i,j)";"ddx.l(i,j)";"ddy.l(i,j)";"ddth0.l(i,j)";"ddth1.l(i,j)";"ddth2.l(i,j)";"/);               lambda_x_p
                 loop(j,
                         put tt.l(i,j)";"h.l(i)";"x.l(i,j)";"y.l(i,j)";"th0.l(i,j)";"th1.l(i,j)";"th2.l(i,j)";"th3.l(i,j)";"th4.l(i,j)";"th5.l(i,j)";"th6.l(i,j)";"t1.l(i)";"t2.l(i)";"t3.l(i)";"t4.l(i)";"t5.l(i)";"t6.l(i)";";
                         put lambda_y1.l(i,j)";"lambda_x1_p.l(i,j)"; "lambda_x1_n.l(i,j)";"lambda_y2.l(i,j)"; "lambda_x2_p.l(i,j)"; "lambda_x2_n.l(i,j)";"lambda_y3.l(i,j)";"lambda_x3_p.l(i,j)"; "lambda_x3_n.l(i,j)";"lambda_y4.l(i,j)"; "lambda_x4_p.l(i,j)"; "lambda_x4_n.l(i,j)";";
                         put dx.l(i,j)"; "dy.l(i,j)"; "dth0.l(i,j)"; "dth1.l(i,j)"; "dth2.l(i,j)";"dth3.l(i,j)";"dth4.l(i,j)";"dth5.l(i,j)"; "dth6.l(i,j)";"ddx.l(i,j)"; "ddy.l(i,j)"; "ddth0.l(i,j)"; "ddth1.l(i,j)"; "ddth2.l(i,j)";"ddth3.l(i,j)";"ddth4.l(i,j)";"ddth5.l(i,j)"; "ddth6.l(i,j)";"/);
                 );
*         );
         putclose;
*abort "Radau Fixed Done!";

*$offtext

*display tt0.l, c0.l, t0.l;
* option nlp = conopt;
* Solve hicks minimizing phi using nlp;

