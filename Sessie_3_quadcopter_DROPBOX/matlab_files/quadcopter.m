%***************************************
% Session 3: Quadcopter
%***************************************
clear all
close all
clc
load('references_14.mat')

%******************
% Parameters
%******************

m = 0.5;
L = 0.25;
k = 3*10^-6;
b = 10^-7;
g = 9.81;
k_d = 0.25;
I_xx = 5*10^-3;
I_yy = 5*10^-3;
I_zz = 10^-2;
c_m = 10^4;

%%
%******************
% 4.1 Linearization
%******************

w = sqrt(m*g/(4*k));
v = w/sqrt(c_m);
vsq = v^2;

% Linearization/ equilibrium point
u_e = [v^2 v^2 v^2 v^2].';
x_e = [0 0 0 0 0 0 0 0 0 0 0 0 ].';

syms x y z v_x v_y v_z phi theta psi w_x w_y w_z v_1 v_2 v_3 v_4
eq1=v_x;
eq2=v_y;
eq3=v_z;
eq4=(-k_d/m)*v_x + (k*c_m/m)*(sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta))*(v_1+v_2+v_3+v_4);
eq5=(-k_d/m)*v_y + (k*c_m/m)*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi))*(v_1+v_2+v_3+v_4);
eq6=(-k_d/m)*v_z -g + (k*c_m/m)*(cos(theta)*cos(phi))*(v_1+v_2+v_3+v_4);
eq7= w_x +w_y*(sin(phi)*tan(theta))+w_z*(cos(phi)*tan(theta));
eq8= w_y*cos(phi)-w_z*sin(phi);
eq9= w_y*sin(phi)/cos(theta) + w_z*cos(phi)/cos(theta);
eq10= (L*k*c_m/I_xx)*(v_1 - v_3) - ((I_yy-I_zz)/I_xx)*w_y*w_z;
eq11= (L*k*c_m/I_yy)*(v_2 - v_4) - ((I_zz-I_xx)/I_yy)*w_x*w_z;
eq12= (b*c_m/I_zz)*(v_1 - v_2 + v_3 - v_4) - ((I_xx-I_yy)/I_zz)*w_x*w_y;

f = [eq1; eq2; eq3;eq4;eq5;eq6;eq7;eq8;eq9;eq10;eq11;eq12];
j_x = [x, y, z, v_x, v_y, v_z, phi, theta, psi, w_x, w_y, w_z]; 
j_u = [v_1, v_2, v_3, v_4];
A_j = jacobian(f, j_x);
B_j = jacobian(f, j_u);

[x, y, z, v_x,v_y,v_z,phi,theta,psi,w_x,w_y,w_z] = deal(0);
[v_1,v_2,v_3,v_4] = deal(vsq);

% System matrices
A = eval(A_j);
B = eval(B_j);
C = zeros(6,12);
C(1:3,1:3) = eye(3);
C(4:6,7:9) = eye(3);
D = zeros(6,4);

%% 
%***************************************
% 4.2 Discretization: bilinear transf.
%***************************************

Ts = 0.05;
A_d = (eye(12)-A*Ts/2)^-1*(eye(12)+ A*Ts/2);
B_d = (eye(12)-A*Ts/2)^-1*B*Ts;
C_d = C*(eye(12)-A*Ts/2)^-1;
D_d = D + C*(eye(12)-A*Ts/2)^-1*B*Ts/2;

% Creation of a state-space object
sys = ss (A_d,B_d,C_d,D_d,Ts);

% Checking the stability: 
disp('Poles open loop:');
eig(A_d)

% Checking whether the system is controlable or not.
CO = ctrb(sys);
disp('Rank of the controllability matrix:');
rank(CO)

% Checking whether the system is observable or not.
Ob = obsv(sys);
disp('Rank of the observability matrix:');
rank(Ob)

% Transmission zeros 
disp('Transmission zeros:');
tzero(A_d,B_d,C_d,D_d) 

%% 
%*******************
% 4.3 a) LQR control
%*******************

[n_states,n_inputs] = size(B_d);
Q = eye(n_states);
Q(3,3)= 500;
R = 1*eye(n_inputs);
Q1 = eye(n_states);
Q1(3,3)= 700;
R1 = 0.01*eye(n_inputs);
% K: feedback gain for LQR without payload.
[K,S,e] = dlqr(A_d,B_d,Q,R);
% K1: feedback gain for LQR with payload.
[K1,S,e] = dlqr(A_d,B_d,Q1,R1);
% Analyse closed loop system (poles, step-response,...)
sys1 = ss(A_d-B_d*K,B_d,C_d-D_d*K,D_d,Ts);
sys2 = ss(A_d-B_d*K1,B_d,C_d-D_d*K1,D_d,Ts);
pzmap(sys1,sys2);
figure; step(sys1,sys2);

% Calculate Nx and Nu
M = [A_d-eye(12) B_d;
          C_d D_d];
V = [zeros(12,6);eye(6)];      
N = M\V;
Nx = N(1:12,:);
Nu = N(13:16,:);

%% 
%****************************************
% 4.3 b) LQR control with integral action
%****************************************

A_int = [eye(3) C_d(1:3,:)*Ts;
         zeros(12,3) A_d];
B_int = [D_d(1:3,:)*Ts;
         B_d];
C_int = [zeros(6,3) C_d];
D_int = D_d;
% Checking whether the augmented system is controllable or not.
CO_int = ctrb(A_int,B_int);
disp('Rank of A_int:');
rank(A_int)
disp('Rank of the controllability matrix of the augmented system:');
rank(CO_int)% = 15 = rank(A_int) -> all modes controllable -> system is also stabilizable

% LQR on the augmented system
[n_states,n_inputs] = size(B_int);
Q_int = eye(n_states);
Q_int(3,3)= 800;
Q_int(6,6)= 150;
R_int = 0.1*eye(n_inputs);
[K_int,S_int,e_int] = dlqr(A_int,B_int,Q_int,R_int);
sysint = ss(A_int-B_int*K_int,B_int,C_int-D_int*K_int,D_int,Ts);
disp('Poles closed loop integral  :');
eig(sysint)

%% 
%*****************
% 4.4 LQG control
%*****************

% Checking whether the system is observable or not
Ob = obsv(A_d,C_d);
disp('LQG: Rank of the observability matrix');
rank(Ob)

var_xyz = 2.5e-5;
var_angle = 7.57e-5;
% measurement noise covariance R:
R = diag([var_xyz, var_xyz, var_xyz, var_angle, var_angle, var_angle]);
% process noise covariance Q: 
Q = 0.1*eye(size(B_d,2));

[M,P,Z,E] = dlqe(A_d,B_d,C_d,Q,R);
L = A_d*M;

