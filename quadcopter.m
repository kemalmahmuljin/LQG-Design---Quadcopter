clear all;
%% Loading references
load references_02.mat
n_states = 12;
n_inputs = 4;
n_outputs = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Linearization
%       First run this before test_linear_model.slx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Parameters
m = 0.5;        % Mass of QC
L = 0.25;       % Radius of Qc
k = 3e-6;       % Propeller lift coefficient
b = 1e-7;       % Propeller drag coefficient
g = 9.81;       % Acceleration of gravity
k_d = 0.25;     % Air friction coefficient 
I_xx = 5e-3;    % QC inertia x-axes
I_yy = 5e-3;    % QC inertia y-axes
I_zz = 1e-2;    % QC inertia z-axes
c_m = 1e4;      % Motor constant


%% Operating point
x = 0;
y = 0;
z = 0;
phi = 0;
theta = 0;
psi = 0;
w_x = 0;                        % From eq (11)
w_y = 0;                        % From eq (12)
w_z = 0;                        % From eq (13)
v_ss = m * g / ( 4 * k * c_m);  % From set of eqs

%% State-Space matrices
% system matrix
A = zeros(n_states);

A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,4) = -k_d / m;
A(4,8) = k * c_m / m * 4 * v_ss;
A(5,5) = -k_d / m;
A(5,7) = k * c_m / m * 4 * v_ss;
A(6,6) = -k_d / m;
A(7,10) = 1;
A(8,11) = 1;
A(9,12) = 1;

% input matrix
B = zeros(n_states,n_inputs);

B(6,1) = k * c_m / m; 
B(6,2) = k * c_m / m; 
B(6,3) = k * c_m / m; 
B(6,4) = k * c_m / m; 
B(10,1) = L * k * c_m / I_xx;
B(10,3) = -L * k * c_m / I_xx;
B(11,2) = L * k * c_m / I_yy;
B(11,4) = -L * k * c_m / I_yy;
B(12,1) = b * c_m / I_zz;
B(12,2) = -b * c_m / I_zz;
B(12,3) = b * c_m / I_zz;
B(12,4) = -b * c_m / I_zz;

% output matrix
C = zeros(n_outputs,n_states);

C(1,1) = 1;
C(2,2) = 1;
C(3,3) = 1;
C(4,7) = 1;
C(5,8) = 1;
C(6,9) = 1;

% NO feedthrough matrix
D = zeros(6,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Discretization
%       Here we discretize using :
%           Trapezoidal rule
%           Zero Order Hold
%           Forward Euler rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating system
sys= ss(A,B,C,D);
G = tf(sys);
% Sampling time
Ts = 0.05;

%% Trapezoidal Rule : bilinear
G_DT_bilinear = c2d(G,Ts,'tustin');

%% Zero order Hold
G_DT_zoh = c2d(G,Ts,'zoh');

%% Forward Euler
A_FE = eye( size(A,1) ) + A * Ts;
B_FE = B * Ts;
C_FE = C;
D_FE = D;

G_DT_FE = tf( ss(A_FE,B_FE,C_FE,D_FE,Ts) );


%%  manual trapezoidal
% Deriving transformation by hand to keep dimensions smaller
Ad = (eye(n_states) - A*Ts/2)\(eye(n_states) + A*Ts/2);
Bd = (eye(n_states) - A*Ts/2)\B*Ts;
Cd = C/(eye(n_states) - A*Ts/2);
Dd = D + C/(eye(n_states) - A*Ts/2) * B*Ts/2;
sysD = ss(Ad, Bd, Cd, Dd, Ts);

%% step and frequency responses for the different discretizations
%{
%% Step response
T_sim = 5;
t_sim_steps = [0:Ts:T_sim];
step(G,G_DT_bilinear,G_DT_zoh,G_DT_FE,t_sim_steps)
legend show

%% Frequency response
figure;
H=bodeplot(G,G_DT_bilinear,G_DT_zoh, G_DT_FE);
legend show
grid

pp = getoptions(H);
pp.YLim{1}=[-25 25];
setoptions(H,pp)

%% Or only visualize one transfer function:
% For 4 input to 1st output

%% Step response
T_sim = 5;
t_sim_steps = [0:Ts:T_sim];
step(G(4,1),G_DT_bilinear(4,1),G_DT_zoh(4,1),G_DT_FE(4,1),t_sim_steps)
legend('G41', 'G_DT_bilinear41', 'G_DT_zoh41', 'G_DT_FE41' );

%% Frequency response
figure;
H=bodeplot(G(4,1),G_DT_bilinear(4,1),G_DT_zoh(4,1), G_DT_FE(4,1));
legend('G41', 'G_DT_bilinear41', 'G_DT_zoh41', 'G_DT_FE41' );
grid
%}

%% Performing checks on discretisation
fprintf("========================================================================");
fprintf(['\n\nThe original continiuous time system itself is unstable (pole = 0). \n', ...
     'The poles of the continiuous system are :\n']);
disp(eig(A))
fprintf('\n\nThe poles of the discretized system using any transformation are: \n');
p_D = eig(sysD.A)
% p_T = pole(G_DT_bilinear)
% p_ZOH = pole(G_DT_zoh)

fprintf('\nNotice the same number of unstable poles.\n');
fprintf('\nNext we check if the discretized system is controllable and observable.');
fprintf(['\n\nWe can see that the rank of the controllability matrix is the same', ...
    'as the order of the system: \n']);
fprintf('\nOrder of the system is : %i \n', n_states);
fprintf('Rank of controllability matrix is : %i \n', rank(ctrb(sysD.A, sysD.B)));
fprintf('Rank of observability matrix is : %i \n', rank(obsv(sysD.A, sysD.C)));
fprintf(['\nThe system is thus controllable and observable.']);
fprintf('\n(Ad, Bd) is stabilizable, because all unstable nodes are controllable.\n');
fprintf('Furthermore, (Ad, Cd) is detectable, because all unstable nodes are observable.\n');
% Confirming that all nodes are controllable and observable using PBH test
[V_FE,~] = eig(sysD.A);
result_C_FE = sysD.C*V_FE;
[V_FET,~] = eig(sysD.A');
result_B_FE = sysD.B'*V_FET;
% Checking transmission zero's of discretized model
% tzero calculates the solution for the matrix-vector product
% at the bottom of page 82 in the course. This does not guarantee a rank
% decrease of G(e) that's why we need to check it for each solution
fprintf('\n\nHere we show the solutions to tzero for the different systems\n');
Continuous_T_zero = tzero(sys)
D_T_zero = tzero(sysD)
% bilinear_T_zero = tzero(ss(G_DT_bilinear))
% zoh_T_zero = tzero(ss(G_DT_zoh))
fprintf(['\n\nTrapezoidal rule and seems to have solutions, however we need to check if the are actually', ...
    'transmission zeros: \n']);
M_D = [D_T_zero(1)*eye(length(sysD.A))-sysD.A -sysD.B; sysD.C sysD.D];
z_D = size(null(M_D))
rank([5*eye(length(sysD.A))-sysD.A -sysD.B; sysD.C sysD.D])
% M_T = [bilinear_T_zero(2)*eye(length(A_FE))-A_FE -B_FE; C_FE D_FE];
% z_T = size(null(M_T))
% M_ZOH1 = [zoh_T_zero(1)*eye(length(ss(G_DT_zoh).A))-ss(G_DT_zoh).A -ss(G_DT_zoh).B; ss(G_DT_zoh).C ss(G_DT_zoh).D];
% z_ZOH1 = size(null(M_ZOH1))
% M_ZOH2 = [zoh_T_zero(2)*eye(length(ss(G_DT_zoh).A))-ss(G_DT_zoh).A -ss(G_DT_zoh).B; ss(G_DT_zoh).C ss(G_DT_zoh).D];
% z_ZOH2 = size(null(M_ZOH2))
fprintf('\n\nWe can thus confirm that the Trapezoidal discretization has 4 transmission zeros.\n');

% Performing Kalman Decomposition
[ADc,BDc,CDc,~,kDc] = ctrbf(sysD.A,sysD.B,sysD.C); 
[ADo,BDo,CDo,~,kDo] = obsvf(sysD.A,sysD.B,sysD.C);

[Abarc,Bbarc,Cbarc,~,kc] = ctrbf(A,B,C); 
[Abaro,Bbaro,Cbaro,~,ko] = obsvf(A,B,C);
fprintf(['\n\nLastly as the system is controllable and observable we know that the system is minimal' , ...
    'This is also confirmed after performing the Kalman decomposition.\n']);
fprintf("=======================================================================\n");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. LQR Control
%       Here we continiue using :
%           Trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating steady state value of output and input in discrete case (no feedback)
x_ss = (eye(n_states) - sysD.A)\sysD.B*v_ss;

%--------------------
%% Picking random diagonal value for Q and R
Q  = [4*eye(3), zeros(3,9); ...                % coordinates
       zeros(3,3), 1*eye(3), zeros(3,6); ...    % velocities
       zeros(3,6), 1*eye(3), zeros(3,3); ...    % angles
       zeros(3,9), 1*eye(3)] ;                  % angular velocity
R = 1*eye(n_inputs);

Q = 1*Q;
R = 100*R;
%% Calculating feedback gain for continuous and discrete cases
% lqr solver for feedback in nonlinear system
[K_c,~, CLP_c] = lqr(sys.A, sys.B, Q,R);

% Performing matlab dlqr solver
[K_d,~,CLP_d] = dlqr(sysD.A,sysD.B,Q,R);

%% calculating Nu & Nx in continuous and discrete case
% (n_states + n_outputs) * ( n_states + n_inputs)
big_A = [sys.A , sys.B;
         sys.C(1:3,:), sys.D(1:3,:)];
big_Y =[zeros(n_states,n_outputs - 3);
         eye(n_outputs - 3) ];
big_N = big_A\big_Y;

fprintf('\nReference-Input  full state feedback continuous case. The matrices Nx and Nu are:\n');
Nx_c = big_N(1:n_states,:)
Nu_c = big_N(n_states+1:end,:)

%----------------------------------------------
% (n_states + n_outputs) * ( n_states + n_inputs)
big_A = [sysD.A - eye(n_states), sysD.B;
         sysD.C(1:3,:), sysD.D(1:3,:)];
big_Y =[ zeros(n_states,n_outputs - 3);
         eye(n_outputs - 3) ];
big_N = big_A\big_Y;


fprintf('\nReference-Input  full state feedback discrete case. The matrices Nx and Nu are:\n');
Nx_short = big_N(1:n_states,:)
Nu_short = big_N(n_states+1:end,:)


%% Delete this part later
big_A = [sysD.A - eye(n_states), sysD.B;
         sysD.C, sysD.D];
big_Y =[ zeros(n_states,n_outputs);
         eye(n_outputs) ];
big_N = big_A\big_Y;


fprintf('\nReference-Input  full state feedback discrete case. The matrices Nx and Nu are:\n');
Nx_d = big_N(1:n_states,:)
Nu_d = big_N(n_states+1:end,:)

%% poles open loop vs closed loop
Open = eig(sysD.A)
Closed_discrete = eig(sysD.A - sysD.B*K_d)
Closed_continuous = eig(sys.A - sys.B*K_c)



%% Saving plot for LaTeX
%print -depsc Bifurcation1.eps


%% Comments

% I assume we're supposed to guess these values initially
% as they represent the relative importance of the speed of the controller
% vs the cost of the input
% So once this is implemented correctly we shoudl do trial and error till
% we're satisfied.




%% Comments

% feedback gain seems to be calculated correctly, however simulinkg
% implementation doesn't work.

%{
Also, i noticed that as references we only get the coordinates
However they only represent 3 of the 6 output values. 
That is why i've appended 3 zeros to the reference input.
I'm not sure if were supposed to do this or recalculate the gain, 
Nu and NX such that only 3 inputs are important
%}


%% steady state output guess
t1 = eye(n_states) - sysD.A;
t1 = t1(1:9,:);
t2 = sysD.B*(ones(4,1)*v_ss);
t2 = t2(1:9,:);
xrr = t1\t2;
yres = sysD.C*xrr + sysD.D*(ones(4,1)*v_ss);

%% guess two
%sym z
G =  tf(sysD)
tfmat=G(,3)
% tf (1-1/z) --> evalueren in z = 1 --> maal v_ss moet steady state error
% yss geven
%transfermat(3)


%% discretized system seems to have 12 poles and 8 transmission zeros

res = tzero(sysD);
tz = res(3)
M_T = [tz*eye(length(sysD.A))-sysD.A, -sysD.B; sysD.C, sysD.D];
rank(M_T)
rank([1.5*eye(length(sysD.A))-sysD.A, -sysD.B; sysD.C, sysD.D])
sizz = size(null(M_T))
F = null(M_T);
F = F(:,2);

Rx = F(1:n_states);
Ru = F(n_states + 1: end)

% if they are transmission zeros the following expressions should be close
% to 0
eq1 = sysD.A*Rx + sysD.B*Ru + tz*eye(n_states)*Rx
eq2 = sysD.C*Rx + sysD.D*Ru

