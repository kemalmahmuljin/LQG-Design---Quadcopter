%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Initialization
%   First run this before any .slx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%% Loading references and setting parameters
load references_02.mat
n_states = 12;
n_inputs = 4;
n_outputs = 6;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Linearization
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
A(5,7) = -k * c_m / m * 4 * v_ss;
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




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Discretization
%       Here we discretize using :
%           Trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating system object
sys= ss(A,B,C,D);
G = tf(sys);
% Sampling time
Ts = 0.05;

%%  manual trapezoidal
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
step(G,t_sim_steps)
legend show

%% Frequency response
figure;
H=bodeplot(G);
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
step(G(4,1),G_DT_bilinear(4,1),t_sim_steps)
legend('G41', 'G_DT_bilinear41');

%% Frequency response
figure;
H=bodeplot(G(4,1),G_DT_bilinear(4,1));
legend('G41', 'G_DT_bilinear41');
grid
%}

%% Performing checks on discretisation
fprintf("========================================================================");
fprintf(['\n\nThe original continiuous time system itself is unstable (pole = 0). \n', ...
     'The poles of the continiuous system are :\n']);
 
disp(eig(A))

fprintf('\n\nThe poles of the discretized system using any transformation are: \n');

p_D = eig(sysD.A)

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
% decrease of G(e) that's why we need to check it for each solution !
fprintf('\n\nHere we show the solutions to tzero for the different systems\n');
Continuous_T_zero = tzero(sys)
D_T_zero = tzero(sysD)
fprintf(['\n\nTrapezoidal rule and seems to have solutions, however we need to check if the are actually', ...
    'transmission zeros: \n']);
M_D = [D_T_zero(1)*eye(length(sysD.A))-sysD.A -sysD.B; sysD.C sysD.D];
z_D = size(null(M_D))
rank([5*eye(length(sysD.A))-sysD.A -sysD.B; sysD.C sysD.D])
fprintf('\n\nWe can thus confirm that the Trapezoidal discretization has 4 transmission zeros.\n');

% Performing Kalman Decomposition to confirm minimal realization
[ADc,BDc,CDc,~,kDc] = ctrbf(sysD.A,sysD.B,sysD.C); 
[ADo,BDo,CDo,~,kDo] = obsvf(sysD.A,sysD.B,sysD.C);

[Abarc,Bbarc,Cbarc,~,kc] = ctrbf(A,B,C); 
[Abaro,Bbaro,Cbaro,~,ko] = obsvf(A,B,C);
fprintf(['\n\nLastly as the system is controllable and observable we know that the system is minimal' , ...
    'This is also confirmed after performing the Kalman decomposition.\n']);
fprintf("=======================================================================\n");




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.a LQR Control: Full state Feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---
% Calculating Nu & Nx in discrete case 
%---

% (n_states + n_references) * ( n_states + n_inputs)
big_A = [ (sysD.A - eye(n_states)), sysD.B;
         sysD.C(1:3,:), sysD.D(1:3,:)];
big_Y =[ zeros(n_states,3);
         eye(3) ];
big_N = big_A\big_Y;

Nx_short = big_N(1:n_states,:);
Nu_short = big_N(n_states+1:end,:);

%---
% Picking value for Q and R
%---

Q= [    1, zeros(1,11); ...                     % x - coordinate
        0, 1, zeros(1,10); ...                  % y - coordinate
        0, 0, 10, zeros(1,9); ...               % z - coordinate
        zeros(1,3), 1, zeros(1,8); ...          % x - velocities
        zeros(1,4), 1, zeros(1,7); ...          % y - velocities
        zeros(1,5), 1, zeros(1,6); ...          % z - velocities        
        zeros(1,6), 1, zeros(1,5); ...          % phi - angles
        zeros(1,7), 1, zeros(1,4); ...          % theta - angles
        zeros(1,8), 1, zeros(1,3); ...          % psi - angles    
        zeros(3,9), 1*eye(3)] ;                 % angular velocity
R = 0.0001*eye(n_inputs);

% when adding a weight
QW = 0.01*eye(n_states);
QW(3,3) = 700;
RW = 0.009*eye(n_inputs);

%---
% Calculating feedback gain for continuous and discrete cases
%---

% Performing matlab dlqr solver
[K_d,~,CLP_d] = dlqr(sysD.A,sysD.B,Q,R);
[K_dW,~,CLP_dW] = dlqr(sysD.A,sysD.B,QW,RW);



% %% poles open loop vs closed loop
% Open = eig(sysD.A)
% Closed_discrete = eig(sysD.A - sysD.B*K_d)




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.b LQR Control: Integral Action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructing the Augmented system
NA = [ eye(3,3)  sysD.C(1:3,:) ;
       zeros(n_states,3)  sysD.A]

NB = [ sysD.D(1:3,:) ; 
       sysD.B]

%checking the controlabillity of the Augmented system
disp('Rank of the controllability matrix of the augmented system:');
rank(ctrb(NA,NB))

%%
% Computing Q and R for 0kg case
Q  = [  350, zeros(1,11); ...                     % x - coordinate
        0, 350, zeros(1,10); ...                  % y - coordinate
        0, 0, 300, zeros(1,9); ...              % z - coordinate
        zeros(1,3), 1, zeros(1,8); ...          % x - velocities
        zeros(1,4), 1, zeros(1,7); ...          % y - velocities
        zeros(1,5), 1, zeros(1,6); ...          % z - velocities         
        zeros(1,6), 300, zeros(1,5); ...        % phi - angles
        zeros(1,7), 300, zeros(1,4); ...        % theta - angles
        zeros(1,8), 1, zeros(1,3); ...          % psi - angles
        zeros(3,9), 1*eye(3)] ;                  % angular velocity
Q = [   4*eye(3,3) zeros(3,12);
        zeros(12,3) 1*Q] 
R = eye(n_inputs);
Q = 1*Q;
R = 0.04*R;

% Computing Q and R for 0.1kg case


% Performing matlab dlqr solver
[K_d_Int,S,CLP_d_Int] = dlqr(NA,NB,Q,R);

Ki = K_d_Int(:,1:3)
Ks = K_d_Int(:,3+1:end)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. LQG control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
Usefull matlab functions
    destim: [Ae,Be,Ce,De] = destim(A,B,C,D,L)
    kalman: [KEST,L,P] = kalman(SYS,QN,RN,NN)
    dlqe: [M,P,Z,E] = dlqe(A,G,C,Q,R)
%}

% Checking whether the system is observable or not
Ob = obsv(sysD.A,sysD.C);
disp('LQG: Rank of the observability matrix');
rank(Ob)

var_xyz = 2.5e-5;
var_angle = 7.57e-5;
% measurement noise covariance R:
R = diag([var_xyz, var_xyz, var_xyz, var_angle, var_angle, var_angle]);
% process noise covariance Q: 
Q = 1*eye(n_states); % video of prof shows that it should be taken to have this form
% Calculating gain
[M,P,Z,E] = dlqe(sysD.A,eye(n_states),sysD.C,Q,R);
LKalman = sysD.A*M;




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. State feedback with pole-placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Controller gain K_pp
% Determination of the desired closed-loop poles
% The maximum admitted error is given by the size of the spheres : 0.08m
% ! might be better to use overshoot as a criterium
% Mp = 12*8/5;                             % Allowing 6.4% overshoot
% xi = log(Mp)/sqrt(log(Mp)^2 + pi^2);
xi = 0.999;
ts = 5.8;                                 % Settling time 
wn = 4.6/(xi*ts);                       % Natural frequency
alpha = -xi*wn;                         % Real part of the dominant poles
beta = wn*sqrt(1-xi^2);                 % Imaginary part of the dominant poles

non_dom = linspace(3.4*alpha,3.6*alpha, 10);
C_list = [alpha + i*beta, alpha - i*beta, non_dom];

% mapping poles to discrete time
P_list = exp(C_list.*Ts);
Kpp = place(sysD.A,sysD.B,P_list);

% Estimator gain L_pp
% Becuase the system is observable we can arbitrarely place the poles of
% the estimator !
F1 = 1.5;
F2 = 1.5;
L_list = [F1*alpha + i*beta, F1*alpha - i*beta, F2*non_dom];
PL_list = exp(L_list.*Ts);

L = place(sysD.A', sysD.C', PL_list)';
% smaller poles lead to faster reaction time
% however they also introduce rapid changes in estimation error

%% Saving plot for LaTeX
%print -depsc ../plots/name.eps























