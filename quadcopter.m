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

C(1,4) = 1;
C(2,5) = 1;
C(3,6) = 1;
C(4,10) = 1;
C(5,11) = 1;
C(6,12) = 1;

% NO feedthrough matrix
D = zeros(6,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating system
sys= ss(A,B,C,D);
G = tf(sys);
%Sampling time
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

% step and frequency responses
%{
%% Step response
T_sim = 5;
t_sim_steps = [0:Ts:T_sim];
step(G,G_DT_bilinear,G_DT_zoh,G_DT_FE,t_sim_steps)
%step(G,G_DT_bilinear,G_DT_zoh,Gpz_no_gain_correction,t)
legend show

%% Frequency response
figure;
H=bodeplot(G,G_DT_bilinear,G_DT_zoh, G_DT_FE);
%H=bodeplot(G,G_DT_bilinear,G_DT_zoh,Gpz_no_gain_correction);
legend show
grid

pp = getoptions(H);
pp.YLim{1}=[-25 25];
setoptions(H,pp)
%}

%% Performing checks on FE discretisation
fprintf("========================================================================");
fprintf(['\n\nThe original continiuous time system itself is unstable (pole = 0). \n', ...
     'The poles of the original matrix are :\n']);
disp(eig(A))
fprintf('\n\nThe poles of the discrete system using FE transformation are: \n');
disp(eig(A_FE))
fprintf('\nNotice the same unstable poles.\n');
fprintf("=======================================================================");
fprintf('\nNext we check if the discretized system is controllable and observable.');
fprintf(['\n\nWe can see that the rank of the controllability matrix is the same', ...
    'as the order of the system: \n']);
fprintf('\nOrder of the system is : %i \n', n_states);
fprintf('Rank of controllability matrix is : %i \n', rank(ctrb(A_FE, B_FE)));
fprintf('Rank of observability matrix is : %i \n', rank(obsv(A_FE, C_FE)));
fprintf(['\nThe system is thus controllable and not observable', ...
    'the lack of observability is because we do not measure the angular and lateral velocity vectors.']);
% This doesn't fully explain it ( rank 8 out of 12, but we drop 6 inputs)
fprintf('\n(A_FE, B_FE) is stabilizable, because all unstable nodes are controllable.\n');
% checking which modes are NOT observable using PBH Test
[V_FE,d_FE] = eig(A_FE);
result_m_FE = C_FE*V_FE;
fprintf('\nThe first 3 modes associated with the location are NOT observable and the yaw is also not observable.\n'); 
fprintf('\nThe system is thus NOT detectable!\n');
fprintf('\nAs the system is not detectable it is also NOT minimal!\n');
% Checking transmission zero's of discretized model
tzero(A_FE, B_FE, C_FE, D_FE)