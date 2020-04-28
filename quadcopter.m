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
p_FE = eig(A_FE)
p_T = pole(G_DT_bilinear)
p_ZOH = pole(G_DT_zoh)

fprintf('\nNotice the same number of unstable poles.\n');
fprintf('\nNext we check if the discretized system is controllable and observable.');
fprintf(['\n\nWe can see that the rank of the controllability matrix is the same', ...
    'as the order of the system: \n']);
fprintf('\nOrder of the system is : %i \n', n_states);
fprintf('Rank of controllability matrix is : %i \n', rank(ctrb(A_FE, B_FE)));
fprintf('Rank of observability matrix is : %i \n', rank(obsv(A_FE, C_FE)));
fprintf(['\nThe system is thus controllable and observable.']);
fprintf('\n(A_FE, B_FE) is stabilizable, because all unstable nodes are controllable.\n');
fprintf('Furthermore, (A_FE, C_FE) is detectable, because all unstable nodes are observable.\n');
% Confirming that all nodes are controllable and observable using PBH test
[V_FE,~] = eig(A_FE);
result_C_FE = C_FE*V_FE;
[V_FET,~] = eig(A_FE');
result_B_FE = B_FE'*V_FET;
% Checking transmission zero's of discretized model
% tzero calculates the solution for the matrix-vector product
% at the bottom of page 82 in the course. This does not guarantee a rank
% decrease of G(e) that's why we need to check it for each solution
fprintf('\n\nHere we show the solutions to tzero for the different systems\n');
Continuous_T_zero = tzero(sys)
FE_T_zero = tzero(ss(G_DT_FE))
bilinear_T_zero = tzero(ss(G_DT_bilinear))
zoh_T_zero = tzero(ss(G_DT_zoh))
fprintf(['\n\nTrapezoidal rule and ZOH seem to have solutions, however we need to check if the are actually', ...
    'transmission zeros: \n']);
M_T = [bilinear_T_zero(2)*eye(length(A_FE))-A_FE -B_FE; C_FE D_FE];
z_T = size(null(M_T))
M_ZOH1 = [zoh_T_zero(1)*eye(length(ss(G_DT_zoh).A))-ss(G_DT_zoh).A -ss(G_DT_zoh).B; ss(G_DT_zoh).C ss(G_DT_zoh).D];
z_ZOH1 = size(null(M_ZOH1))
M_ZOH2 = [zoh_T_zero(2)*eye(length(ss(G_DT_zoh).A))-ss(G_DT_zoh).A -ss(G_DT_zoh).B; ss(G_DT_zoh).C ss(G_DT_zoh).D];
z_ZOH2 = size(null(M_ZOH2))
fprintf('\n\nWe can thus confirm that the ZOH discretization has 2 transmission zeros.\n');

% Performing Kalman Decomposition
[Abarc,Bbarc,Cbarc,~,kc] = ctrbf(A,B,C); 
[Abaro,Bbaro,Cbaro,~,ko] = obsvf(A,B,C);
fprintf(['\n\nLastly as the system is controllable and observable we know that the system is minimal' , ...
    'This is also confirmed after performing the Kalman decomposition.\n']);
fprintf("=======================================================================\n");


