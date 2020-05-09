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

Ts = 0.05;


Ad = (eye(n_states) - A*Ts/2)\(eye(n_states) + A*Ts/2);
Bd = (eye(n_states) - A*Ts/2)\B*Ts;
Cd = C/(eye(n_states) - A*Ts/2);
Dd = D + C/(eye(n_states) - A*Ts/2) * B*Ts/2;
sysD = ss(Ad, Bd, Cd, Dd, Ts);


%% applying final value theorem
tff = tf(sysD);

syms z;

%-----------
G2 = tff(2,1);

f2T = (5.677e-06*z^4 + 2.271e-05*z^3 + 3.406e-05*z^2 + 2.271e-05*z + 5.677e-06)*(z-1)/z;
f2N = z^4 - 3.975*z^3 + 5.926*z^2 - 3.926*z + 0.9753;

f2Td = diff(f2T);
f2Nd = diff(f2N);

f2 = subs(f2Td, 1)/subs(f2Nd, 1);
% f2 = 0.090834000000000
%------------
G3 = tff(3,1);

f3T = (3.704e-05*z^2 + 7.407e-05*z + 3.704e-05)*(z-1)/z;
f3N = z^2 - 1.975*z + 0.9753;

f3Td = diff(f3T);
f3Nd = diff(f3N);

f3 = subs(f3Td, 1)/subs(f3Nd, 1);
% f3 = 0.005926000000000
%----
G4 = tff(4,1);
f4T = (0.0009375*z^2 + 0.001875*z + 0.0009375)*(z-1)/z;
f4N = z^2 - 2*z + 1;

f4Td = diff(f4T);
f4Nd = diff(f4N);

f4Tdd = diff(f4Td);
f4Ndd = diff(f4Nd);

f4 = subs(f4Tdd, 1)/subs(f4Ndd, 1);
%---------------

G6 = tff(6,1);
f6T = (6.25e-05*z^2 + 0.000125*z + 6.25e-05)*(z-1)/z;
f6N = z^2 - 2*z + 1;

f6Td = diff(f6T);
f6Nd = diff(f6N);

f6Tdd = diff(f6Td);
f6Ndd = diff(f6Nd);

f6 = subs(f6Tdd, 1)/subs(f6Ndd, 1);