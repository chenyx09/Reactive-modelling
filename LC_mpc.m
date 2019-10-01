% function u = LC_MPC_con(x,ar)
%%
x = [4;4;18;3;0];
x0 = x;
ar = 0;
global N Ts v_des mpc
%  x = [Delta_X,Y,v,v-vr,psi]
Delta_X_max = 4;
Delta_Y_max = 3.5;
bigM = 1e4;
%%
v0 = x(3);
vr = x(4);
Th = 1.4;
Ac = [0 0 0 1 0;...
    0 0 0 0 v0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0];
Bc = [0 0;0 0;1 0;1 0;0 1];
Ec = [0;0;0;-1;0];
n = size(Ac,1);
m = size(Bc,2);
[dyn.A,dyn.B] = c2d(Ac,Bc,Ts);
[~,dyn.E] = c2d(Ac,Ec,Ts);

dyn.K = zeros(5,1);
tau = 2;
d_pred = (exp(-(0:N-1)*Ts/tau)*ar)';
[Lx, Lu, Ld, Lk] = mpc_matrices(dyn,N);
Q0 = diag([0.2,2,0.3,0,1]);

R0 = diag([0.02,0.2]);
u_max = [a_max;0.2];
u_min = [a_min;-0.2];
Q_mpc = kron(eye(N),Q0);
R_mpc = kron(eye(N),R0);

x_des = [Th*vr;0;v_des;v_des;0];


k = Lx*x+Ld*d_pred;
% Lu_xf = Lu(5(N-1)+1:end,:);
% k_xf = k(5*(N-1)+1:end);

% H_xf = Lu_xf'*Q_xf*Lu_xf;
% f_xf = (k_xf-x_des)'*Q_xf*Lu_xf;


f_x = 2*(k-kron(ones(N,1),x_des))'*Q_mpc*Lu;
H_x = Lu'*Q_mpc*Lu+R_mpc;

H = blkdiag(H_x,zeros(2*N+1,2*N+1));
f = [f_x zeros(1,2*N) 1e8];

H = 0.5*(H+H');
P1 = [0 -1 0 0 0];
q1 = -Delta_Y_max;
P1_mpc = kron(eye(N),P1);
q1_mpc = kron(ones(N,1),q1);
P2 = [-1 0 0 0 0];
q2 = -Delta_X_max;
P2_mpc = kron(eye(N),P2);
q2_mpc = kron(ones(N,1),q2);

Aineq = [P1_mpc*Lu -eye(N)*bigM zeros(N,N) -ones(N,1);...
    P2_mpc*Lu zeros(N,N) -eye(N)*bigM -ones(N,1);...
    zeros(N,m*N) eye(N) eye(N) zeros(N,1)];
bineq = [q1_mpc-P1_mpc*k;...
    q2_mpc-P2_mpc*k;...
    ones(N,1)];

mpc.A = sparse(Aineq);
mpc.obj = f;
mpc.Q=sparse(H);
mpc.rhs = bineq;
mpc.sense = repmat('<',1,length(bineq));
mpc.vtype = [repmat('C',m*N,1);repmat('B',2*N,1);'C'];
mpc.modelsense = 'min';
mpc.lb=[kron(ones(1,N),u_min') zeros(1,2*N+1)];
mpc.ub=[kron(ones(1,N),u_max') ones(1,2*N) inf];
mpc.Q_mpc = Q_mpc;
mpc.R_mpc = R_mpc;
mpc.P1_mpc = P1_mpc;
mpc.q1_mpc = q1_mpc;
mpc.P2_mpc = P2_mpc;
mpc.q2_mpc = q2_mpc;
%     model.varnames = names;

gurobi_write(mpc, 'mip.lp');

% clear params;
% params.outputflag = 0;
% 
% result = gurobi(mpc, params);
% u=result.x(1:m*N);
% con = @LC_MPC_con;
% dinput = @(t)0;
% f = @(t,x)lane_change_dyn(t,x,dinput,con);
% [T,Y] = ode23(f,[0,6],x0);


function dxdt = lane_change_dyn(t,x,d,con)

u = con(x,d(t));
dxdt = [x(4);...
    x(3)*sin(x(5));...
    u(1);...
    u(1)-d(t);...
    u(2);
    0];
end

function u = LC_MPC_con(x,d)
global N Ts v_des mpc
v0 = x(3);
delta_v = x(4);
vr = v0-delta_v;
Th = 1.4;
bigM = 1e3;
Ac = [0 0 0 1 0;...
    0 0 0 0 v0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0];
Bc = [0 0;0 0;1 0;1 0;0 1];
Ec = [0;0;0;-1;0];
m = size(Bc,2);
[dyn.A,dyn.B] = c2d(Ac,Bc,Ts);
[~,dyn.E] = c2d(Ac,Ec,Ts);

dyn.K = zeros(5,1);
tau = 2;
d_pred = (exp(-(0:N-1)*Ts/tau)*d)';
[Lx, Lu, Ld, ~] = mpc_matrices(dyn,N);

x_des = [Th*vr;0;v_des;v_des;0];
k = Lx*x+Ld*d_pred;

f_x = 2*(k-kron(ones(N,1),x_des))'*mpc.Q_mpc*Lu;
H_x = Lu'*mpc.Q_mpc*Lu+mpc.R_mpc;

H = blkdiag(H_x,zeros(2*N+1,2*N+1));
f = [f_x zeros(1,2*N) 1e6];

H = 0.5*(H+H');


Aineq = [mpc.P1_mpc*Lu -eye(N)*bigM zeros(N,N) -ones(N,1);...
    mpc.P2_mpc*Lu zeros(N,N) -eye(N)*bigM -ones(N,1);...
    zeros(N,m*N) eye(N) eye(N) zeros(N,1)];
bineq = [mpc.q1_mpc-mpc.P1_mpc*k;...
    mpc.q2_mpc-mpc.P2_mpc*k;...
    ones(N,1)];
mpc.A = sparse(Aineq);
mpc.obj = f;
mpc.Q=sparse(H);
mpc.rhs = bineq;
param.outputflag = 0;

result = gurobi(mpc, param);
u=result.x(1:m);
end
