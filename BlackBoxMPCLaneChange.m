function [T, XT, YT, LT, CLG, Guards] = BlackBoxMPCLaneChange(X0,simT,TU,U,filter)

LT = [];
CLG = [];
Guards = [];



% Change the parameter values in the model
% set_param([model,'/Pedal Angle (deg)'],'Amplitude',num2str(X0(1)));
% set_param([model,'/Pedal Angle (deg)'],'Period',num2str(X0(2)));

% Run the model
global N Ts mpc
if isempty(Ts)||isempty(mpc)||isempty(N)
    app = 1;
    params
end
NN = ceil(simT/Ts);
T = [];
XT = [];
YT = [];
x = X0;
for i=0:NN
%% use ODE solver, higher accuracy    
%     d = @(t)interp1(TU,U,t+i*Ts);
%     u = LC_MPC_con(x,d(0));
%     f = @(t,x)lane_change_dyn(t,x,d,u);
%     [tt,xx] = ode23(f,[0,Ts],x);
%     x = interp1(tt,xx,Ts)';
%     T = [T;tt+i*Ts];
%     XT = [XT;xx];
%     YT = [YT;xx repmat(u',length(tt),1) d(tt)];

%% use discrete integration, faster
    d = @(t)interp1(TU,U,i*Ts);
    u = LC_MPC_con(x,d(0));
    dxdt = lane_change_dyn(0,x,d,u,filter);
    T = [T;i*Ts];
    XT = [XT;x'];
    YT = [YT;x' u' d(0)];
    x = x + dxdt*Ts;
    
    
end

end

function dxdt = lane_change_dyn(t,x,d,u,filter)
a_min = -8;
dd = d(t);
if filter
    d_max = LC_reactive_bound(x([1,2,4,5]));
    if d_max<a_min
        d_max = a_min;
    end
    dd = min(d_max,dd);
end
    
dxdt = [x(4);...
    x(3)*sin(x(5));...
    u(1);...
    u(1)-dd;...
    u(2)];
end


function u = LC_MPC_con(x,d)
% x=[Delta X,Delta Y,v1,Delta v,psi]

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
% mpc.A = sparse([zeros(N,m*N) eye(N) eye(N) zeros(N,1)]);
mpc.obj = f;
mpc.Q=sparse(H);
mpc.rhs = bineq;
% mpc.rhs = ones(N,1);
% mpc.sense = repmat('<',1,length(mpc.rhs));
param.outputflag = 0;

result = gurobi(mpc, param);
u=result.x(1:m);
end
