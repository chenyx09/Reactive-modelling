function phi_negative = two_robot_staliro_shortest_time_reactive(filter)
global app bdry Ts x_norm
params
if nargin<1
    filter=0;
end

% r2_con = @(pos,t1,t2)robot_con_w(pos,t1,t2,2);







input_range = [-vmax,vmax;-vmax,vmax;-vmax,vmax;-vmax,vmax];

cp_array = [10;10;10;10];



% phi = '!<>col';


dis_max = 1.5;

n_plane = 20;
[A,b]=inner_approx_circle(dis_max,n_plane,[0;0]);
state_const.P=A;
state_const.q=b;

phi = '[]connect';
preds(1).str='connect';
preds(1).A = [zeros(n_plane,2) A -A];
preds(1).b = b;

r2_con = @(r1,r2,t2)robot_con_w(r2,t2,r1,2);
% r2_con = @(r1,r2,t2)robot_con_safety(r1,r2,t2,state_const,2);
model = @(t,x,d)two_robot_sim(t,x,r2_con,d,filter);




%
% disp(' ')
% disp('Total Simulation time:')
% time = 10


opt = staliro_options();
opt.ode_solver = @ode23;

% disp(' Press any key to continue ... ')
% pause

opt.optimization_solver = 'CE_Taliro';
opt.runs = 1;

opt.optim_params.n_tests = 10;
opt.optim_params.num_subdivs = 10;
opt.optim_params.num_iteration = 10;
opt.spec_space='X';

% opt.n_workers = input('How many workers (concurrent simulations)? ');
opt.n_workers = 1;
opt.dispinfo = 0;
% opt

% time=5;
%
%     delta_pos = randn(2,1);
%     if norm(delta_pos)>0.8
%         delta_pos = delta_pos/norm(delta_pos)*0.8;
%     end
%     delta_pos = delta_pos * dis_max;
%     r1_pos0 = [bdry(1,1)+dis_max+(bdry(1,2)-bdry(1,1)-2*dis_max)*rand;bdry(1,1)+dis_max+(bdry(2,2)-bdry(2,1)-2*dis_max)*rand];
%     r2_pos0 = r1_pos0 + delta_pos;
%
% init_cond = [-4,4;-4,4;r1_pos0(1),r1_pos0(1);r1_pos0(2),r1_pos0(2);r2_pos0(1),r2_pos0(1);r2_pos0(2),r2_pos0(2)];
% [results, history] = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);
% [T,X,Y,IT] = SimFunctionMdl(model,init_cond,input_range,cp_array,results.run(results.optRobIndex).bestSample,time,opt);
% T_bar=T(1):Ts:T(end);
%         X_bar = interp1(T,X,T_bar);
%         U_bar=get_input_2_robot(T_bar,X_bar,r2_con);




t_max = 10;
t_min = 0.5;
tol =0.2;
err = 5;

delta_t = 0.1;
feas = 0;
while err>tol
    time = 0.5*(t_max+t_min);
    r1_pos0 = [bdry(1,1)+dis_max+(bdry(1,2)-bdry(1,1)-2*dis_max)*rand;bdry(1,1)+dis_max+(bdry(2,2)-bdry(2,1)-2*dis_max)*rand];
    delta_pos = randn(2,1);
    if norm(delta_pos)>0.8
        delta_pos = delta_pos/norm(delta_pos)*0.8;
    end
    delta_pos = delta_pos * dis_max;
    r2_pos0 = r1_pos0 + delta_pos;
    init_cond = [-4,4;-4,4;r1_pos0(1),r1_pos0(1);r1_pos0(2),r1_pos0(2);r2_pos0(1),r2_pos0(1);r2_pos0(2),r2_pos0(2)];
    
    % disp(' ')
    % disp('Ready to run S-TaLiRo with Cross Entropy ...')
    % disp(' Press any key to continue ... ')
    
    [results, history] = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);
    if results.run(results.optRobIndex).bestRob<-1e-4
        t_max = time;
        feas=1;
    else
        t_min = time;
    end
    err = abs(t_max-t_min);
    
end
time1 = time+delta_t
N_trace = 10;
Trace = [];
i=1;
if feas
    counter=0;
while i<=N_trace
    
    r1_pos0 = [bdry(1,1)+dis_max+(bdry(1,2)-bdry(1,1)-2*dis_max)*rand;bdry(1,1)+dis_max+(bdry(2,2)-bdry(2,1)-2*dis_max)*rand];
    delta_pos = randn(2,1);
    if norm(delta_pos)>0.8
        delta_pos = delta_pos/norm(delta_pos)*0.8;
    end
    delta_pos = delta_pos * dis_max;
    r2_pos0 = r1_pos0 + delta_pos;
    init_cond = [-4,4;-4,4;r1_pos0(1),r1_pos0(1);r1_pos0(2),r1_pos0(2);r2_pos0(1),r2_pos0(1);r2_pos0(2),r2_pos0(2)];
    [results, history] = staliro(model,init_cond,input_range,cp_array,phi,preds,time1,opt);
    
    if results.run(results.optRobIndex).bestRob<0
        
        [T,X,Y,IT] = SimFunctionMdl(model,init_cond,input_range,cp_array,results.run(results.optRobIndex).bestSample,time,opt);
        T_bar=T(1):Ts:T(end);
        X_bar = interp1(T,X,T_bar);
        U1_bar = interp1(IT(:,1),IT(:,4:5),T_bar);
        %         U_bar=get_input_2_robot(T_bar,X_bar,r2_con);
        h = min(kron(ones(1,length(T_bar)),preds(1).b)- preds(1).A*X_bar');
        vio_idx = min(find(h<0));
        if ~isempty(vio_idx)
            crit_idx = max(1,vio_idx-5):min(length(T_bar),vio_idx+5);
        else
            crit_idx=[];
        end
        X1 = X_bar(crit_idx,3:4);
        X2 = X_bar(crit_idx,5:6);
        Trace=[Trace;X1 X2 U1_bar(crit_idx,:)];
        i = i+1;
    end
    counter=counter+1;
    if counter>10*N_trace
        break
    end
end

feature_order=2;
p1=mss_asd(size(Trace,2)-1,feature_order);
n_p = size(p1,1);
p = [zeros(n_p,2) p1;
    ones(n_p,1) zeros(n_p,1) p1;
    zeros(n_p,1) ones(n_p,1) p1];
x_neg = Trace./kron(x_norm',ones(size(Trace,1),1));

M=size(p,1);
phi_negative=zeros(size(x_neg,1),M);
for i=1:size(x_neg,1)
    counter=1;
    for j=1:M
        %         if ~all(p(j,2:end)==0)
        phi_negative(i,counter)=prod([x_neg(i,5:6) x_neg(i,1:4) 1].^p(j,:));
        counter=counter+1;
        %         end
    end
    %     phi_negative(i,M+1:M+2)=x_neg(i,5:6);
    %     phi_positive(i,M+1)=SVM_data_positive(i,6);
    %     phi_positive(i,M+2)=SVM_data_positive(i,7);
end
else
    phi_negative=[];
end

end

function dxdt = two_robot_sim(t,x,r2_con,d,filter)
global bdry vmax u_poly
if nargin<5
    filter=0;
end
t2_pos = x(1:2);
r1_pos = x(3:4);
r2_pos = x(5:6);


dxdt      = zeros(6,1);

t2_vel = d(1:2);

if filter
    [Pu1_rea,qu1_rea] = robot_reactive_bound([r1_pos;r2_pos]);
    Pu1 = [Pu1_rea;u_poly.P];
    qu1 = [qu1_rea;u_poly.q];
%     opts1=  optimset('display','off');
%     r1_vel = quadprog(eye(2),-d(3:4),Pu1,qu1,[],[],[],[],[],opts1);
    options = qpOASES_options('printLevel',0);
%     options.printLevel = PL_LOW;
    [sol,fval,exitflag] =qpOASES( eye(2),-d(3:4),Pu1,[],[],[],qu1,options);
    if exitflag==0
    r1_vel = sol;
    else
        r1_vel = d(3:4);
    end
%     params.P = [Pu1;zeros(50-length(qu1),2)];
%     params.q = [qu1;ones(50-length(qu1),1)];
%     settings.verbose=0;
%     params.H = eye(2);
%     params.f = -2*d(3:4);
%     params.ub = [vmax;vmax];
%     params.lb = -[vmax;vmax];
%     [vars,status]=csolve(params,settings);
%     
%     
%     if strcmp(status,'Solved')
%     r1_vel = vars.x;
%     else
%         r1_vel = d(3:4);
%     end
else
    r1_vel = d(3:4);
end




if t2_pos(1)<bdry(1,1)
    t2_vel(1)=1;
elseif t2_pos(1)>bdry(1,2)
    t2_vel(1)=-1;
end


if t2_pos(2)<bdry(2,1)
    t2_vel(2)=1;
elseif t2_pos(2)>bdry(2,2)
    t2_vel(2)=-1;
end


if norm(r1_pos-r2_pos)>1.5
    disp('')
end
if t>2.8 &&t<2.9
    disp('')
end


r2_vel = r2_con(r1_pos,r2_pos,t2_pos);

if norm(r1_vel)>vmax
    r1_vel = r1_vel/norm(r1_vel)*vmax;
end

if norm(t2_vel)>vmax
    t2_vel = t2_vel/norm(r1_vel)*vmax;
end

dxdt = [t2_vel;r1_vel;r2_vel];
% dxdt (5:6)=zeros(2,1);

end

function vel = robot_con_w(pos,t1,t2,weight_choice)
if nargin<4
    weight_choice = 1;
end
global priority vmax
switch weight_choice
    case 1
        if pos(1)<=0 && pos(2)>0
            weight = priority(1,:);
        elseif pos(1)>0&&pos(2)>0
            weight = priority(2,:);
        elseif pos(1)<=0&&pos(2)<=0
            weight = priority(3,:);
        elseif pos(1)>0&&pos(2)<=0
            weight = priority(4,:);
        end
    case 2
        weight = [2 8];
    case 3
        weight = [10 0];
    case 4
        weight = [0 10];
end
max_dis = 5;
delta_x1 = t1-pos;
delta_x2 = t2-pos;
delta_x1 = delta_x1/norm(delta_x1)*min(max_dis,norm(delta_x1));
delta_x2 = delta_x2/norm(delta_x2)*min(max_dis,norm(delta_x2));
vel = (weight(1)*delta_x1 + weight(2)*delta_x2);
if norm(vel)>vmax
    vel = vel/norm(vel)*vmax;
end
end

function u2 = robot_con_safety(r1,r2,t2,const,weight_choice)
if nargin<6
    weight_choice = 1;
end
global vmax u_poly
u2=robot_con_w(r2,t2,r1,weight_choice);
if norm(u2)>vmax
    u2 = u2/norm(u2)*vmax;
end

Px=const.P;
qx=const.q;
h = qx-Px*(r1-r2);
[h_min,idx]=min(h);
if h_min<0
    disp('')
end


[Pu1_rea,qu1_rea] = const_para([r1;r2]);
if norm(Pu1_rea)>0.1
    disp('')
end
Pu1 = [Pu1_rea;u_poly.P];
qu1 = [qu1_rea;u_poly.q];
Pu1 = [Pu1;zeros(50-length(qu1),2)];
qu1 = [qu1;ones(50-length(qu1),1)];
params.H = zeros(2,2);
params.f = -Px(idx,:)';
params.P = Pu1;
params.q = qu1;
params.ub = [vmax;vmax];
params.lb = -[vmax;vmax];
settings.verbose=0;
[vars,status]=csolve(params,settings);
u1_worst=vars.x;

params.H = eye(2);
params.f = -2*u2;
Pu2 = [-Px;u_poly.P];
qu2 = [h-Px*u1_worst;u_poly.q];
params.P = [Pu2;zeros(50-length(qu2),2)];
params.q = [qu2;ones(50-length(qu2),1)];
[vars,status]=csolve(params,settings);
u2 = vars.x;



end


function U=get_input_2_robot(T,X,con)
t2_pos = X(:,1:2);
r1_pos = X(:,3:4);
r2_pos = X(:,5:6);
for i=1:length(T)
    U(i,:)=con(r1_pos(i,:)',r2_pos(i,:)',t2_pos(i,:)')';
end
end

