params;

r1_con = @(pos,t1,t2)robot_con_2by2(pos,t1,t2);
r2_con = @(pos,t1,t2)robot_con_w(pos,t1,t2,2);
f = @(t,x)two_robot_dyn(t,x,r2_con);

N_trial = 200;
show_plot = 0;
x_rec = [];

for n = 1:N_trial
    n
    t2_pos0=bdry(:,1)+(bdry(:,2)-bdry(:,1)).*rand(2,1);
    t2_vel0 = randn(2,1);
    t2_vel0 = t2_vel0/norm(t2_vel0)*vmax;
    r1_pos0=bdry(:,1)+(bdry(:,2)-bdry(:,1)).*rand(2,1);
%     r2_pos0=bdry(:,1)+(bdry(:,2)-bdry(:,1)).*rand(2,1);
    r2_pos0=r1_pos0+0.5*rand(2,1);
    r2_pos0 = max(bdry(:,1),r2_pos0);
    r2_pos0 = min(bdry(:,2),r2_pos0);
    x0 = [zeros(2,1);t2_pos0;t2_vel0;r1_pos0;r2_pos0];
    Ts = 0.1;
    t = 0:Ts:20;
    x_rec_trial = [];
    % [t,x_rec] = ode23(f,tspan,x0);
    x = x0;
    % x_rec = [x0' zeros(1,4)];
    for i=1:length(t)
        dxdt = two_robot_dyn(t(i),x,r2_con);
        deltav1 = x(1:2);
        t2_pos  = x(3:4);
        t2_vel = x(5:6);
        r1_pos = x(7:8);
        r2_pos = x(9:10);
        r1_vel = dxdt(7:8);
        r2_vel = dxdt(9:10);
        x = x + dxdt*Ts;
        
        if norm(r1_pos-r2_pos)>2
            disp('')
        end
        x_rec_trial = [x_rec_trial;x([3:4,7:10])' r1_vel' r2_vel'];
        
    end
    
    if show_plot
        figure(1)
        
        for i=1:length(t)
            clf
            hold on
            plot(x_rec_trial(i,1),x_rec_trial(i,2),'rO')
            plot(x_rec_trial(i,3),x_rec_trial(i,4),'bs')
            plot(x_rec_trial(i,5),x_rec_trial(i,6),'rs')
            plot(x_rec_trial(i,[3,5]),x_rec_trial(i,[4,6]));
            dim = [(0.5*sum(x_rec_trial(i,[3,5]))-bdry(1,1))/(bdry(1,2)-bdry(1,1))-0.05,(0.5*sum(x_rec_trial(i,[4,6]))-bdry(2,1))/(bdry(2,2)-bdry(2,1)),0.1,0.1];
            dim = min(ones(1,4),dim);
            dim=max(zeros(1,4),dim);
            annotation('textbox',dim,'String',num2str(norm(x_rec_trial(i,3:4)-x_rec_trial(i,5:6))),'linestyle','none')
            axis([bdry(1,1)-1,bdry(1,2)+1,bdry(2,1)-1,bdry(2,2)+1]);
            drawnow
            pause(0.05)
        end
    end
    x_rec = [x_rec;x_rec_trial];
end
% [r1_x,r1_y,r2_x,r2_y,r1_vx,r1_vy]
feature_data = [x_rec(:,[3:8])./kron(x_norm',ones(size(x_rec,1),1))];


M=size(p,1);
phi_positive=zeros(size(feature_data,1),M);
for i=1:size(feature_data,1)
    counter=1;
    for j=1:M
%         if ~all(p(j,2:end)==0)
            phi_positive(i,counter)=prod([feature_data(i,5:6) feature_data(i,1:4) 1 ].^p(j,:));
            counter=counter+1;
%         end
    end
%     phi_positive(i,M+1:M+2)=feature_data(i,5:6);
    %     phi_positive(i,M+1)=SVM_data_positive(i,6);
    %     phi_positive(i,M+2)=SVM_data_positive(i,7);
end


function dxdt = two_robot_dyn(t,x,r2_con)
%% x=[deltav1_x;deltav1_y;t2_x;t2_y;t2_vx;t2_vy;r1_x;r1_y;r2_x;r2_y]
global bdry vmax
dxdt      = zeros(10,1);
delta_acce1   = randn(2,1)*0.5;
t2_acce   = randn(2,1);

deltav1 = x(1:2);
t2_pos  = x(3:4);
t2_vel = x(5:6);
r1_pos = x(7:8);
r2_pos = x(9:10);

if norm(deltav1)>0.2*vmax
    delta_acce1=-deltav1;
end
if t2_pos(1)<bdry(1,1)
    t2_acce(1)=1;
elseif t2_pos(1)>bdry(1,2)
    t2_acce(1)=-1;
end

if t2_pos(2)<bdry(2,1)
    t2_acce(2)=1;
elseif t2_pos(2)>bdry(2,2)
    t2_acce(2)=-1;
end

if norm(t2_vel)>vmax
    t2_acce = t2_acce+t2_vel*(vmax-norm(t2_vel));
end
if r1_pos(1)<0
    r1_vel = [-0.4 1 ;-1 -0.4]*(r1_pos-r2_pos)+deltav1;
else
    r1_vel = [-0.4 -1 ;1 -0.4]*(r1_pos-r2_pos)+deltav1;
end
if norm(r1_vel)>vmax
    r1_vel = r1_vel/norm(r1_vel)*vmax;
end
r2_vel = r2_con(r2_pos,t2_pos,r1_pos);

dxdt = [delta_acce1;t2_vel;t2_acce;r1_vel;r2_vel];

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

function vel = robot_con_2by2(pos,t1,t2)

global vmax
weight = [1;9];
max_dis = 4;
delta_x1 = t1-pos;
delta_x2 = t2-pos;
delta_x1 = delta_x1/norm(delta_x1)*min(max_dis,norm(delta_x1));
delta_x2 = delta_x2/norm(delta_x2)*min(max_dis,norm(delta_x2));
% if pos(1)<=0 && pos(2)>0
%     vel = [0.3*vmax;0];
% elseif pos(1)>0 && pos(2)<=0
%     vel = [-0.3*vmax;0];
% else
%     vel = weight(1)*delta_x1 + weight(2)*delta_x2;
% end
if pos(1)<=0
    vel = [0.3*vmax;0];
elseif pos(1)>0
    vel = weight(1)*delta_x1 + weight(2)*delta_x2;
end

% vel = weight(1)*delta_x1 + weight(2)*delta_x2;
if norm(vel)>vmax
    vel = vel/norm(vel)*vmax;
    
end
end