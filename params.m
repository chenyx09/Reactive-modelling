global a_min a_max rmin rmax v0 vmax vmin psimax Ts N v_des mpc a_norm x_norm x_min ...
bdry priority p feature_order feature_dim x_idx d_idx TIMEOUTE_LIMIT EPS u_poly cost_type slack_max 

switch app
    case 1
%% Lane change
load LC_data
TIMEOUTE_LIMIT = 60;
EPS            = 1e-5;
cost_type = 'mixed';
slack_max = -1e-2;
% load LC_data
vy_max = 1.5;
vy_min = -1.5;
ax_max = 3;
ax_min = -5;

v0=32;
vmax = 40;
vmin = 0;
psimax = 0.1;

Ts=0.1;
% x_min=[-2;-1;-5;-0.0873];
% x_max=[15;4;5;0.0873];
% x_norm=x_max-x_min;
a_max=4;
a_min=-4;
a_norm = a_max-a_min;
rmax=0.1;
rmin=-0.1;
rnorm=rmax-rmin;



N = 15;
v_des = 35;
% clear model
% model = LTISystem('A', Ad, 'B', Bd);
% model.x.min = [-inf; x_min(4)-0.1];
% model.x.max = [inf; x_max(4)+0.1];
% model.u.min = [rmin];
% model.u.max = [rmax];
% 
% model.x.with('reference');
% model.x.reference = 'free';
% 
% model.x.penalty = OneNormFunction( diag([1 1]) );
% model.u.penalty = OneNormFunction( 1 );
% 
% ctrl = MPCController(model, N);

LC_mpc;
feature_order = p(1,end);
M = size(p,1);
for i=1:4
    x_idx(i) = find(p(:,i)==1&p(:,end)==feature_order-1);
end
d_idx = M+1;
n_class = 2;
n_plane = 1;

    case 2
%% 2 robot surveillance
TIMEOUTE_LIMIT = 60;
% CEGAR_SIZE     = 200;
% CEGAR_threshold = 1e-3;
EPS            = 1e-5;
bdry     = [-5,5;-5,5];
vmax     = 2;
cost_type = 'mixed';
slack_max = -1e-2;
[P,q]=inner_approx_circle(vmax,20,[0;0]);
u_poly.P=P;
u_poly.q=q;
priority = [0.2 0.8;0.8 0.2;0.5 0.5;0.5 0.5];
Ts = 0.1;
feature_dim = 6;
x_norm = [abs(bdry(:));vmax;vmax];
feature_order=2;
p1=mss_asd(feature_dim-1,feature_order);
n_p = size(p1,1);
p = [zeros(n_p,2) p1;
     ones(n_p,1) zeros(n_p,1) p1;
     zeros(n_p,1) ones(n_p,1) p1];
 
x1_row = p(1,:)*0;
x1_row(3)=1;
x1_row(end)=feature_order-1;
[tf, ~]=ismember(p,x1_row,'rows');
x1_idx = find(tf);

y1_row = p(1,:)*0;
y1_row(4)=1;
y1_row(end)=feature_order-1;
[tf, ~]=ismember(p,y1_row,'rows');
y1_idx = find(tf);

x2_row = p(1,:)*0;
x2_row(5)=1;
x2_row(end)=feature_order-1;
[tf, ~]=ismember(p,x2_row,'rows');
x2_idx = find(tf);

y2_row = p(1,:)*0;
y2_row(6)=1;
y2_row(end)=feature_order-1;
[tf, ~]=ismember(p,y2_row,'rows');
y2_idx = find(tf);
x_idx = [x1_idx;y1_idx;x2_idx;y2_idx];

u1x_row = p(1,:)*0;
u1x_row(1)=1;
u1x_row(end)=feature_order;
[tf, ~]=ismember(p,u1x_row,'rows');
u1x_idx = find(tf);

u1y_row = p(1,:)*0;
u1y_row(2)=1;
u1y_row(end)=feature_order;
[tf, ~]=ismember(p,u1y_row,'rows');
u1y_idx = find(tf);
n_class = 2;
n_plane = 3;
d_idx=[u1x_idx;u1y_idx];
end