function phi_negative = lane_change_staliro(isfilter)
if nargin<1
    isfilter = 0;
end
global app
params

% controller = @LC_MPC_con;
model = @(X0,simT,TU,U)BlackBoxMPCLaneChange(X0,simT,TU,U,isfilter);




init_cond = [-2,3;3.8,4.4;28,32;-3,3;-0.02,0.02];


input_range = [3,a_max];
% input_range = [a_max-0.5,a_max];
cp_array = 5;

disp(' ')
disp('The specification:')

time = 6;

% phi = '!<>col';
phi = ['([]!col)/\(<>_[0,',num2str(time),'] c)/\([]k)'];

min_Delta_X = 3.5;
min_Delta_Y = 3;


preds(1).str='col';
preds(1).A = [1 0 0 0 0;-1 0 0 0 0;0 1 0 0 0;0 -1 0 0 0];
preds(1).b = [min_Delta_X;min_Delta_X;min_Delta_Y;min_Delta_Y];

preds(2).str='c';
preds(2).A = [0 1 0 0 0;0 -1 0 0 0];
preds(2).b = [0.4;0.4];

preds(3).str='k';
preds(3).A = [0 -1 0 0 0];
preds(3).b = [0.9];



opt = staliro_options();

% disp(' Press any key to continue ... ')
% pause

opt.optimization_solver = 'CE_Taliro';
opt.runs = 1;

opt.optim_params.n_tests = 100;
opt.optim_params.num_subdivs = 25;
opt.optim_params.num_iteration = 25;
opt.spec_space='X';
opt.black_box = 1;
opt.interpolationtype = {'pconst'};
opt.loc_traj = 'end';

opt.map2line = 1;
% opt.n_workers = input('How many workers (concurrent simulations)? ');
opt.n_workers = 4;


% disp(' ')
% disp('Ready to run S-TaLiRo with Cross Entropy ...')
% disp(' Press any key to continue ... ')

tic
[results, history] = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);
toc




% for i=1:length(results.run)
%     res(i)=results.run(i).bestRob;
% end
if results.run(results.optRobIndex).bestRob<0

% figure
[T1,XT1,YT1,IT1,LT1] = SimBlackBoxMdl(model,init_cond,input_range,cp_array,results.run(results.optRobIndex).bestSample(:,1),time,opt);
x_negative = [];
critical_idx = find(~(XT1(:,1)>8 & XT1(:,2)<0.4));

x_negative = [x_negative;YT1(critical_idx,[1,2,4,5,8])];

M=size(p,1);
phi_negative=zeros(size(x_negative,1),M);
for i=1:size(x_negative,1)

    for j=1:M
        phi_negative(i,j)=prod([(x_negative(i,1:4)-x_min')./x_norm' 1].^p(j,:));
    end
    phi_negative(i,M+1)=(x_negative(i,5)-a_min)/a_norm;
end
else
    phi_negative=[];
end
end
