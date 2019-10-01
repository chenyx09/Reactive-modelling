
clear w_val
params
global gamma
gamma =100;

% x = load('robot_data');
% phi_p = x.phi_positive;
% phi_n = x.phi_negative;
phi_p = phi_positive;
phi_n = phi_negative;

% threshold_x = 0;
n_p = size(phi_p,1);
n_n = size(phi_n,1);








% [status, sol, obj] = cegar_SVM(phi_p_bar,phi_n_bar,n_class,n_plane);
threshold_x = 0;
old_threshold_x = threshold_x;
for i=1:20
    %     [status,sol,obj]=Cplex_SVM_soft_bdry(phi_p,phi_n,n_class,threshold_x);
    [w_val,slack,flag]=Cplex_SVM_soft_bdry_polytope(phi_p,phi_n,n_class,n_plane,threshold_x);
    % [status,sol,obj]=Cplex_SVM(phi_p,phi_n,n_class,n_plane,[]);
    if flag
        sol.threshold_x = threshold_x;
        
        deriv = 1;
        tol = 1e-3;
        counter=0;
        NN = 20;
        alpha = 0.1/gamma;
        x_p = phi_p(:,x_idx(1));
        x_n = phi_n(:,x_idx(1));
        if n_plane >1
            vp1 = min(w_val{1}*phi_p')';
            vp2 = min(w_val{2}*phi_p')';
            vn1 = min(w_val{1}*phi_n')';
            vn2 = min(w_val{2}*phi_n')';
        else
            vp1 = (w_val{1}*phi_p')';
            vp2 = (w_val{2}*phi_p')';
            vn1 = (w_val{1}*phi_n')';
            vn2 = (w_val{2}*phi_n')';
        end
        
        vp1 = max(0,vp1)*1e-6 + min(0,vp1)*10;
        vp2 = max(0,vp2)*1e-6 + min(0,vp2)*10;
        %         vn1 = min(slack_max,vn1) + max(slack_max,vn1)*100;
        %         vn2 = min(slack_max,vn2) + max(slack_max,vn2)*100;
        while counter <NN &&abs(deriv)>tol
            
            sigma_p1 = 1./(1+exp(gamma*(x_p-threshold_x)));
            sigma_p2 = 1-sigma_p1;
            sigma_n1 = 1./(1+exp(gamma*(x_n-threshold_x)));
            sigma_n2 = 1-sigma_n1;
            if n_plane>1
                vp = min(w_val{1}*phi_p'.*sigma_p1' + w_val{2}*phi_p'.*sigma_p2')';
                vn = min(w_val{1}*phi_n'.*sigma_n1' + w_val{2}*phi_n'.*sigma_n2')';
            else
                vp = (w_val{1}*phi_p'.*sigma_p1' + w_val{2}*phi_p'.*sigma_p2')';
                vn = (w_val{1}*phi_n'.*sigma_n1' + w_val{2}*phi_n'.*sigma_n2')';
            end
            %         sigma_n1 = 1./(1+exp(gamma*(x_n-threshold_x)));
            %         sigma_n2 = 1-sigma_n1;
            grad_p = gamma*exp(gamma*(x_p-threshold_x))./((1+exp(gamma*(x_p-threshold_x))).^2);
            grad_n = gamma*exp(gamma*(x_n-threshold_x))./((1+exp(gamma*(x_n-threshold_x))).^2);
            deriv =  sum(grad_p.*(vp1-vp2))*1e-2 - sum(grad_n.*(vn1-vn2));
            threshold_x = threshold_x+deriv*alpha;
            counter=counter+1;
        end
        threshold_x = max(-1.1,threshold_x);
        threshold_x = min(1.1,threshold_x);
        if norm(old_threshold_x-threshold_x)<2e-2
            break
        else
            old_threshold_x = threshold_x;
        end
    else
        threshold_x = -1+2*rand;
        old_threshold_x = threshold_x;
        
    end
    disp(['threshold_x = ',num2str(threshold_x)]);
end
extract_reactive_bound_2region_soft_bdry(w_val,gamma,threshold_x)

pos_mis = [];
neg_mis = [];
switch app
    case 1
        for i=1:n_p
            x = phi_p(i,x_idx)'.*x_norm + x_min;
            d = phi_p(i,d_idx)'.*a_norm + a_min;
            d_ub=LC_reactive_bound(x);
            if ~(d<=d_ub)
                pos_mis=[pos_mis;i];
            end
        end
        
        for i=1:n_n
            x = phi_n(i,x_idx)'.*x_norm + x_min;
            d = phi_n(i,d_idx)'.*a_norm + a_min;
            d_ub=LC_reactive_bound(x);
            if (d<=d_ub)
                neg_mis=[neg_mis;i];
            end
        end
    case 2
        for i=1:n_p
            x = phi_p(i,x_idx)'.*x_norm(1:4);
            d = phi_p(i,d_idx)'.*x_norm(5:6);
            [A,b]=const_para(x);
            if ~all(A*d<=b)
                pos_mis=[pos_mis;i];
            end
        end
        neg_mis = [];
        for i=1:n_n
            x = phi_n(i,x_idx)'.*x_norm(1:4);
            d = phi_n(i,d_idx)'.*x_norm(5:6);
            [A,b]=const_para(x);
            if all(A*d<=b)
                neg_mis=[neg_mis;i];
            end
        end
end
% disp(['number of positive mismatch = ',num2str(length(pos_mis))])
% disp(['number of negative mismatch = ',num2str(length(neg_mis))])
% if ~isempty(sol)
% extract_reactive_bound_2region(sol);
% else
%     disp('failure')
% end
% pos_mis = [];
% for i=1:n_p
%     x = phi_p(i,x_idx)'.*x_norm(1:4);
%     u = phi_p(i,u1_idx)'.*x_norm(5:6);
%     [A,b]=const_para(x);
%     if ~all(A*u<=b)
%         pos_mis=[pos_mis;i];
%     end
% end
% neg_mis = [];
% for i=1:n_n
%     x = phi_n(i,x_idx)'.*x_norm(1:4);
%     u = phi_n(i,u1_idx)'.*x_norm(5:6);
%     [A,b]=const_para(x);
%     if all(A*u<=b)
%         neg_mis=[neg_mis;i];
%     end
% end


function [w_val,slack,flag]=Cplex_SVM_soft_bdry_polytope(phi_p,phi_n,n_class,n_plane,threshold)
global x_idx gamma
flag = 0;
[n_n,d] = size(phi_n);
% n_p = size(phi_p,1);
slack = inf*ones(n_plane,n_n);
cut_threshold = round(0.2*n_n);
cut_slack = -0.01;
for n=1:n_plane
    if n==1
        [status,sol,obj]=Cplex_SVM_soft_bdry(phi_p,phi_n,n_class,threshold);
    else
        min_slack = min(slack);
        [sort_slack,sort_idx]=sort(min_slack);
        if sort_slack(end-cut_threshold)<cut_slack
            phi_n_sample = phi_n(sort_idx(end-cut_threshold:end),:);
        else
            phi_n_sample = phi_n(min_slack>=cut_slack,:);
        end
        [status,sol,obj]=Cplex_SVM_soft_bdry(phi_p,phi_n_sample,n_class,threshold);
    end
    if ~isempty(sol)
        flag = 1;
        x_n=phi_n(:,x_idx(1));
        sigma_n = 1./(1+exp(gamma*(x_n-threshold)));
        
        slack(n,:)=(phi_n*sol.w(:,1).*sigma_n+phi_n*sol.w(:,2).*(1-sigma_n))';
        w_val{1}(n,:)=sol.w(:,1)';
        w_val{2}(n,:)=sol.w(:,2)';
    else
        w_val{1}(n,:)=zeros(1,d);
        w_val{2}(n,:)=zeros(1,d);
        slack(n,:) = zeros(1,n_n);
        break
    end
end
end

function [status,sol,obj]=Cplex_SVM_soft_bdry(phi_p,phi_n,n_class,threshold)
global TIMEOUTE_LIMIT x_idx cost_type slack_max x_norm gamma


x_p=phi_p(:,x_idx(1));
sigma_p = 1./(1+exp(gamma*(x_p-threshold)));
% phi_p_bar = [sigma_p.*phi_p (1-sigma_p).*phi_p];

x_n=phi_n(:,x_idx(1));
sigma_n = 1./(1+exp(gamma*(x_n-threshold)));

d = size(phi_p,2);
n_p = size(phi_p,1);
n_n = size(phi_n,1);
numberOfCores = 4;

num.w = d*n_class;

% num.slack_p = n_p;

num.slack_np = n_n;
num.slack_nn = n_n;

num.slack_max = 1;

var_name = {'w', 'slack_np','slack_nn' ,'slack_max'};
nvar = 0;

for v = var_name
    var = v{1};
    ndx.(var) = nvar + (1:num.(var));
    nvar = nvar + num.(var);
end
% ndx.b_s_p = reshape(ndx.b_s_p,n_class,n_p);
% ndx.b_s_n = reshape(ndx.b_s_n,n_class,n_n);
ndx.w = reshape(ndx.w,d,n_class);
% ndx.slack = reshape(ndx.slack,n_class,n_n);

% options = cplexoptimset('Display', 'off');
cplex = Cplex('SVM');
cplex.DisplayFunc=[];
obj = zeros(nvar,1);
if strcmp(cost_type,'L1')
    obj(reshape(ndx.slack_nn,[],1))=ones(num.slack_nn,1);
    obj(reshape(ndx.slack_np,[],1))=100*ones(num.slack_np,1);
elseif strcmp(cost_type,'Linf')
    obj(ndx.slack_max)=5000;
elseif strcmp(cost_type,'mixed')
    obj(reshape(ndx.slack_nn,[],1))=ones(num.slack_nn,1);
    obj(reshape(ndx.slack_np,[],1))=50*ones(num.slack_np,1);
    obj(ndx.slack_max)=300*n_n;
else
    error('cost type not supported')
end
lb = -inf*ones(nvar,1);
ub = inf*ones(nvar,1);
lb(ndx.slack_np)=zeros(num.slack_np,1);
ub(ndx.slack_nn)=zeros(num.slack_nn,1);

% if strcmp(cost_type,'L1')
%     ub(ndx.slack_np)=slack_max*ones(num.slack_np,1);
% end
for i=1:n_class
    lb(ndx.w(1,i))=-10;
    ub(ndx.w(1,i))=10;
end

% lb(reshape(ndx.lambda,[],1))=zeros(num.lambda,1);
cplex.addCols(obj,[],lb,ub);
cplex.Param.threads.Cur = numberOfCores;
cplex.Param.parallel.Cur =1;
% cplex.Param.mip.display.Cur=0;
cplex.Param.timelimit.Cur = TIMEOUTE_LIMIT;
% cplex.Param.timelimit.Cur = max(TIMEOUTE_LIMIT,n_p/5);
cplex.Param.output.clonelog.Cur=0;

for i=1:n_class
    Q = sparse(nvar,nvar);
    Q(ndx.w(2:end,i),ndx.w(2:end,i))=eye(d-1);
    cplex.addQCs(sparse(nvar,1), Q, 'L', 1.0);
end


A = sparse(n_p,nvar);
A(:,ndx.w(:,1))=sigma_p.*phi_p;
A(:,ndx.w(:,2))=(1-sigma_p).*phi_p;
cplex.addRows(zeros(n_p,1),A,inf*ones(n_p,1));

A = sparse(n_n,nvar);
A(:,ndx.w(:,1))=sigma_n.*phi_n;
A(:,ndx.w(:,2))=(1-sigma_n).*phi_n;
A(:,ndx.slack_np) = -speye(n_n);
A(:,ndx.slack_nn) = -speye(n_n);
cplex.addRows(zeros(n_n,1),A,zeros(n_n,1));

% if strcmp(cost_type,'Linf')
A = sparse(num.slack_np,nvar);
A(:,ndx.slack_np)=speye(num.slack_np);
A(:,ndx.slack_nn)=speye(num.slack_np);
A(:,ndx.slack_max) = -ones(num.slack_nn,1);
cplex.addRows(-inf*ones(num.slack_nn,1),A,zeros(num.slack_nn,1));
% end




cplex.solve();
status = cplex.Solution.status;
if ~isempty(cplex.Solution.x)
    sol0 = cplex.Solution.x;
    
    for v = var_name
        var = v{1};
        sol.(var) = sol0(ndx.(var));
        obj=cplex.Solution.objval;
    end
    %     sol.slack_p = (phi_p*sol.w(:,1).*sigma_p+phi_p*sol.w(:,2).*(1-sigma_p))';
    %     sol.slack_nnnnnnn = (phi_n*sol.w(:,1).*sigma_n+phi_n*sol.w(:,2).*(1-sigma_n))';
    
    %     sol.class_sense=class_sense;
else
    sol=[];
    obj=inf;
end
end














