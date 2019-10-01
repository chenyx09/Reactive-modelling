function []=extract_reactive_bound_2region_soft_bdry(w_val,gamma,threshold_x)
global app
params
switch app
    case 1
M = size(p,1);
tol = 1e-10;
for i=1:2
    for j=1:size(w_val{i},1)
        for k=1:size(w_val{i},2)
            if abs(w_val{i}(j,k))<tol
                w_val{i}(j,k)=0;
            end
        end
    end
end

n_class = 2;
n_plane = size(w_val{1},1);


filename='LC_reactive_bound.m';
f=fopen(filename,'w');
fprintf(f,'function au = LC_reactive_bound(x)\n');

fprintf(f,['x_norm=',mat2str(x_norm),';\n']);
fprintf(f,['x_min=',mat2str(x_min),';\n']);
fprintf(f,'y=(x-x_min)./x_norm(1:4);\n');
fprintf(f,'m1=1/(1+exp(%d*(y(1)-%d)));\n',[gamma,threshold_x]);
fprintf(f,'m2=1-m1;\n');
fprintf(f,'ub = zeros(%d,1);\n',n_plane);
for i=1:n_plane
    fprintf(f,'ub(%d)=(',i);
    for n=1:n_class
        fprintf(f,'m%d*(',n);
        for j=1:M-1
            if w_val{n}(i,j)
                fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d+...\n',[p(j,1:end-1) w_val{n}(i,j)]);
            end
        end
        fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d)...\n',[p(M,1:end-1) w_val{n}(i,M)]);
        if n<n_class
            fprintf(f,'+');
        end
    end
    fprintf(f,')*%d/(-%d*m1-%d*m2)+%d;\n',a_norm,w_val{1}(i,M+1),w_val{2}(i,M+1),a_min);
end
fprintf(f,'au = min(ub);\n');


fprintf(f,'end');

fclose(f);

    case 2
A1_idx  = find(p(:,1));
A2_idx = find(p(:,2));
b_idx = find(~p(:,1).*~p(:,2));

% digits(4);
% sol.w=vpa(sol.w);
tol = 1e-10;
for i=1:2
    for j=1:size(w_val{i},1)
        for k=1:size(w_val{i},2)
            if abs(w_val{i}(j,k))<tol
                w_val{i}(j,k)=0;
            end
        end
    end
end

n_class = 2;
n_plane = size(w_val{1},1);

% K=length(w_val);


% D=5;

filename='robot_reactive_bound.m';
f=fopen(filename,'w');
fprintf(f,'function [A,b] = robot_reactive_bound(x)\n');
%
% fprintf(f,'n_hidden=[')
% dlmwrite(filename,n_hidden,'-append');
% fseek(f, 0, 'eof');
% fprintf(f,'];\n');
% fprintf(f,['n_input=',num2str(n_input),';\n']);
% fprintf(f,['n_output=',num2str(n_output),';\n']);
% fprintf(f,['output=zeros(',num2str(n_output),',1);\n'])
fprintf(f,['A1=zeros(',num2str(n_plane),',2);\n']);
fprintf(f,['b1=zeros(',num2str(n_plane),',1);\n']);

fprintf(f,['A2=zeros(',num2str(n_plane),',2);\n']);
fprintf(f,['b2=zeros(',num2str(n_plane),',1);\n']);

fprintf(f,['A=zeros(',num2str(n_plane),',2);\n']);
fprintf(f,['b=zeros(',num2str(n_plane),',1);\n']);

fprintf(f,['x_norm=',mat2str(x_norm),';\n']);
fprintf(f,'y=x./x_norm(1:4);\n');
fprintf(f,'m1=1/(1+exp(%d*(y(1)-%d)));\n',[gamma,threshold_x]);
fprintf(f,'m2=1-m1;\n');


for n=1:n_class
    for i=1:n_plane
        fprintf(f,'A%d(%d,1)=(',[n,i]);
        for j=1:length(A1_idx)-1
            if w_val{n}(i,A1_idx(j))
                fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d+...\n',[p(A1_idx(j),3:end-1) -w_val{n}(i,A1_idx(j))]);
            end
        end
        fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d)/%d;\n',[p(A1_idx(end),3:end-1) -w_val{n}(i,A1_idx(end)) x_norm(5)]);
    end
    
    for i=1:n_plane
        fprintf(f,'A%d(%d,2)=(',[n,i]);
        for j=1:length(A2_idx)-1
            if w_val{n}(i,A2_idx(j))
                fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d+...\n',[p(A2_idx(j),3:end-1) -w_val{n}(i,A2_idx(j))]);
            end
        end
        fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d)/%d;\n',[p(A2_idx(end),3:end-1) -w_val{n}(i,A2_idx(end)) x_norm(6)]);
    end
    
    for i=1:n_plane
        fprintf(f,'b%d(%d)=',[n,i]);
        for j=1:length(b_idx)-1
            if w_val{n}(i,b_idx(j))
                fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d+...\n',[p(b_idx(j),3:end-1) w_val{n}(i,b_idx(j))]);
            end
        end
        fprintf(f,'y(1)^%d*y(2)^%d*y(3)^%d*y(4)^%d*%d;\n',[p(b_idx(end),3:end-1) w_val{n}(i,b_idx(end))]);
    end
end
fprintf(f,'A=A1*m1+A2*m2;\n');
fprintf(f,'b=b1*m1+b2*m2;\n');

fprintf(f,'end');

fclose(f);

end
end