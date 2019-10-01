function [x_bar,dx_bar] = bezier_regression(x,t)
%% use bezier polynomial to fit x locally, and generate filtered x and dxdt.

global bez_reg Bez_matr Bez_matr_der Bez_matr_dder N_bezier dB ddB
n = 9;
N = length(x);
Ts = (t(end)-t(1))/(length(t)-1);
if size(x,1)==1
    x = x';
end
for i=1:length(x)
    idx_low  = max(1,i-n);
    idx_high = min(N,i+n);
    idx = idx_low:idx_high;               % regression window
    idx0 = find(~isnan(x(idx)));          % indices that are not NAN
    idx = idx(idx0);
    M = idx_high-idx_low+1;
    T=(M-1)*Ts;
    if length(idx)< 1.5*N_bezier
        x_bar(i) = nan;
        dx_bar(i) = nan;
    else
        if length(idx)==M
            reg=bez_reg{M};
        else
            A = Bez_matr{M}(idx0,:);
            reg = (A'*A)\(A');
        end 
        alpha=reg*x(idx);                % get the bezier coefficient with pre-computed regressor
        
        dalpha=dB*alpha;                 % get the bezier coefficient for the derivative
        x_temp = Bez_matr{M}*alpha;
        dx_temp = Bez_matr_der{M}*dalpha/T;
        x_bar(i)=x_temp(1+i-idx_low);
        dx_bar(i)=dx_temp(1+i-idx_low);
    end
    if abs(x_bar(i))>10
        disp('')
    end
end
if size(x_bar,1) ==1
    x_bar = x_bar';
end
if size(dx_bar,1) ==1
    dx_bar = dx_bar';
end