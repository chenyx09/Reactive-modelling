function [P,q]=inner_approx_circle(R,n,center)
if size(center,2)>1
    center=center';
end

theta = (2*pi/n:2*pi/n:2*pi)';
R1=R*cos(2*pi/n);

P=[cos(theta) sin(theta)];
q=R1*ones(n,1)+P*center;
