function [A,b]=polygon_circle(x,R,n)
theta=pi/n:2*pi/n:2*pi-pi/n;
theta = theta';
A=zeros(n,2);
b=zeros(n,1);
% for i=1:n
%     if theta(i)==0
%         H=[1 0];
%         h=R+x(1);
%     elseif theta(i)==pi
%         H=[-1 0];
%         h=R-x(1);
%     else
%         if sin(theta(i))>0
%             H=[cot(theta(i)) 1];
%             h=x(2)+R/sin(theta(i))+cot(theta(i))*x(1);
%         else
%             H=[-cot(theta(i)) -1];
%             h=-x(2)-R/sin(theta(i))-cot(theta(i))*x(1);
%         end
%     end
%     A(i,:)=H;
%     b(i)=h;
% end
A = [cos(theta) sin(theta)];
b = x(1)*cos(theta)+x(2)*sin(theta)+R*ones(n,1);
% A=A./b;
% b=b./b;