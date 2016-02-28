function [Jdomain,qnode,radius] = Jdomain(tip_elem,xTip)

global node element

numnode = size(node,1);
% -------------------------------------
% calculation of the area of the tip element
x = node(element(tip_elem,:),:);
% Area = sum of areas of each sub-triangle
x0 = x(1,1);
y0 = x(1,2);

x1 = x(2,1);
y1 = x(2,2);

x2 = x(3,1);
y2 = x(3,2);

x3 = x(4,1);
y3 = x(4,2);

A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;
A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
area = A1 + A2 ;

% J radius = fac * sqrt(area);
fac    = 1.0;
radius = fac * sqrt(area);
center = xTip;
%radius=4.53;

r=[];
% Distance from the center of tip element
for i = 1 : numnode
    sctr = node(i,:);
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
    r    = [r,rho];
end

test = r-radius;
test = test(element)';
test = max(test).*min(test);
Jdomain = find(test<=0);
test1 = r-radius;
test1 = test1(element(Jdomain,:))';
test1 = (test1<=0);
qnode = test1';

