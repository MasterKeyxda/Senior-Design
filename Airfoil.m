t = 0.07;
c = 1;
dx = 0.01;
%% Segment 1 - Biconvex
x1 = 0:dx:0.1;
y1 = 2*(t/c).*x1.*(1-x1);
dydx1 = 2*(t/c) - 4*(t/c)*x1(end);
%% Segment 2 - Ellipse
t2 = t - 0.5*t;
x2start = (-dydx1 + 2*(t2/c))/(4*(t2/c));
x2 = (0.1+x2start):dx:(0.9+x2start);
% b = sqrt((dydx1*(1/x2(1)))^2*(1-x2(1)^2));
% y2 = sqrt((1- (x2.^2/(c^2))*b^2));
x2 = x2 - x2start;
y2 = 2*(t2/c).*x2.*(1-x2);
yup = y1(end) - y2(1);
y2 = y2 + yup;

%% Segment 3 - Biconvex
x3 = 0.9:dx:1;
y3 = 2*(t/c).*x3.*(1-x3);


%% Segment 4 - Bottom Biconvex

x4 = 1:-dx:0;
y4 = -2*(t/c).*x4.*(1-x4);

plot([x1,x2,x3,x4],[y1,y2,y3,y4])
ylim([-.2 .2])
hold on
camber = [y1(1:end-1),y2(1:end-1),y3]+y4;
plot(x4,camber)
%% Two Biconvex

t = 0.07;

x5 = 0:dx:1;
y5 = 2*((t-0.5*t)/c).*x5.*(1-x5);
figure

plot(x4,y4)
hold on
plot(x5,y5)
ylim([-.2 .2])
plot(x4,y4+y5)