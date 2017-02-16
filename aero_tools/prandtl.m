%% Prandtl Lifting Line Theory (Sadraey pg 245)
function [CL_wing] = prandtl(N, S_w, AR, lambda, alpha_twist, i_w, a_2d, alpha_0)

b = sqrt(AR*S_w); % wing span
MAC = S_w/b; % Mean Aerodynamic Chord
Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord (m)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w+alpha_twist:-alpha_twist/(N-1):i_w;
% segment’s angle of attack
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % Mean Aerodynamics
Chord at each segment (m)
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1))/sin(theta(i)));
    end
end
A=B\transpose(LHS);
for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1 : N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end
CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];
plot(y_s,CL1,'-o')
grid on;
title('Lift distribution')
xlabel('Semi-span location (m)')
ylabel ('Lift coefficient')
CL_wing = pi * AR * A(1);
end