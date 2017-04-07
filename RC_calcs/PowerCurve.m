clc
clear all
close all
Tav = 21000*3;
W = 84366.76;
rho = 0.002378*32.174;
S = 797.5934;
CL = [.5:0.01:1.4];
CD0 = 0.0474;
A = 3;
e = 0.8;
% Tav = 0.75*Tav;
RC = sqrt(2.*W./(rho.*S)).*( (Tav/W).*(CL).^(-0.5) - ((CD0+ (CL.^2./(pi*e*A)))./CL.^(1.5)));
max(RC)
% RC = RC.*60;
V = sqrt(2*W./(rho.*CL*S));
plot(V,RC)
ylabel('Rate of climb (fpm)')
xlabel('Velocity (ft/s)')
title('Takeoff/Landing')
%% ========================= SERVICE CEILING =========================== %%
CD0 = 0.0214;
h = 40:1:75;
[a,b,sigma,d] = AltTable(h,'h');
CL = [0.01:0.01:2];
V = sqrt(2*W./(rho.*CL*S));
TSL = 21000*3;
K = 0.21;
for i = 1:length(h)
rho1 = rho*sigma(i);
T = TSL.*sigma(i).*(1+K.*(V./(d(i).*1116.2)));
CD = CD0+1./(pi*e*A).*CL.^2;
Treq = 0.5.*rho1.*(V.^2).*S.*CD;
% Treq = (CD0+ (CL.^2./(pi.*e.*A))).*1481.3.*a(i).*(V./(d(i).*1116.2)).^2*S;
% Treq = (CD0+ (CL.^2./(pi.*e.*A))).*1481.3.*(V./(d(i).*1116.2)).^2*S;
% Wod = CL.*1481.3.*(V./(d(i).*1116.2)).^2*S;
% RC = sqrt(2.*W./(rho1.*S)).*( (Tav/W).*(CL).^-0.5 - ((CD0+ (CL.^2./(pi*e*A)))./CL.^(1.5)));
RC = ((T-Treq).*V)./W;

% RC = RC.*60;
RCmax(i) = max(RC);
end
figure(2)
hold on
plot(RCmax,h)




xlabel('Rate of climb (fpm)','Fontname','Times New Roman')
ylabel('Altitude (kft)','Fontname','Times New Roman')
title('Service Ceiling','Fontname','Times New Roman')
% axis([0 180 0 35])
plot([100 100], ylim)
hold on
y = spline(RCmax,h,100);
plot(xlim,[y y])
str = 'Service Ceiling (kft) = %0.2f';
text(105,55,sprintf(str,y))
    
% rc_sc = 100*ones(1,length(RC_vec));
% [~,SC] = polyxpoly(RC_vec,h,rc_sc,h);  % Service Ceiling (kft)

% fprintf('Service Ceiling = %0.1f kft\n',SC)