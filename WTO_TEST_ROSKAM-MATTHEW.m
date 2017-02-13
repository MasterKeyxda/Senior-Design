% Inital Guess
W_PL = 1960;
%W_F = 40000;
W_Crew = 400;
W_TOguess = 74600;
A = 0.4221; B = 0.9876;
W_OE_tent = W_TOguess - W_F - W_PL;
W_E_tent = W_OE_tent - W_Crew;
W_E = 10.^((log10(W_TOguess)-A)./B);
fprintf('W_E: %0.2f \n',W_E) 
fprintf('W_TO: %0.2f \n',W_TOguess) 

%%
W_FTO = 0.49;
A = 1.02;
C = -0.06;
y = linspace(80000,190000,100);
x1 = A*y.^C;
x2 = 0.51 - 2360./y;
plot(y,x1,y,x2)

