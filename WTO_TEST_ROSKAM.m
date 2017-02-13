% % Inital Guess
% W_PL = 1720;
% W_F_WTO = 0.5293*1.05;
% W_Crew = 340;
% W_TOguess = 20000;
% W_F = W_F_WTO * W_TOguess;
% A = 0.4221; B = 0.9876;
% res = 1e7;
% i = 1;
% while res > 1 && i < 1e5
%     W_OE_tent = W_TOguess - W_F - W_PL;
%     W_E_tent = W_OE_tent - W_Crew;
%     W_E = 10.^((log10(W_TOguess)-A)./B);
%     res = abs(W_E_tent - W_E);
%     W_TOOO = W_TOguess;
%     W_F = W_F_WTO * W_TOguess;
%     W_TOguess = W_PL + W_F + W_OE_tent;
%     i = i + 1;
% end
% fprintf('W_E: %0.2f \n',W_E)
% fprintf('W_TO: %0.2f \n',W_TOOO)
% 
% %%
% W_FTO = 0.49;
% A = 1.59;
% C = -0.1;
% y = linspace(80000,170000,100);
% x1 = A*y.^C;
% x2 = 0.51 - 2360./y;
% plot(y,x1,y,x2)
% 
% % WE/WTO = 0.50
% % WF/WTO = 0.49
%%
W_PL = 2160;
W_F_WTO = 0.5293*1.05;
W_Crew = 400;
W_TOguess = 100000;
W_F = W_F_WTO * W_TOguess;
A = 0.4221; B = 0.9876;

    W_OE_tent = W_TOguess - W_F - W_PL;
    W_E_tent = W_OE_tent - W_Crew;
    W_E = 10.^((log10(W_TOguess)-A)./B);
    res = abs(W_E_tent - W_E);
    W_TOguess = W_PL + W_F + W_OE_tent;
    W_TOOO = W_TOguess;
    W_F = W_F_WTO * W_TOguess;
   
fprintf('W_E: %0.2f \n',W_E)
fprintf('W_TO: %0.2f \n',W_TOOO)