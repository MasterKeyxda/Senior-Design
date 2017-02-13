clc
clear all
W_TOg = 60000;
Wfuelratio = 0.5293*1.05;
W_Pl = 1720;
W_Crew = 340;

W_fg = Wfuelratio*W_TOg;

A = 0.4221;
B = 0.9876;
res = 1e6;
i = 0
while res > 1 && i < 1e7
    W_OEt = W_TOg - W_fg - W_Pl;
    W_E = 10^((log10(W_TOg) - A)/B);
    W_Etent = W_OEt - W_Crew;
    
    res = abs(W_E - W_Etent);
    W_TT = W_TOg;
    W_TOg = W_Pl + W_fg + W_Etent; 
    W_fg = Wfuelratio*W_TOg;
    W_TOg = W_Pl + W_fg + W_OEt;
    i = i +1;
end