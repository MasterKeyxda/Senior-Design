% TF33 P-7
TF33.mdot = 498; % lb/s
TF33.SFC = 0.56;
TF33.Tt4 = 2210; % R
TF33.Tt0 = 518.69; % R
TF33.Force = 21000; %lbf
TF33.BPR = 1.21; % bypass ratio
TF33.FPR = 1.9; % Fan Pressure Ratio
TF33.OPR = 16.0;

% F100-PW-229
F100_PW_229 = TF33;
F100_PW_229.Force = 17800;
F100_PW_229.SFC = 0.74;
F100_PW_229.mdot = 248;
F100_PW_229.Tt4 = 3160; % R
F100_PW_229.BPR = 0.4;
F100_PW_229.FPR = 3.8;
F100_PW_229.OPR = 23; % total compressor pressure ratio

% F101-GE-102
F101_GE_102 = TF33;
F101_GE_102.Force = 17390;
F101_GE_102.SFC = 0.562;
F101_GE_102.mdot = 356;
F101_GE_102.Tt4 = 3010; % R
F101_GE_102.BPR = 1.91;
F101_GE_102.FPR = 2.31;
F101_GE_102.OPR = 26.8;