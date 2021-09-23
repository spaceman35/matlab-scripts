clc; clear all; clearvars; close all;
EbNoDEdB = 1:11;
EbNoDE = 10.^(EbNoDEdB/10);
Pb = erfc(sqrt(EbNoDE));
EbNoBPSK = (erfcinv(2*Pb)).^2;
EbNoBPSKdB = 10*log10(EbNoBPSK);
%Eb=[(erfcinv(2*1.157e-6)).^2]*10e-10
%EbNo=Eb/(10e-10)
%EbNodec=10*log(EbNo)


DeltadB = EbNoDEdB - EbNoBPSKdB;
fprintf('\n\n')
fprintf('Eb/No comparison, all values in decibels (dB)\n')
fprintf('w/ DE       w/o DE       Delta\n')
for k=1:length(EbNoDE)
   fprintf('%4.1f         %4.1f         %0.1f\n',...
       EbNoDEdB(k), EbNoBPSKdB(k), DeltadB(k))
end
figure
semilogy(DeltadB,Pb)
ylabel('P_b')
xlabel('\Delta (E_b/N_0) in dB')
title('Additional E_b/N_0 required for Differential Encoding in BPSK')
grid

