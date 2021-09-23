%  coherent demodulation of differentially encoded binary phase shift keying (DBPSK)
clear
close all
clc;
N = 10^6; 
rand('state',100); 
randn('state',200);
ip = rand(1,N)>0.5;

ipD = mod(filter(1,[1 -1],ip),2);
s = 2*ipD-1;
n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)];

Eb_N0_dB = [-3:10];
for ii = 1:length(Eb_N0_dB)
y = s + 10^(-Eb_N0_dB(ii)/20)*n;

ipDHat_coh = real(y) > 0; 
ipHat_coh = mod(filter([1 -1],1,ipDHat_coh),2);
nErr_dbpsk_coh(ii) = size(find([ip - ipHat_coh]),2);
end
simBer_dbpsk_coh = nErr_dbpsk_coh/N;

theoryBer_dbpsk_coh = erfc(sqrt(10.^(Eb_N0_dB/10))).*(1 - 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))));

close all
figure
semilogy(Eb_N0_dB,theoryBer_dbpsk_coh,'b.-');
hold on
semilogy(Eb_N0_dB,simBer_dbpsk_coh,'mx-');
axis([-2 10 10^-6 0.5])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for coherent demodulation of DBPSK')