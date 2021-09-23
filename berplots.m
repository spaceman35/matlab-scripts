clear all
close all
EbNodB = 0:0.05:25;
EbNo   = 10.^(0.1*EbNodB);
dBticks = 0:1:25;

%% QAM plots
M=4;    L=sqrt(M); Pb4QAM   = (2*(1-1/L)/log2(L))*qfunc(sqrt((3*log2(L)/(L^2-1))*2*EbNo));
M=16;   L=sqrt(M); Pb16QAM  = (2*(1-1/L)/log2(L))*qfunc(sqrt((3*log2(L)/(L^2-1))*2*EbNo));
M=64;   L=sqrt(M); Pb64QAM  = (2*(1-1/L)/log2(L))*qfunc(sqrt((3*log2(L)/(L^2-1))*2*EbNo));
M=256;  L=sqrt(M); Pb256QAM = (2*(1-1/L)/log2(L))*qfunc(sqrt((3*log2(L)/(L^2-1))*2*EbNo));

figure
semilogy(EbNodB,Pb4QAM,'b','LineWidth',2)
hold on
semilogy(EbNodB,Pb16QAM,'r-.','LineWidth',2)
semilogy(EbNodB,Pb64QAM,'k--','LineWidth',2)
semilogy(EbNodB,Pb256QAM,'g:','LineWidth',2)
grid on
grid minor
ylim([1e-7 1])
xticks(dBticks)
xlabel('E_b/N_0 (dB)')
ylabel('Probability of bit error, P_b')
legend('4QAM','16QAM','64QAM','256QAM')
title('Rectangular QAM BER Curves')

%% PSK plots
PbBPSKQPSK = qfunc(sqrt(2*EbNo));
PbDEBPSK = 2*qfunc(sqrt(2*EbNo));
PbDPSK = 0.5*exp(-EbNo);

m=3; M=2^m; Pb8PSK  = 2*qfunc(sqrt(2*m*EbNo)*sin(pi/M))/m;
m=4; M=2^m; Pb16PSK = 2*qfunc(sqrt(2*m*EbNo)*sin(pi/M))/m;
m=5; M=2^m; Pb32PSK = 2*qfunc(sqrt(2*m*EbNo)*sin(pi/M))/m;

figure
semilogy(EbNodB,PbBPSKQPSK,'b','LineWidth',2)
hold on
semilogy(EbNodB,PbDEBPSK,'r-.','LineWidth',2)
semilogy(EbNodB,PbDPSK,'k--','LineWidth',2)
semilogy(EbNodB,Pb8PSK,'b:','LineWidth',2)
semilogy(EbNodB,Pb16PSK,'r-','LineWidth',2)
semilogy(EbNodB,Pb32PSK,'g-.','LineWidth',2)

grid on
grid minor
ylim([1e-7 1])
xticks(dBticks)
xlabel('E_b/N_0 (dB)')
ylabel('Probability of bit error, P_b')
legend('BPSK/QPSK','DE-BPSK','DPSK','8PSK','16PSK','32PSK')
title('MPSK BER Curves')

%% mixed plots
figure
semilogy(EbNodB,PbBPSKQPSK,'b','LineWidth',2)
hold on
semilogy(EbNodB,PbDPSK,'k--','LineWidth',2)
semilogy(EbNodB,Pb8PSK,'b:','LineWidth',2)
semilogy(EbNodB,Pb16PSK,'r-.','LineWidth',2)
semilogy(EbNodB,Pb16QAM,'k-.','LineWidth',2)

grid on
grid minor
ylim([1e-7 1])
xticks(dBticks)
xlabel('E_b/N_0 (dB)')
ylabel('Probability of bit error, P_b')
legend('BPSK/QPSK','DPSK','8PSK','16PSK','16QAM')
title('Comparison of Selected BER Curves')
