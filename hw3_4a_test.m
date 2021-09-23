%HW3 4a
%% EENG 571 Matlab Simulation Project
% Dr Mills
% Winter 2021
%
% The purpose of this project is to reinforce understanding of how 
% digital signals are demodulated in the presence of additive noise.
% 
% The code below will implement a Monte Carlo simulation of a 
% quadrature phase shift keying (QPSK) communication system. You will 
% need to complete the code and then simulations using different 
% parameters in order to provide answers to several questions posed in
% the project assignment. Once you have the code working and calibrated, 
% Id suggest keeping a copy in case you make a change that breaks 
% everything.)
%
% See slides 31-35 in slide set 04b for more info on QPSK implementation.
%

clearvars; close all; clc; clear all;

% plot switches - turn off/on as desired 
Plot.DataSignal  = 0;               % plot I/Q data waveforms
Plot.ModSignal   = 0;               % plot modulated MPSK signal
Plot.Noise       = 0;               % plot noise signal
Plot.SigNoise    = 0;               % plot RF signal + noise
Plot.Correlator  = 0;               % plot correlator outputs
Plot.Histograms  = 0;               % plot correlator histograms
Plot.Histogram3D = 0;               % plot 3D histogram
Plot.Scatterplot = 0;               % plot scatterplot w/o errors
Plot.SymErrors   = 0;               % plot scatterplot w/ errors
Plot.CondPDFs    = 0;               % plot cond'l PDFs using truth data

% simulation/system variables
fs    = 100;                        % Sample rate [samples/sec]
                                    % fs should be at least 2(fs+2Rb)
Nfft  = 2^16;                       % # FFT pts--chg freq res in spec plots
Nsyms = 1000;                        % # symbols to use in simulation

% signal parameters -- OK to change 
Rsym    = 1;                        % Symbol rate [sym/sec]
f0      = 20;                       % Carrier frequency [Hz]
EbNo_dB = 0;                   % Eb/No Ratio [dB] 

for EbNo_dB=0:10
    
% receiver synch errors
RcvError.phi_deg = 0;               % rcvr phase error in degrees
RcvError.phi_rad = RcvError.phi_deg*pi/180;
RcvError.freq_eps = 0;              % f0*1e-6 = 1 PPM relative to f0

% Fixed / dependent parameters... Be careful adjusting anything here
NoiseSigma = 1;                     % Variance of AWGN samples
M     = 4;                          % QPSK -- code won't work if changed
m     = log2(M);                    % Number of bits/symbol
Nbits = Nsyms * m;                  % Number of bits 
Tsym  = 1/Rsym;                     % Period of 1 symbol [sec/sym]
No    = 2*NoiseSigma^2/fs;          % Noise power spectral density
Ns    = fs*Tsym;                    % # samples per symbol [samples/sym]
Rb    = m*Rsym;                     % Bit rate [bits/sec]

% Time vector
t = (0:1/fs:(Nsyms*Tsym-(1/fs)));   % Time vector from 0 to max [sec]
t = t(:);                           % Force column vector

% Set amplitude of carrier signal to get desired Eb/No
% You'll need to complete all 3 lines -- i.e., convert Eb/No from dB, then
% get bit energy Eb, and then determine A that will produce that Eb based
% on the bit rate.
%
% Note it's possible to get the simulation to "run" ok even if you don't 
% have correct entries below. You'll need to do some simulations 
% to verify the simulated results agree with theory. 
EbNo   = 10.^(EbNo_dB/10);          % Eb/No (linear ratio)
C      = EbNo*Rb*No;
Eb     = C/Rb;                     % Energy per bit [J/bit]
A      = sqrt(2*C);                % Amplitude [Volts]

%% generate vector of input bits to send
Input.Bits = randi(2,Nbits,1)-1;

% then reshape into m rows where m=log2(M) bits/symbol
Input.Bits = reshape(Input.Bits,m,Nsyms);
Input.B1   = Input.Bits(1,:);       % row 1
Input.B2   = Input.Bits(2,:);       % row 2
Input.Bits = Input.Bits(:);         % orig bits back to single vector

% the following lines do a gray code bit-to-symbol mapping
% first bit is from B1, 2nd bit is from B2
% see slide 62 in slide set #4
% higher values of M would need more rows and angle points
% 11 = upper right -- phase angle num 0
% 01 = upper left  -- phase angle num 1 
% 00 = lower left  -- phase angle num 2 
% 10 = lower right -- phase angle num 3
% actual angle = (2*pi)*angle/M + pi/M
Input.Angle_num(Input.B1==1 & Input.B2==1) = 0;
Input.Angle_num(Input.B1==0 & Input.B2==1) = 1;
Input.Angle_num(Input.B1==0 & Input.B2==0) = 2;
Input.Angle_num(Input.B1==1 & Input.B2==0) = 3;
Input.Angle_rad = (Input.Angle_num*2*pi/M) + pi/M;

%% for QPSK, we can use the above angle mapping, or we can just apply 
% the even / odd bits directly to I & Q (normalize to +1/-1)
pulse = ones(Ns,1);                 % NRZ pulse shape
Data.I = pulse*(2*Input.B1-1);      % pulses normalized to +1/-1
Data.Q = pulse*(2*Input.B2-1);
Data.Phase = pulse*Input.Angle_rad;

% force all to column vectors so they align with time vector
Data.I = Data.I(:); 
Data.Q = Data.Q(:);
Data.Phase = Data.Phase(:);

if Plot.DataSignal
    Plot_IQData(Data,fs);
end

%% modulate onto carrier -- conventional QPSK
% In following line, construct the QPSK signal
% Note that Data.I and Data.Q represent vI(kTs) and vQ(kTs), respectively
% You can implement either using I&Q modulation (2 mixing operations) or
% by adjusting the phase angle of the carrier using Data.Phase. 
% In fact, you might want to implement both approaches to see that they 
% are equivalent. 
Signal.Data =  Data.I.*(A/sqrt(2)).*cos(2*pi*f0*t)-Data.Q.*(A/sqrt(2)).*sin(2*pi*f0*t);

Signal.PSD = ComputePSD(Signal.Data,Nfft);

if Plot.ModSignal
    Plot_Signal(Signal,fs,'Modulated Signal');
end

% AWGN - zero mean, variance = NoiseSigma^2
Noise.Data = NoiseSigma*randn( size(t) );
Noise.PSD = ComputePSD(Noise.Data,Nfft);
if Plot.Noise
    Plot_Signal(Noise,fs,'AWG Noise');
end

% received signal + noise
SigNoise.Data = Signal.Data + Noise.Data;
SigNoise.PSD = ComputePSD(SigNoise.Data,Nfft);
if Plot.SigNoise
    Plot_Signal(SigNoise,fs,sprintf('Signal+Noise,Eb/No=%0.1f dB',EbNo_dB));
end

%% correlation receiver
% local oscillator signals used for correlation in I and Q 

% shift received signal by timing error ...noise is shifted too, but 
% it's random so it doesn't matter
CorrOutput.I = SigNoise.Data .* cos(2*pi*f0*(1+RcvError.freq_eps)*t ...
        + RcvError.phi_rad);
CorrOutput.Q = - SigNoise.Data .* sin(2*pi*f0*(1+RcvError.freq_eps)*t ...
        + RcvError.phi_rad);

% reshape correlator outputs into matrices 
% each column is a symbol, Nsyms = # of columns
% then integrate columns 

%Sum.I = cumsum(reshape(CorrOutput.I,Ns,Nsyms));
%Sum.Q = cumsum(reshape(CorrOutput.Q,Ns,Nsyms));
Sum.I = cumsum(reshape(CorrOutput.I,Ns,Nsyms));
Sum.Q = cumsum(reshape(CorrOutput.Q,Ns,Nsyms));
if Plot.Correlator
    Plot_Correlator(Sum,fs,EbNo_dB)
end

% last row is the output of integrator at end of symbol
Test.I = Sum.I(end,:);
Test.Q = Sum.Q(end,:);

if Plot.Histograms
    Plot_Histograms(Test,EbNo_dB);
end

if Plot.CondPDFs
    Plot_CondPDFs(Input,Test,EbNo_dB);
end

if Plot.Scatterplot
    Plot_Scatterplot(Test,EbNo_dB);
end

if Plot.Histogram3D
    Plot_Histogram3D(Test,EbNo_dB);
end


%% recover symbols (phase angles) based on correlator outputs
% use same bit-to-symbol mapping as before, i.e., positive correlator 
% output --> "1" otherwise "0", thus "11" is upper right corner, and 
% "00" is lower left corner
% for higher level MPSK, need to determine angle then recover the bits
Output.Angle_num = zeros(size(Input.Angle_num));
Output.Angle_num(Test.I>=0 & Test.Q>=0) = 0;
Output.Angle_num(Test.I<0  & Test.Q>=0) = 1;
Output.Angle_num(Test.I<0  & Test.Q<0)  = 2;
Output.Angle_num(Test.I>=0 & Test.Q<0)  = 3;

% find any symbol errors and count them
Error.Symbols = Output.Angle_num ~= Input.Angle_num;
Error.Psym = sum(Error.Symbols)/Nsyms;

if Plot.SymErrors
    Plot_SymErrors(Test,Error,EbNo_dB);
end

% get the recovered bits and reconstruct the bit stream
Output.B1 = (Test.I >= 0);
Output.B2 = (Test.Q >= 0);
Output.Bits = [Output.B1; Output.B2];
Output.Bits = Output.Bits(:);

Error.Bits = Output.Bits ~= Input.Bits;
Error.Pbit = sum(Error.Bits)/Nbits;
Error.Theory = qfunc(sqrt(2*EbNo));

%% Summary output
fprintf('MPSK with M=%d, Eb/No=%6.1f dB\n',M,EbNo_dB);
fprintf('Synch errors --  freq: %0.2e           phase: %0.1f deg\n', ...
        RcvError.freq_eps, RcvError.phi_deg);
fprintf('%1.0e syms    %1.0e sym errors     Est Psym: %5.1e\n', ...
        Nsyms, sum(Error.Symbols), Error.Psym);
fprintf('%1.0e bits    %1.0e bit errors     Est Pbit: %5.1e\n', ...
        Nbits, sum(Error.Bits), Error.Pbit );
fprintf('                                   Thr Pbit: %5.1e\n',Error.Theory);


semilogy(EbNo_dB,Error.Theory,'b.-');
hold on
semilogy(EbNo_dB,Error.Psym,'mx-');
axis([0 10 10^-7 0.5])
grid on
legend('Theory', 'Simulated');
xlabel('E_b/N_0 (dB)')
ylabel('Probability of bit error, P_b')
title('Pb versus E_b/N_0')
hold on
end


function PSD = ComputePSD(input,Nfft)
    PSD = 20*log10(abs(fftshift(fft(input,Nfft))));
end



