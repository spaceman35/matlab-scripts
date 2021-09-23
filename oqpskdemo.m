% QPSK, OQPSK, MSK demo
% Bob Mills
% Feb 2020
%
clc; clearvars; close all

%% simulation parameters
% plot control variables
p.discrete = 0;
p.spectrum = 0;

% simulation parameters
Nsyms = 20;                     % number of symbols in simulation
Nbits = Nsyms*2;                % number of bits in simulation
fsamp = 100;                    % sampling freq
Nfft = 4096;                    % FFT size for freq plots

% signal parameters
f0 = 3;                        % carrier freq
Rsym = 1;                       % symbol rate 
Tsym = 1/Rsym;                  % symbol dwell time

% derived parameters
Ns = Tsym*fsamp;                % num samples / symbol
Nsamp = Nsyms*Ns;               % total num samples in simulation
A = 1;                          % carrier amplitude to get desired Eb

% time vectors
t = (0:Nsamp-1)'/fsamp;         % time base vector (column)
ts = (0:Ns-1)'/fsamp;           % time vector for a single symbol (column)

%% generate the modulated signals 
% first generate the I/Q bit streams 
bits = 2*randi(2,Nbits,1)-3;    % total # of bits, values set to +/- 1
bits = reshape(bits,2,Nbits/2);
Ibits = bits(1,:);              % I channel bits
Qbits = bits(2,:);              % Q channel bits
bits = bits(:);

% build the pulse for the I&Q bit streams
pulse = ones(Ns,1);                         % NRZ pulse
% pulse  = [ ones(Ns/2,1); zeros(Ns/2,1) ];   % RZ pulse
dI = pulse * Ibits; dI = dI(:);
dQ = pulse * Qbits; dQ = dQ(:);
theta_t = atan2(dQ,dI) * 180/pi;

% modulate onto carrier -- conventional QPSK
QPSK = A * ( dI .* cos(2*pi*f0*t) + dQ .* sin(2*pi*f0*t) );

% OQPSK has a 1/2 symbol delay in the Q channel
% or 1/2 symbol advance in the I channel
% for simplicity, just circular wrap from L to R end of dI vector
dIo = circshift(dI,-Ns/2);
thetao_t = atan2(dQ,dIo) * 180/pi;
OQPSK = A * ( dIo .* cos(2*pi*f0*t) + dQ .* sin(2*pi*f0*t) );

% MSK / CPFSK uses half sine pulses
% just multiply the dI and dQo pulses from OQPSK with half sine pulse shape
COS = cos(pi*t/(Tsym));
SIN = sin(pi*t/(Tsym));
dIm = COS .* dIo;
dQm = SIN .* dQ;
thetam_t = atan2(dQm,dIm) * 180/pi;
MSK = A * ( dIm .* cos(2*pi*f0*t) + dQm .* sin(2*pi*f0*t) );


%% plot time signals 
% Conventional QPSK %%%%%%%%%%%%%
figure
subplot(411)
if p.discrete
    plot(t,dI,'.')
else
    plot(t,dI)
end
ylim([-1.5 1.5])
ylabel('In-phase d_I(t)')
title('Conventional QPSK')

subplot(412)
if p.discrete
    plot(t,dQ,'.')
else
    plot(t,dQ);
end
ylim([-1.5 1.5])
xlabel('time, s')
ylabel('Quadrature d_Q(t)')

subplot(413)
if p.discrete
    plot(t,theta_t,'.')
else
    plot(t,theta_t);
end
ylim([-180 180])
ylabel('Phase angle \theta(t)^o')

subplot(414)
if p.discrete
    plot(t,QPSK,'.')
else
    plot(t,QPSK)
end
title('Modulated signal')
xlabel('time, s')
ylabel('Signal s(t)')
ylim([ -1.8*A 1.8*A ])

% Offset QPSK %%%%%%%%%%%%%
figure
subplot(411)
if p.discrete
    plot(t,dIo,'.')
else
    plot(t,dIo)
end
ylim([-1.5 1.5])
ylabel('In-phase d_I(t)')
title('Offset QPSK')

subplot(412)
if p.discrete
    plot(t,dQ,'.')
else
    plot(t,dQ);
end
ylim([-1.5 1.5])
xlabel('time, s')
ylabel('Quadrature d_Q(t)')

subplot(413)
if p.discrete
    plot(t,thetao_t,'.')
else
    plot(t,thetao_t);
end
ylim([-180 180])
ylabel('Phase angle \theta_O(t)^o')

subplot(414)
if p.discrete
    plot(t,OQPSK,'.')
else
    plot(t,OQPSK)
end
title('Modulated signal')
xlabel('time, s')
ylabel('Signal s(t)')
ylim([ -1.8*A 1.8*A ])

% MSK %%%%%%%%%%%%%
figure
subplot(411)
if p.discrete
    plot(t,dIm,'.')
else
    plot(t,dIm)
end
ylim([-1.5 1.5])
ylabel('In-phase d_I(t)')
title('MSK/CPFSK')

subplot(412)
if p.discrete
    plot(t,dQm,'.')
else
    plot(t,dQm);
end
ylim([-1.5 1.5])
xlabel('time, s')
ylabel('Quadrature d_Q(t)')

subplot(413)
if p.discrete
    plot(t,thetam_t,'.')
else
    plot(t,thetam_t);
end
ylim([-180 180])
ylabel('Phase angle \theta_m(t)^o')

subplot(414)
if p.discrete
    plot(t,MSK,'.')
else
    plot(t,MSK)
end
title('Modulated signal')
xlabel('time, s')
ylabel('Signal s(t)')
ylim([ -1.5*A 1.5*A ])

figure %% more MSK
MSK_I = dIm.*cos(2*pi*f0*t);
MSK_Q = dQm.*sin(2*pi*f0*t);
subplot(311)
plot(t,MSK_I,t,COS,'k:',t,-COS,'k:')
ylim([-1.5 1.5])
ylabel('dI cos(\pi t/T_S)cos(2\pi f_0 t)')
title('Modulated MSK signal')

subplot(312)
plot(t,MSK_Q,t,SIN,'k:',t,-SIN,'k:')
ylim([-1.5 1.5])
xlabel('time, s')
ylabel('dQ sin(\pi t/T_S)sin(2\pi f_0 t)')

subplot(313)
plot(t,MSK)
ylabel('MSK signal')
xlabel('time, s')
ylim([ -1.5*A 1.5*A ])

%% plot freq domain results
if p.spectrum
    QPSKf = db( fftshift( abs( fft(QPSK,Nfft) ) ), 'power');
    dIf = db( fftshift( abs( fft(dI,Nfft) ) ), 'power');
    dQf = db( fftshift( abs( fft(dQ,Nfft) ) ), 'power');
    fk = fsamp*((0:Nfft-1)/Nfft-.5);         % freq vector for FFT plots

    figure
    subplot(311)
    plot(fk,QPSKf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Mod Sig Amp, dB')
    title('Power Spectrum, QPSK')

    subplot(312)
    plot(fk,dIf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('In-Phase Amp, dB')

    subplot(313)
    plot(fk,dQf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Quadrature Amp, dB')

    OQPSKf = db( fftshift( abs( fft(OQPSK,Nfft) ) ), 'power');
    dIof = db( fftshift( abs( fft(dI,Nfft) ) ), 'power');
    figure
    subplot(311)
    plot(fk,OQPSKf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Mod Sig Amp, dB')
    title('Power Spectrum, OQPSK')

    subplot(312)
    plot(fk,dIof)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('In-Phase Amp, dB')

    subplot(313)
    plot(fk,dQf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Quadrature Amp, dB')

    MSKf = db( fftshift( abs( fft(MSK,Nfft) ) ), 'power');
    dImf = db( fftshift( abs( fft(dIm,Nfft) ) ), 'power');
    dQmf = db( fftshift( abs( fft(dQm,Nfft) ) ), 'power');
    figure
    subplot(311)
    plot(fk,MSKf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Mod Sig Amp, dB')
    title('Power Spectrum, MSK/CPFSK')

    subplot(312)
    plot(fk,dImf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('In-Phase Amp, dB')

    subplot(313)
    plot(fk,dQmf)
    axis([-fsamp/2 fsamp/2 -20 30])
    xlabel('freq, Hz')
    ylabel('Quadrature Amp, dB')
end
