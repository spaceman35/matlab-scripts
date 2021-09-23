clear,close,clc all
function [pb,ps]=snr2p(snr_in_dB)

N=10000; % #symbols 
Es=1;
snr=10^(snr_in_dB/10);	 	  	
sgma=sqrt(Es/(4*snr)); % noise rms

% signal mapping
s00=[1 0];
s01=[0 1];
s11=[-1 0];
s10=[0 -1];

% generating signal source
for i=1:N
    temp=rand;			  	
    if (temp<0.25)            % s00	  
      dsource1(i)=0;
      dsource2(i)=0;		   
    elseif (0.25<=temp<0.5)   % s01
      dsource1(i)=0;
      dsource2(i)=1;
    elseif (0.5<=temp<0.75)   % s10
      dsource1(i)=1;	
      dsource2(i)=0;
    else			                % s11
      dsource1(i)=1;
      dsource2(i)=1;
    endif
endfor

numofsymbolerror=0;
numofbiterror=0;

% SER calculate using Box–Muller transform
for i=1:N
    n(1)=bmgauss(sgma);	  	  
    n(2)=bmgauss(sgma);
    if ((dsource1(i)==0) & (dsource2(i)==0))
        r=s00+n;    
    elseif ((dsource1(i)==0) & (dsource2(i)==1))
        r=s01+n;  
    elseif ((dsource1(i)==1) & (dsource2(i)==0))
        r=s10+n;
    else
        r=s11+n;  
    endif

    c00=dot(r,s00); % dot product
    c01=dot(r,s01);
    c10=dot(r,s10);
    c11=dot(r,s11);
    
    % decision
    c_max=max([c00 c01 c10 c11]);
    if (c00==c_max)
        decis1=0; decis2=0;
    elseif (c01==c_max)
        decis1=0; decis2=1;
    elseif (c10==c_max),
        decis1=1; decis2=0;
    else
        decis1=1; decis2=1;
    endif

    symbolerror=0;
    if (decis1~=dsource1(i))
        numofbiterror=numofbiterror+1;
        symbolerror=1;
    endif
    if (decis2~=dsource2(i))
        numofbiterror=numofbiterror+1;
        symbolerror=1;
    endif
    if (symbolerror==1)
        numofsymbolerror = numofsymbolerror+1;
    endif
endfor
ps=numofsymbolerror/N;	       
pb=numofbiterror/(2*N);   
endfunction

function [y1,y2]=bmgauss(m,sigma,N)
    if nargin == 0
        m=0; sigma=1; N=1;
    elseif nargin == 1
        sigma=m; m=0; N=1;
    elseif nargin == 2
        N=1;
    endif
    u1 = rand(1,N);
    u2 =rand(1,N);
    r = sigma*sqrt(-2*log(u1)); % 瑞利分佈
    y1 = m+r.*cos(2*pi*u2);
    y2 = m+r.*sin(2*pi*u2);
endfunction

function y=Qfunc(x)		
    y=(1/2)*erfc(x/sqrt(2));
endfunction

%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%

SNRindB1=0:2:10;
SNRindB2=0:0.1:10;
for i=1:length(SNRindB1)
    [pb,ps]=snr2p(SNRindB1(i));    	
    simu_BER(i)=pb; 
    simu_SER(i)=ps;
endfor

for i=1:length(SNRindB2)
    %SNR=exp(SNRindB2(i)*log(10)/10);     	
    SNR=10.^(SNRindB2(i)/10);
    theo_BER(i)=Qfunc(sqrt(2*SNR)); 	
end

semilogy(SNRindB1,simu_BER,'*');
hold
semilogy(SNRindB1,simu_SER,'o');
semilogy(SNRindB2,theo_BER);
xlabel('Eb/N0 (dB)');
ylabel('Prob of Error');
legend('simulation BER','simulation SER','theoretic BER');