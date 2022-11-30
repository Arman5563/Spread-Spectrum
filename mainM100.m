clc ; clear all; close all ; 
N=200;
M = 100;
%---------------------------------------------------
% index 1 for normal and index 2 for spread spectrum 
%--------------------------------------------------- 
n1 = 1:N;
n2 = 1:N*M;
Tb = M;
wgn_power=0.3;

omega_sin=0.3.*2.*pi;% omega for sin
fsinc=0.3; % for sinc. sinc is used for low pass fiter
fm=0.3;
% for low pass filters
nsinc1=(n1-1-N./2.*1);
nsinc2=(n2-1-N./2.*M);
nsync_p=-10*N*M:1:10*N*M;
mysinc1= fsinc.*sinc(fsinc.*nsinc1);
mysinc2= fsinc.*sinc(fsinc.*nsinc2);

unit_spread=ones(1,M);%for sperading
pn=2.*randi([0,1],1,M)-1;%pn generator

initial_data=randi([0,1],1,N);

%BPSk : convert form 0 & 1 to 1 & -1
initial_data_BPSK = 2*initial_data - 1;

array1= initial_data_BPSK;
array2=transpose((transpose(initial_data_BPSK)*unit_spread)); %%for spreading before modulation

modulator1=cos(2.*pi.*fm.*n1);
modulator2=cos(2.*pi.*fm.*n2);
% modulation-------------------
modulated1=initial_data.*modulator1;
modulated2=(array2(:)').*modulator2;
%------------------------------
% PN---------------------------
modulated_p2 = modulated2.*repmat(pn,1,N);
%------------------------------
% gaussian noise and sine-------
wgn_noise1=wgn(1,length(n1),wgn_power,'linear');
sin_n1=sine_noise(M,n1,omega_sin);
noise1=sin_n1+wgn_noise1;

wgn_noise2=wgn(1,length(n2),wgn_power,'linear');
sin_n2=sine_noise(M,n2,omega_sin);
noise2=sin_n2+wgn_noise2;

%adding noise
noised_signal1=noise1+modulated1;
noised_signal2=noise2+modulated_p2;
%-------------------------------
%PN again----------------------------
pn_removed2=noised_signal2.*repmat(pn,1,N);
%-------------------------------
% demodulation -----------------
demodulated1=noised_signal1.*modulator1;
demodulated2=pn_removed2.*modulator2;
% low pass filter
filtered1=2.*conv(mysinc1,demodulated1,'same');
filtered2=2.*conv(mysinc2,demodulated2,'same');
% ------------------------------

% detector => calculate the sum and decide by the sum (accumulator block)
data1=  filtered1 ;
data2=  sum(reshape(filtered2,M,N)) ; % sum of M samples
%thereshold of the detector
thereshold1 = (max(data1) + min(data1))/2;
thereshold2 = (max(data2) + min(data2))/2;

for i=1:N
    if data1(i) <= thereshold1 
        data1(i) = -1 ;
    else
        data1(i) = 1;
    end
    if data2(i) <= thereshold2 
        data2(i) = -1 ;
    else
        data2(i) = 1;
    end
end

% converting from -1 & 1 to 0 & 1
data_out1 = (data1 + 1)/2;
data_out2 = (data2 + 1)/2;

% calculating error = abs of difference

error_normal=sum(abs(data_out1-initial_data))/N

error_spread_spectrum=sum(abs(data_out2-initial_data))/N






