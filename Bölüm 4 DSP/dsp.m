%%  Alçak Geçiren Filtre kulalnýmý Frekans Örnekleme Metodu

clear all;
clc;
close all;
M=63;
Wp=0.5*pi; %örnek sayýsý ve geçiþ bandý kesme frekansý
m=0:(M+1)/2;
Wm=2*pi*m./(M+1);%Örnekleme sayýsý ve durduran bant kesme frekansý
mtr=floor(Wp*(M+1)/(2*pi))+2; %negatife yuvarlama floor(3.5)=3 Floor(-3.2)=4
Ad=[Wm<=Wp];
Ad(mtr)=0.38;
Hd=Ad.*exp(-j*0.5*M*Wm); %frekans etki alaný örnekleme vektörünü tanýmlayýn H(k)
Hd=[Hd conj(fliplr(Hd(2:(M+1)/2)))];
h=real(ifft(Hd)); %h(n)=IDFT H(k)
w=linspace(0,pi,1000);
H=freqz(h,[1],w);
xlabel("Normalize Frekans")
ylabel("Kazanç/dB")
title("Kazanç Yanýtý Düþük geçiren filter")
axis([0 1 -50 0.5]);
f1=100;f2=300;f3=700;fs=2000;%the frequencies of sines signal that needs filtered and the sample frequency
figure(2)
subplot(211)
t=0:1/fs:0.25;%define the time domain and steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t);%signal before filtering
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512); AFs=abs(Fs);%transform to the frequency domain
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
sf=filter(h,1,s);%use function filter
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering')
axis([0.2 0.25 -2 2]);%set the range of image coordinates
subplot(212)
Fsf=fft(sf,512); AFsf=abs(Fsf);%frequency-domain and the amplitude diagram
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% Hamming Pencere

clear all;
clc;
close all;
m=80; % The length of the Karnel%
n=0:1:m-1; %Defines Range of position value%
p=n-(m-1)/2; %Angle%
fc=0.1; %Define Cutoff frequency%
Z=sin(2*pi*fc*p)./(pi*p); %Define truncated Sinc function%
stem(n,Z);grid %To represent the discrete signal value & draw the grid lines%
title('Unit Sample Response of the sin function') %Define the title of the figure%
xlabel('n') %Define the label of on x axis%
ylabel('Z') %Define the label of on y axis%
figure; %To draw a figure%
[h,w]=freqz(Z); %get Frequency Response%
plot(w/pi,abs(h));grid
title('Frequency response of the sin function')
xlabel('Frequency')
ylabel('Amplitude')
figure;
s=2*pi*(n/(m-1)); %Angle%
w=0.54-0.46*cos(s); %Define Hamming window function%
stem(n,w);grid
title('Hamming Window')
xlabel('n')
ylabel('w')
figure;
t=Z.*w; %Multiplication of Hamming Window and sin function%
stem(n,t);grid
title('Multiplication of Hamming Window and sin function')
xlabel('n')
ylabel('t')
figure;
[h,w]=freqz(t);
plot(w/pi,abs(h));grid
xlabel('Frequency')
ylabel('Magnitude')
title('Frequency response of the windowed sin function')
figure;grid
freqz(t)
title('Frequency Response of the windowed sin function in dB')

%%  Blackman Pencere

clear all;
clc;
close all;
m=80; % The length of the Karnel%
n=0:1:m-1; %Defines Range of position value%
p=n-(m-1)/2; %Angle%
fc=0.1; %Define Cutoff frequency%
Z=sin(2*pi*fc*p)./(pi*p); %Define truncated Sinc function%
stem(n,Z);grid %To represent the discrete signal value & draw the grid lines%
title('Unit Sample Response of the sin function') %Define the title of the figure%
xlabel('n') %Define the label of on x% axis%
ylabel('Z') %Define the label of on y% axis%
figure; %To draw a figure%
[h,w]=freqz(Z); %get Frequency Response%
plot(w/pi,abs(h));grid %Plot figure & draw the grid lines%
title('Frequency response of the sin function')
xlabel('Frequency')
ylabel('Amplitude')
figure;
s=2*pi*(n/(m-1));%Angle%
w=0.42-0.5*cos(s)+.08*cos(2*s); %Define Blackmanwindow function%
stem(n,w);grid
title('Blackman Window')
xlabel('n')
ylabel('w')
figure;
t=Z.*w; %Mulplying truncated Sinc function byBlackman window%
stem(n,t);grid
title('Multiplication of Blackman Window and sin function')
xlabel('n')
ylabel('t')
figure;
[h,w]=freqz(t);
plot(w/pi,abs(h));grid
title('frequency response of the windowed sin function')
xlabel('frequency')
ylabel('Magnitude')
figure;
freqz(t)
title('frequency Response of the windowed sin function indB')

%% Pencereleme Fonksiyon  Metod

clear all;
clc;
close all;
f1=100;f2=200;%the frequencies of sines signal that needs filtered
fs=2000;%sampling frequency
m=(0.3*f1)/(fs/2);%define tansition bandwidth
M=round(8/m);%define the window length
N=M-1;%define the order of filter
b=fir1(N,0.5*f2/(fs/2));%use the firl function to design a filter
%Input parameters are respectively the order number and the cutoff
%frequency of filter
figure(1)
[h,f]=freqz(b,1,512);%amplitude-frequency characteristic graph
plot(f*fs/(2*pi),20*log10(abs(h)))%parameters are respectively frequency and amplitude
xlabel('frequency/Hz');ylabel('gain/dB');title('The gain response of lowpass filter');
figure(2)
subplot(211)
t=0:1/fs:0.2;%define the time domain and the steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t);%signal before filtering
plot(t,s);%plot the signal graph before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
axis([0 0.1 -2 2]);
subplot(212)
Fs=fft(s,512);%transform the signal to frequency domain
AFs=abs(Fs);%take the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
sf=filter(b,1,s);%use filter function to filter
subplot(211)
plot(t,sf)%plot the signal graph after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering');
axis([0.1 0.2 -2 2]);
subplot(212)
Fsf=fft(sf,512);%frequency-domain diagram after filtering
AFsf=abs(Fsf);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram after filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% Band geçiren filtre Pencer efonksiyon metoduyla

clear all;
clc;
close all;
Rs=0.01;fs=8000;%sample frequency
fcuts=[1000 1300 2210 2410];a=[0 1 0];
dev=Rs*ones(1,length(a));
[M,Wc,beta,ftype]=kaiserord(fcuts,a,dev,fs);
%M is the minimum order of filter that meets the requirements
%Wc is cutoff frequency
b = fir1(48,[0.35 0.65]);
[h,f]=freqz(b,1,512);%amplitude-frequency characteristic diagram
%[H,W]=freqz(B,A,N) when N is an integer, function returns to N frequency
%vector and amplitude-frequency response vector
figure(1)
plot(f*fs/(2*pi),20*log10(abs(h)))% parameters are respectively frequecy and amplitude
xlabel('frequency/Hz');ylabel('gain/dB');title('The gain response of bandpass filter');
f1=500;f2=1500;f3=2000;f4=3000;%frequencies of sines signal that needs filtered
t=(0:200)/fs;%define the time steplength
t1=(0.002:0.00001:0.006);
s=sin(2*f1*pi*t)+sin(2*f2*pi*t)+sin(2*f3*pi*t)+sin(2*f4*pi*t);
s1=sin(2*f1*pi*t1)+sin(2*f2*pi*t1)+sin(2*f3*pi*t1)+sin(2*f4*pi*t1);
sf=filter(b,1,s);%use function filter
figure(2)
subplot(211)
plot(t1,s1);%plot the diagram before filtering
xlabel('time/s');
ylabel('amplitude');
title('Time-domain diagram beforefiltering');
subplot(212)
Fs=fft(s,512);
AFs=abs(Fs);f=fs/512*(0:255);
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering');
axis([0.005 0.025 -4 4]);
subplot(212)
Fsf=fft(sf,512);%frequency-domain diagram after filtering
AFsf=abs(Fsf);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram after filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% Çoklu bant geçiren filtre tasarýmý

clear all;
clc;
close all;
Rs=0.01;fs=200;%sample frequency
fcuts=[10 20 40 50 60 70 80 90];a=[0,1,0,1,0];
dev=Rs*ones(1,length(a));
[M,Wc,beta,ftype]=kaiserord(fcuts,a,dev,fs);
b = fir1(48,[0.35 0.65]);
[h,f]=freqz(b,1,512);%amplitude-frequency characteristic graph
figure(1)
plot(f*fs/(2*pi),20*log10(abs(h)))%parameters are respectively frequency and amplitude
xlabel('frequency/Hz');ylabel('gain/dB');title('The gain response of multi-passband filter');
f1=5;f2=20;f3=30;f4=55;f5=75;f6=95;%frequencies of sines signal that needs filtered
t=(0:200)/fs;%define the time steplength
s1=sin(2*f1*pi*t)+sin(2*f2*pi*t)+sin(2*f3*pi*t);
s=s1+sin(2*f4*pi*t)+sin(2*f5*pi*t)+sin(2*f6*pi*t);%the signal before filtering
sf=filter(b,1,s);%use function filter
figure(2)
subplot(211)
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512);AFs=abs(Fs);f=fs/512*(0:255);
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering');
subplot(212)
Fsf=fft(sf,512); AFsf=abs(Fsf);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram after filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% Yüksek geçiren filtre frekans örnekleme metoduyla

M=32;%the number of samples
Wp=0.6*pi;%passband cutoff frequency
m=0:M/2;%the sampling points
Wm=2*pi*m./(M+1);%stopband cutoff frequency
mtr=ceil(Wp*(M+1)/(2*pi));%round to positive part,i.e.ceil(3.5)=4;ceil(-3.2)=-3;
Ad=[Wm>=Wp];
Ad(mtr)=0.28;
Hd=Ad.*exp(-j*0.5*M*Wm);%define frequency-domain sampling vector H(k))
Hd=[Hd conj(fliplr(Hd(2:M/2+1)))];
%fliplr is to realize the fliplr of matrix and conj is the conjugate
h=real(ifft(Hd));%h(n)=IDFT[H(k)]
w=linspace(0,pi,1000);%get 1000 row vectors between 0 and pi 
H=freqz(h,[1],w);%the amplitude -frequency characteristic diagram of the filter
figure(1)
plot(w/pi,20*log10(abs(H)));%parameters are respectively the normalized frequency and amplitude
xlabel('the normailzed frequency');ylabel('gian/dB');title('The gain response of highpass filter');
axis([0 1 -50 0]);
f1=200;f2=700;f3=800;%the frequencies of sines signal that needs filtered
fs=2000;%the sample frequency
figure(2)
subplot(211)
t=0:1/fs:0.25;%define the time domain and steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t);%signal before filtering
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512);%transform to the frequency domain
AFs=abs(Fs);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
sf=filter(h,1,s);%use function filter
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering')
axis([0.2 0.25 -2 2]);%set the range of image coordinates
subplot(212)
Fsf=fft(sf,512);AFsf=abs(Fsf);
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% funciton remez kullanarak eþ dalgalý alçak geçiren filtre tasarlama


fs=2000;%the sample frequency40
rp=3;%passband ripple
rs=40;%stopband ripple
f=[500 600];%cutoff frequency
a=[1 0];%desired amplitude
dev=[(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
[n,fo,ao,w]=remezord(f,a,dev,fs);
%function remezord returns to parameter n and it stand for the order of
%filter
%when the number of frequency band is B, f,a,dev are vectors with 2B-2,B,B
%elements respectively.
b=remez(n,fo,ao,w);%The return number of function Remez is the order of n-order filter.
%fo,ao are vectors with 2B elements, represent boundary frequencies
%and amplitude in B frequencies band.
%w is a vector with B elements, represents weighted value of each
%frequency band.
figure(1)
freqz(b,1,1024,fs);%the characteristic diagram of the filter
f1=400;f2=700;%the frequencies of sines signal that needs filtered
t=0:1/fs:0.1;%define the time domain and steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t);%signal before filtering
figure(2)
subplot(211)
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512);%transform to the frequency domain
AFs=abs(Fs);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
figure(3)
sf=filter(b,1,s);%use function filter
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering')
subplot(212)
Fsf=fft(sf,512);AFsf=abs(Fsf);41
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% funciton remez kullanarak eþ dalgalý bant geçiren filtre tasarlama

fs=2000;%the sample frequency
rp=3;%passband ripple
rs=40;%stopband ripple
f=[200 300 600 700];%cutoff frequency
a=[0 1 0];%desired amplitude
dev=[10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
[n,fo,ao,w]=remezord(f,a,dev,fs);
b=remez(n,fo,ao,w);
figure(1)
freqz(b,1,1024,fs);%the characteristic diagram of the filter
f1=100;f2=400;f3=500;f4=800;%the frequencies of sines signal that needs filtered
t=0:1/fs:0.1;%define the time domain and steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t)+sin(2*pi*f4*t);%signal before filtering
figure(2)
subplot(211)
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512);%transform to the frequency domain
AFs=abs(Fs);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
sf=filter(b,1,s);%use function filter
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering')
subplot(212)
Fsf=fft(sf,512);AFsf=abs(Fsf);
f=(0:255)*fs/512;%frequency sampling42
plot(f,AFsf(1:256))%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% funciton remez kullanarak eþ dalgalý bant durduran filtre tasarlama

fs=2000;%the sample frequency
rp=3;%passband ripple
rs=40;%stopband ripple
f=[200 300 600 700];%cutoff frequency
a=[1 0 1];%desired amplitude
dev=[10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
[n,fo,ao,w]=remezord(f,a,dev,fs);
b=remez(n,fo,ao,w);
figure(1)
freqz(b,1,1024,fs);%the characteristic diagram of the filter
f1=100;f2=400;f3=500;f4=800;%the frequencies of sines signal that needs filtered
t=0:1/fs:0.1;%define the time domain and steplength
s=sin(2*pi*f1*t)+sin(2*pi*f2*t)+sin(2*pi*f3*t)+sin(2*pi*f4*t);%signal before filtering
figure(2)
subplot(211)
plot(t,s);%plot the diagram before filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram before filtering');
subplot(212)
Fs=fft(s,512);%transform to the frequency domain
AFs=abs(Fs);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram before filtering');
figure(3)
sf=filter(b,1,s);%use function filter
subplot(211)
plot(t,sf)%plot the diagram after filtering
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram after filtering')
subplot(212)
Fsf=fft(sf,512);AFsf=abs(Fsf);
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram before filtering43
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domain diagram after filtering');

%% FIR filtre tasarýmý

clear all;
close all;
clc;

Fs=8e3; %Specify Sampling Frequency
Ts=1/Fs; %Sampling period.
Ns=512; %Number of time samples to be plotted.
t=[0:Ts:Ts*(Ns-1)]; %Make time array that contains Ns elements
 %t = [0, Ts, 2Ts, 3Ts,..., (Ns-1)Ts]
f1=500;
f2=1800;
f3=2000;
f4=3200;
x1=sin(2*pi*f1*t); %create sampled sinusoids at different frequencies
x2=sin(2*pi*f2*t);
x3=sin(2*pi*f3*t);
x4=sin(2*pi*f4*t);
x=x1+x2+x3+x4; %Calculate samples for a 4-tone input signal
grid on;
N=16; %FIR1 requires filter order (N) to be EVEN
 %when gain = 1 at Fs/2.
W=[0.4 0.6]; %Specify Bandstop filter with stop band between
 %0.4*(Fs/2) and 0.6*(Fs/2)
B=fir1(N,W,'DC-1');%Design FIR Filter using default (Hamming window.
B %Leaving off semi-colon causes contents of
 %B (the FIR coefficients) to be displayed.
A=1; %FIR filters have no poles, only zeros.
freqz(B,A); %Plot frequency response - both amp and phase response.
pause; %User must hit any key on PC keyboard to go on.
figure; %Create a new figure window, so previous one isn't lost.
subplot(2,1,1); %Two subplots will go on this figure window.
Npts=200;
plot(t(1:Npts),x(1:Npts)) %Plot first Npts of this 4-tone input signal
title('Time Plots of Input and Output');
xlabel('time (s)');
ylabel('Input Sig');
 %Now apply this filter to our 4-tone test sequence
y = filter(B,A,x);
subplot(2,1,2); %Now go to bottom subplot.
plot(t(1:Npts),y(1:Npts)); %Plot first Npts of filtered signal.
xlabel('time (s)');
ylabel('Filtered Sig');
pause;
figure; %Create a new figure window, so previous one isn't lost.
subplot(2,1,1);
xfftmag=(abs(fft(x,Ns))); %Compute spectrum of input signal.
xfftmagh=xfftmag(1:length(xfftmag)/2);
 %Plot only the first half of FFT, since second half is mirror imag
 %the first half represents the useful range of frequencies from
 %0 to Fs/2, the Nyquist sampling limit.
f=[1:1:length(xfftmagh)]*Fs/Ns; %Make freq array that varies from
 %0 Hz to Fs/2 Hz.
plot(f,xfftmagh); %Plot frequency spectrum of input signal
title('Input and Output Spectra');
xlabel('freq (Hz)');
ylabel('Input Spectrum');
subplot(2,1,2);
yfftmag=(abs(fft(y,Ns)));
yfftmagh=yfftmag(1:length(yfftmag)/2);
 %Plot only the first half of FFT, since second half is mirror image
 %the first half represents the useful range of frequencies from
 %0 to Fs/2, the Nyquist sampling limit.
plot(f,yfftmagh); %Plot frequency
xlabel('freq (Hz)');
