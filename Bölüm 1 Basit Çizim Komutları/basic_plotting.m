%% Sin�s Dalgas�

Fs=150; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
f=5; %f frekansl� Sin�s dalgas� olu�turuyoruz
x=sin(2*pi*t*f);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Sin�s Dalgas� Sinyali');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Sin�s Dalgas� Sinyali');
xlabel('Frekans(Hz)');
ylabel('G��');
%%% Sin�s Dalgas�

%% Cosin�s Dalgas�
Fs=150; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
f=5; %f frekansl� Sin�s dalgas� olu�turuyoruz
x=cos(2*pi*t*f);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Cosin�s Dalgas� Sinyali');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Cosin�s Dalgas� Sinyali');
xlabel('Frekans(Hz)');
ylabel('G��');
%%% Cosin�s Dalgas�

%% Cosin�s Dalgas� Faz Kayd�rma

Fs=150; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
f=5; %f frekansl� Sin�s dalgas� olu�turuyoruz
pha=1/3*pi;%faz kayd�rma
x=cos(2*pi*t*f+pha);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Cosin�s Dalgas� Sinyali ve Faz�');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Cosin�s Dalgas� Sinyali ve Faz�');
xlabel('Frekans(Hz)');
ylabel('G��');
%%% Cosin�s Dalgas� Faz Kayd�rma

%% Kare Dalga

Fs=150; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
f=5; %f frekansl� Sin�s dalgas� olu�turuyoruz
x=square(2*pi*t*f);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Kare Dalga Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Kare Dalgas� Sinyal');
xlabel('Frekans(Hz)');
ylabel('G��');

%%+ Kare Dalga


%% Kare Darbe

Fs=150; %�rnekleme Frekans� Sampling Frequency
t=-0.5:1/Fs:0.5; %Zaman vekt�r� 1 saniye
w=.2; %kare geni�li�i
x=rectpuls(t,w);
nfft=512; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Kare  Darbe Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Kare Darbe Sinyal');
xlabel('Frekans(Hz)');
ylabel('G��');

%%% Kare Darbe

%% Gaussian Darbe

Fs=60; %�rnekleme Frekans� Sampling Frequency
t=-0.5:1/Fs:0.5; %Zaman vekt�r� 1 saniye
x=1/(sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Gaussian  Darbe Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Gaussian Darbe ');
xlabel('Frekans(Hz)');
ylabel('G��');

%%%Gaussian Darbe

%% Eksponansiyel Azalma

Fs=150; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
x=2*exp(-5*t);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Eksponansiyel  Azalma Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Eksponansiyel  Azalma Sinyali ');
xlabel('Frekans(Hz)');
ylabel('G��');

%%% Eksponansiyel Azalma

%% Chirp Sinyali C�v�lt�

Fs=200; %�rnekleme Frekans� Sampling Frequency
t=0:1/Fs:1; %Zaman vekt�r� 1 saniye
x=chirp(t,0,1,Fs/6);
nfft=1024; %FFT uzunlu�u
datafft=fft(x,nfft); %%fft ald�k
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vekt�r�
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Chirp Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('G�� Spektrumu Chirp Sinyali ');
xlabel('Frekans(Hz)');
ylabel('G��');

%%



