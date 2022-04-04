%% Sinüs Dalgasý

Fs=150; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
f=5; %f frekanslý Sinüs dalgasý oluþturuyoruz
x=sin(2*pi*t*f);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Sinüs Dalgasý Sinyali');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Sinüs Dalgasý Sinyali');
xlabel('Frekans(Hz)');
ylabel('Güç');
%%% Sinüs Dalgasý

%% Cosinüs Dalgasý
Fs=150; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
f=5; %f frekanslý Sinüs dalgasý oluþturuyoruz
x=cos(2*pi*t*f);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Cosinüs Dalgasý Sinyali');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Cosinüs Dalgasý Sinyali');
xlabel('Frekans(Hz)');
ylabel('Güç');
%%% Cosinüs Dalgasý

%% Cosinüs Dalgasý Faz Kaydýrma

Fs=150; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
f=5; %f frekanslý Sinüs dalgasý oluþturuyoruz
pha=1/3*pi;%faz kaydýrma
x=cos(2*pi*t*f+pha);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Cosinüs Dalgasý Sinyali ve Fazý');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Cosinüs Dalgasý Sinyali ve Fazý');
xlabel('Frekans(Hz)');
ylabel('Güç');
%%% Cosinüs Dalgasý Faz Kaydýrma

%% Kare Dalga

Fs=150; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
f=5; %f frekanslý Sinüs dalgasý oluþturuyoruz
x=square(2*pi*t*f);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Kare Dalga Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Kare Dalgasý Sinyal');
xlabel('Frekans(Hz)');
ylabel('Güç');

%%+ Kare Dalga


%% Kare Darbe

Fs=150; %Örnekleme Frekansý Sampling Frequency
t=-0.5:1/Fs:0.5; %Zaman vektörü 1 saniye
w=.2; %kare geniþliði
x=rectpuls(t,w);
nfft=512; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2);
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Kare  Darbe Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Kare Darbe Sinyal');
xlabel('Frekans(Hz)');
ylabel('Güç');

%%% Kare Darbe

%% Gaussian Darbe

Fs=60; %Örnekleme Frekansý Sampling Frequency
t=-0.5:1/Fs:0.5; %Zaman vektörü 1 saniye
x=1/(sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Gaussian  Darbe Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Gaussian Darbe ');
xlabel('Frekans(Hz)');
ylabel('Güç');

%%%Gaussian Darbe

%% Eksponansiyel Azalma

Fs=150; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
x=2*exp(-5*t);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Eksponansiyel  Azalma Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Eksponansiyel  Azalma Sinyali ');
xlabel('Frekans(Hz)');
ylabel('Güç');

%%% Eksponansiyel Azalma

%% Chirp Sinyali Cývýltý

Fs=200; %Örnekleme Frekansý Sampling Frequency
t=0:1/Fs:1; %Zaman vektörü 1 saniye
x=chirp(t,0,1,Fs/6);
nfft=1024; %FFT uzunluðu
datafft=fft(x,nfft); %%fft aldýk
datafft=datafft(1:nfft/2); %% FFT simetrik
mx=abs(datafft);
%frekans vektörü
f=(0:nfft/2-1)*Fs/nfft;
figure(1);
plot(t,x);
title('Chirp Sinyali ');
xlabel('Zaman(s)');
ylabel('Genlik');
figure(2);
plot(f,mx);
title('Güç Spektrumu Chirp Sinyali ');
xlabel('Frekans(Hz)');
ylabel('Güç');

%%



