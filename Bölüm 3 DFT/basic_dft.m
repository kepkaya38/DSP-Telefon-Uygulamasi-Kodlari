%% DFT

am=4;%Genlik
fm=100;%Frekans
fs=32*fm;
t=0:1/fs:1;
signal=am*sin(2*pi*fm*t); %Sinyalimiz
N=3201;
tic

signal=[signal zeros(1,N-length(signal))]; % N tane s�f�r ekledik
sigF=zeros(1,N); %toplam saya�lar�n� resetledik
for k=0:N-1
    for n=0:N-1
        sigF(k+1)=sigF(k+1)+signal(n+1)*exp(-1i*2*pi*k*n/N);
    end
end
siF=fftshift(abs(sigF));
fh=linspace(-fs/2,fs/2,length(siF));
plot(fh,siF)
title('DFT Sin�soidal')
xlabel('Frekans(Hz)')
ylabel('Genlik')

%% DTFT

pp=-4*pi:8*pi/511:4*pi;
num=[2 1];
den=[1 -0.6];
X=freqz(num,den,pp);
subplot(2,1,1)
plot(pp/pi,real(X))
title('Ger�ek K�s�m X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')
subplot(2,1,2)
plot(pp/pi,imag(X))
title('Sanal K�s�m X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')

figure
subplot(2,1,1)
plot(pp/pi,abs(X))
title('G�� Spektrumu X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')
subplot(2,1,2)
plot(pp/pi,angle(X))
title('Faz Spektrumu X(\omega')
xlabel('\omega / \pi')
ylabel('Faz Derecesi')

%% DTFT

L=input('Dizi Boyutunu girin=');
N=input('DFT uzunlu�unu girin=');
%L burda zaman domeni
x=[ones(1,L)];
% N noktal� DFT
X=fft(x,N);
clf;
t=0:1:L-1;
stem(t,x)
title('Orjinal Zaman Dizisi')
xlabel('Zaman �ndeksi n')
ylabel('Genlik')

figure
subplot(2,1,1)
k=0:1:N-1;
stem(k,abs(X))
title('DFT �rnek G�c�')
xlabel('Frekans �ndeksi k')
ylabel('Genlik')
subplot(2,1,2)
stem(k,angle(X))
title('DFT �rnek A��s�')
xlabel('Frekans �ndeksi k')
ylabel('Faz�')

%% IDFT

K=input('DFT uzunlu�u girin=');
L=input('IDFT uzunlu�unu girin=');
% K noktal� DFT
k=1:K;
X=(k-1)/K;
x=ifft(X,L);
clf;
l=1:K;
stem(k-1,X)
xlabel('Frekans �ndeksi k')
ylabel('Genlik')
title('Orijinal DFT �rneklemesi')

figure
subplot(2,1,1)
n=0:1:L-1;
stem(n,real(x))
title('Ger�ek K�s�m Zaman Domeni �rneklemesi')
xlabel('Zaman �ndeksi n')
ylabel('Genlik')
subplot(2,1,2)
stem(n,imag(x))
title('Sanal K�s�m Zaman Domeni �rneklemesi')
xlabel('Zaman �ndeksi n')
ylabel('Faz�')







