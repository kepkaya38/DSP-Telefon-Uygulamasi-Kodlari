%% DFT

am=4;%Genlik
fm=100;%Frekans
fs=32*fm;
t=0:1/fs:1;
signal=am*sin(2*pi*fm*t); %Sinyalimiz
N=3201;
tic

signal=[signal zeros(1,N-length(signal))]; % N tane sýfýr ekledik
sigF=zeros(1,N); %toplam sayaçlarýný resetledik
for k=0:N-1
    for n=0:N-1
        sigF(k+1)=sigF(k+1)+signal(n+1)*exp(-1i*2*pi*k*n/N);
    end
end
siF=fftshift(abs(sigF));
fh=linspace(-fs/2,fs/2,length(siF));
plot(fh,siF)
title('DFT Sinüsoidal')
xlabel('Frekans(Hz)')
ylabel('Genlik')

%% DTFT

pp=-4*pi:8*pi/511:4*pi;
num=[2 1];
den=[1 -0.6];
X=freqz(num,den,pp);
subplot(2,1,1)
plot(pp/pi,real(X))
title('Gerçek Kýsým X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')
subplot(2,1,2)
plot(pp/pi,imag(X))
title('Sanal Kýsým X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')

figure
subplot(2,1,1)
plot(pp/pi,abs(X))
title('Güç Spektrumu X(\omega')
xlabel('\omega / \pi')
ylabel('Genlik')
subplot(2,1,2)
plot(pp/pi,angle(X))
title('Faz Spektrumu X(\omega')
xlabel('\omega / \pi')
ylabel('Faz Derecesi')

%% DTFT

L=input('Dizi Boyutunu girin=');
N=input('DFT uzunluðunu girin=');
%L burda zaman domeni
x=[ones(1,L)];
% N noktalý DFT
X=fft(x,N);
clf;
t=0:1:L-1;
stem(t,x)
title('Orjinal Zaman Dizisi')
xlabel('Zaman Ýndeksi n')
ylabel('Genlik')

figure
subplot(2,1,1)
k=0:1:N-1;
stem(k,abs(X))
title('DFT Örnek Gücü')
xlabel('Frekans Ýndeksi k')
ylabel('Genlik')
subplot(2,1,2)
stem(k,angle(X))
title('DFT Örnek Açýsý')
xlabel('Frekans Ýndeksi k')
ylabel('Fazý')

%% IDFT

K=input('DFT uzunluðu girin=');
L=input('IDFT uzunluðunu girin=');
% K noktalý DFT
k=1:K;
X=(k-1)/K;
x=ifft(X,L);
clf;
l=1:K;
stem(k-1,X)
xlabel('Frekans Ýndeksi k')
ylabel('Genlik')
title('Orijinal DFT Örneklemesi')

figure
subplot(2,1,1)
n=0:1:L-1;
stem(n,real(x))
title('Gerçek Kýsým Zaman Domeni Örneklemesi')
xlabel('Zaman Ýndeksi n')
ylabel('Genlik')
subplot(2,1,2)
stem(n,imag(x))
title('Sanal Kýsým Zaman Domeni Örneklemesi')
xlabel('Zaman Ýndeksi n')
ylabel('Fazý')







