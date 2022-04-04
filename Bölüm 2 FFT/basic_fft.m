%% Dalga Formu ve Genlik Spektrumu

dt=1/100; %�rnekleme de�eri
et=4; % Aral�k Sonu
t=0:dt:et; %�rnekleme uzunlu�u�
y=3*sin(4*2*pi*t)+5*sin(2*2*pi*t); %�rneklenen Sinyal
subplot(2,1,1)
plot(t,y);
grid on
axis([0 et -8 8]); %�l�ek ayar�
xlabel('Zaman(s)');
ylabel('Genlik');

%g�� spektrumuna ge�i� fft
Y=fft(y);
n=size(y,2)/2; % 2nci yar� complex conjuge
amp_spec=abs(Y)/n; %kesin de�er normalize
subplot(2,1,2)
freq=(0:79)/(2*n*dt); % pencere
plot(freq,amp_spec(1:80));
grid on
xlabel('Frekans(Hz)'); % 1 Herz=1 saniyedeki �rnekleme cycles/second
ylabel('Genlik');

%% sin�s mod�lasyonlu bir Gauss fonksiyonunun 
%% FFT'sini hesaplamak ve �al��ma s�resini kar��la�t�rmak i�in

fr=3200;%Mod�lasyon Frekans�(Hz)
dt=5; %zaman basama�� time step
N=1024;
w=2*pi*fr;
tic
for k=1:N
    X(k)=sin(w*dt*k)*exp(-(dt*N/(10*dt*k))^2);
end
X2=fft(X,N);
X2=fftshift(X2); % sa�a sola yar�m da��t
X2=abs(X2);
% �rnekleme frekans� ayarka
fmax=1/dt; 
df=1/(N*dt);
for k=1:N
    F(k)=-fmax/2+(k-1)*df;
end

plot(F,X2);
title('FFT zaman Sinyali')
xlabel('Frekans(Hz)')
ylabel('Genlik')
toc

%% G�r�lt� Filtrelenmi� Sinyali Fourier D�n���m�

%�nce G�r�lt� �retioruz

dt=1/100;%�rnekleme de�eri
et=4; %�rnekleme sonu
t=0:dt:et; %�rnekleme Aral���
y=3*sin(4*2*pi*t)+5*sin(2*2*pi*t); %�rneklenen Sinyal
Y=fft(y);
noise=randn(1,size(y,2)); %random noise
ey=y+noise; %�rneklenen g�r�lt�
eY=fft(ey);
n=size(ey,2)/2; %�l�ek kullan�yoruz
amp_spec=abs(eY)/n;

figure
subplot(2,1,1)
plot(t,ey) 
grid on
axis([0 et -8 8]); %�l�ek ayar�
xlabel('Zaman(s)');
ylabel('Genlik');

subplot(2,1,2)
freq=(0:79)/(2*n*dt); % pencere
plot(freq,amp_spec(1:80));
grid on
xlabel('Frekans(Hz)'); % 1 Herz=1 saniyedeki �rnekleme cycles/second
ylabel('Genlik');

figure
plot(Y/n,'r+')
hold on
plot(eY/n,'bx')

figure
fY=fix(eY/100)*100;
ifY=ifft(fY);
cy=real(ifY);
plot(t,cy) 
grid on
axis([0 et -8 8]); %�l�ek ayar�
xlabel('Zaman(s)');
ylabel('Genlik');


