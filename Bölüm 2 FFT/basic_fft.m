%% Dalga Formu ve Genlik Spektrumu

dt=1/100; %Örnekleme deðeri
et=4; % Aralýk Sonu
t=0:dt:et; %Örnekleme uzunluðuð
y=3*sin(4*2*pi*t)+5*sin(2*2*pi*t); %Örneklenen Sinyal
subplot(2,1,1)
plot(t,y);
grid on
axis([0 et -8 8]); %ölçek ayarý
xlabel('Zaman(s)');
ylabel('Genlik');

%güç spektrumuna geçiþ fft
Y=fft(y);
n=size(y,2)/2; % 2nci yarý complex conjuge
amp_spec=abs(Y)/n; %kesin deðer normalize
subplot(2,1,2)
freq=(0:79)/(2*n*dt); % pencere
plot(freq,amp_spec(1:80));
grid on
xlabel('Frekans(Hz)'); % 1 Herz=1 saniyedeki örnekleme cycles/second
ylabel('Genlik');

%% sinüs modülasyonlu bir Gauss fonksiyonunun 
%% FFT'sini hesaplamak ve çalýþma süresini karþýlaþtýrmak için

fr=3200;%Modülasyon Frekansý(Hz)
dt=5; %zaman basamaðý time step
N=1024;
w=2*pi*fr;
tic
for k=1:N
    X(k)=sin(w*dt*k)*exp(-(dt*N/(10*dt*k))^2);
end
X2=fft(X,N);
X2=fftshift(X2); % saða sola yarým daðýt
X2=abs(X2);
% örnekleme frekansý ayarka
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

%% Gürültü Filtrelenmiþ Sinyali Fourier Dönüþümü

%önce Gürültü üretioruz

dt=1/100;%Örnekleme deðeri
et=4; %örnekleme sonu
t=0:dt:et; %Örnekleme Aralýðý
y=3*sin(4*2*pi*t)+5*sin(2*2*pi*t); %Örneklenen Sinyal
Y=fft(y);
noise=randn(1,size(y,2)); %random noise
ey=y+noise; %Örneklenen gürültü
eY=fft(ey);
n=size(ey,2)/2; %ölçek kullanýyoruz
amp_spec=abs(eY)/n;

figure
subplot(2,1,1)
plot(t,ey) 
grid on
axis([0 et -8 8]); %ölçek ayarý
xlabel('Zaman(s)');
ylabel('Genlik');

subplot(2,1,2)
freq=(0:79)/(2*n*dt); % pencere
plot(freq,amp_spec(1:80));
grid on
xlabel('Frekans(Hz)'); % 1 Herz=1 saniyedeki örnekleme cycles/second
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
axis([0 et -8 8]); %ölçek ayarý
xlabel('Zaman(s)');
ylabel('Genlik');


