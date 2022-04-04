%% DFT Katlanýr Özellik 

n=0:10;
x=10*(0.8).^n;
y=x(mod(-n,11)+1);
subplot(2,1,1)
stem(n,x); grid on
title('Orjinal Dizi')
xlabel('n')
ylabel('x(n)')
subplot(2,1,2)
stem(n,y); grid on
title('Katlanmýþ Dizi')
xlabel('n')
ylabel('y(n)')

X=fft(x);
Y=fft(y);

figure
k=0:10;
subplot(2,1,1)
stem(k,abs(X)); grid on
title('Genlik Spektrumu Orjinal Dizi')
xlabel('k')
ylabel('X(k)')
subplot(2,1,2)
stem(k,abs(Y)); grid on
title('Genlik Spekturmu Katlanmýþ Dizi')
xlabel('k')
ylabel('Y(k)')

figure
subplot(2,1,1)
stem(k,real(X)); grid on
title('Gerçek Komponent  X(k)')
xlabel('k')
ylabel('Re{X(k)}')
subplot(2,1,2)
stem(k,real(Y)); grid on
title('Gerçek Komponent Y(k)')
xlabel('k')
ylabel('Re{Y(k)}')


figure
subplot(2,1,1)
stem(k,imag(X)); grid on
title('Sanal Komponent  X(k)')
xlabel('k')
ylabel('Im{X(k)}')
subplot(2,1,2)
stem(k,imag(Y)); grid on
title('Sanal Komponent Y(k)')
xlabel('k')
ylabel('Im{Y(k)}')

%% Dairesel Konvolüsyon

a=input('dizi sayýsý girin x(n)');
b=input('dizi sayýsý girin h(n)');
n1=length(a);
n2=length(b);
N=max(n1,n2);
x=[a zeros(1,(N-n1))];
for i=1:N
    k=i;
    for j=1:n2
        H(i,j)=x(k)*b(j);
        k=k-1;
        if (k==0)
            k=N;
        end
    end
end
y=zeros(1,N);
M=transpose(H); % transpoz
for j=1:N
    for i=1:n2
        y(j)=M(i,j)+y(j);
    end
end
disp('y(n) çýkýþ dizisi=');
disp(y);
stem(y);
title('Dairesel Konvolüsyon çýkýþ')
xlabel('n')
ylabel('y(n)')

%% Blok Konvolüsyon üst üste gelme

x=input('dizi sayýsý girin x(n)');
h=input('dizi sayýsý girin h(n)');
n1=length(x);
n2=length(h);
N=n1+n2-1;
h1=[h zeros(1,(N-n1))];
n3=length(h1);
y=zeros(1,N);
x1=[zeros(1,n3-n2) x zeros(1,n3)];
H=fft(h1);
for i=1:n2:N
    y1=x1(i:i+(2*(n3-n2)));
    y2=fft(y1);
    y3=y2.*H;
    y4=round(ifft(y3));
    y(i:(i+n3-n2))=y4(n2:n3);
end
disp('Çýkýþ Dizisi y(n)')
disp(y(1:N));
stem(y(1:N));
title('Üst Üste Getirme Metodu')
xlabel('n')
ylabel('y(n)')












