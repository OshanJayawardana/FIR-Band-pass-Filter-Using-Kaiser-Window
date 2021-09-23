%defining parameters
A=4;
B=3;
C=7;
L=200;

Ap=0.03+(0.01*A)
Aa=45+B
Op1=C*100+300
Op2=C*100+700
Oa1=C*100+150
Oa2=C*100+800
Ws=2*(C*100+1200)
Bt=min([Oa2-Op2,Op1-Oa1])
Wc1=Op1-Bt/2;
Wc1
Wc2=Op2+Bt/2;
Wc2
T=2*pi/Ws
x=-L:L;
%***********************************************************************************%

%plotting ideal filter impulse response
y=hnT(x,Wc1,Wc2,Ws,T); %calling the ideal filter function
f1=figure;
stem(x,y);
xlabel("n")
ylabel("Amplitude h_i_d_e_a_l[n]")
title("Ideal filter: Impulse respomse")
saveas(f1,'1.png')
%***********************************************************************************%

%calculating specifications of the kaiser window
deltap=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1); 
deltaa=10^(-0.05*48);
delta=min([deltap deltaa]); %finding delta
delta
Aaf=(-20)*log10(delta)      %calculating the actual value for Aa
alfa=findalfa(Aaf);         %calculating alfa
alfa
D=findD(Aaf);               %calculating D
D
Ntemp=(Ws*D/Bt)+1;
Nint=round(Ntemp);          %calculating N
if Nint<Ntemp
    Nint=Nint+1;
end
if rem(Nint,2)==0
    Nint=Nint+1;
end
Nint
M=Nint-1;
Tou=M/2;                    %calculating Tou
%***********************************************************************************%

wkn=wknT(x,alfa,Tou,Nint);  % calling kaiser window function
f2=figure;
stem(x,wkn)                 % plotting kaiser window impulse response
xlabel("n")
ylabel("Amplitude w[n]")
title("Kaiser winndow: Impulse response")
saveas(f2,'2.png')
%***********************************************************************************%

hn=y.*wkn; % windowed filter function
f3=figure;
stem(x,hn) %ploting filter
xlabel("n")
ylabel("Amplitude h[n]")
title("Windowed filter: Impulse response")
saveas(f3,'3.png')
%***********************************************************************************%

%shifting filter to make it causal
hdash=zeros(size(hn));
for c3=1:(length(hn)-Tou)
    hdash(c3+Tou)=hn(c3);
end
f4=figure;
stem(x,hdash)

xlim([0 Nint]) %plotting causal filter
xlabel("n")
ylabel("Amplitude h'[n]")
title("Non-recursive FIR filter: Impulse response")
saveas(f4,'4.png')
%***********************************************************************************%

[h,w]=freqz(hdash,1); % frequency response of filter
f5=figure;
plot(w/pi,db(abs(h)))
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Magnitude (dB) H'(w)")
title("Non-recursive FIR filter: Magnitude response")
saveas(f5,'5.png')

f51=figure;
plot(w/pi,abs(h))
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Magnitude H'(w)")
title("Non-recursive FIR filter: Magnitude response")
saveas(f51,'5.1.png')
%***********************************************************************************%

%plotting magnitude response of pass band
f6=figure;
plot(w/pi,db(h))
xlim([Op1*2/Ws Op2*2/Ws])
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Magnitude (dB) H'(w)")
title("Non-recursive FIR filter: Magnitude of the passband")
saveas(f6,'6.png')
%***********************************************************************************%

%creating the excitation
xn=sin(x*Oa1/2*T)+sin(x*(Op1+Op2)/2*T)+sin(x*(Oa2+Ws/2)/2*T);
yfilt=conv(hdash,xn,"same"); %appying the filter to the signal
f7=figure;
stem(x,yfilt);
xlabel("n")
ylim([-1.5 1.50])
xlim([-200 -100])
ylabel("Amplitude y[n]")
title("Filtered signal using the new filter")
saveas(f7,'7.png')
%***********************************************************************************%

yfiltid=conv(hnT(x,Wc1,Wc2,Ws,T),xn,"same"); %applying the ideal filter to the signal
f8=figure;
stem(x,yfiltid)
xlim([-200 -100])
xlabel("n")
ylabel("Amplitude y_i_d_e_a_l[n]")
title("Filtered signal using ideal filter")
saveas(f8,'8.png')
%***********************************************************************************%

f=(0:2*L)/L; %creating the x axis for the frequency domain
filtfft=fft(hdash); %taking DTFT of new filter
infft=fft(xn); %taking DTFT of input signal
f9=figure;
yyaxis left
plot(f,abs(infft)) %plottting magnitude response of input signal
ylabel("Amplitude : magnitude response of the input signal X(w)")
hold on
yyaxis right
plot(f,abs(filtfft)) %plotting magnitude response of new filter
ylabel("Amplitude : magnitude response of the new filter H'(w)")
hold off
xlim([0 1])
xlabel("Normalized Frequency (\times\pi rad/sample)")
title("Filter and the input signal: Magnitude response")
saveas(f9,'9.png')
%***********************************************************************************%

f10=figure;
outfft=fft(yfilt); %taking DTFT of filtered signal
plot(f,abs(outfft)) %plotting magnitude response of filtered signal
ylabel("Amplitude Y(w)")
xlabel("Normalized Frequency (\times\pi rad/sample)")
title("Filtered signal: Magnitude response")
xlim([0 1])
saveas(f10,'10.png')
%***********************functions***********************************%
function [alfa] = findalfa(Aaf) % calculate alfa
    if Aaf<=21
        alfa=0;
    else 
        if Aaf>50
            alfa = 0.1102*(Aaf - 8.7);
        else
            alfa=0.5842*(Aaf-21)^0.4+0.07886*(Aaf-21);
        end
    end
end
%*******************************************************************%

function [h] = hnT(n,Wc1,Wc2,Ws,T) % ideal filter function
h=zeros(size(n));
   for c1 = 1:length(n)
       if n(c1)==0
           h(c1)=2*(Wc2-Wc1)/Ws;
       else
           h(c1)=(sin(n(c1)*Wc2*T)-sin(Wc1*n(c1)*T))/(n(c1)*pi);
       end
   end
   %if n==0
       %h=2*(Wc2-Wc1)/Ws;
   %else
       %h=(sin(n*Wc2*T)-sin(Wc1*n*T))
   %end
end
%*******************************************************************%

function [D] = findD(Aaf) % calculate D
if Aaf<=21
    D=0.9222;
else
    D=(Aaf-7.95)/14.36;
end
end
%*******************************************************************%

function [w] = wknT(n,alfa,Tou,Nint) % kaiser window function
beta=alfa*(1-(2*n/(Nint-1)).^2).^0.5;
w=Io(beta)/Io(alfa);
for c2 = 1:length(n)
   if (-Tou>=n(c2))||(n(c2)>= Tou)
       w(c2)=0;
   end    
end
end
%*******************************************************************%

function [i] = Io(x) %bessel function 
i=ones(size(x));
for k = 1:10
    i=i+(((x/2).^k)*(1/factorial(k))).^2;
end
end
%********************************************************************%
