clear all
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Max_Range=200;
Max_Velocity=100;
Range_Resolution=1;
c=3e8;
Range_target=110;
Velocity_target=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=77e9;
sweep_time=5.5;
B=c/(2*Range_Resolution);
Tc=(sweep_time*2*Max_Range)/c;
slope=B/Tc;
Nd=128;
Nr=1024;
t=linspace(0,Nd*Tc,Nr*Nd);

Tx=zeros(1,length(t));
Rx=zeros(1,length(t));
Mix=zeros(1,length(t));

rt=zeros(1,length(t));
td=zeros(1,length(t));
for i=1:length(t)         
    rt(i) = Range_target + (Velocity_target*t(i));
    td(i) = (2*rt(i)) / c; 
    Tx(i) = cos( 2*pi*( fc*(t(i)) + ( 0.5 * slope * t(i)*t(i))) );
    Rx(i) = cos( 2*pi*( fc*(t(i)-td(i)) + ( 0.5 * slope * (t(i)-td(i))^2) ) );
    Mix(i) = Tx(i).*Rx(i); 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mix=reshape(Mix,[Nr,Nd]);
sig_fft1=fft(Mix,Nr);
sig_fft1=sig_fft1./Nr;
sig_fft1=abs(sig_fft1);

single_fft1=sig_fft1(1:Nr/2);
figure('Name','Range From FFT');
plot(single_fft1);
axis([0 200 0 1]);
title('Range from First FFT');
ylabel('Normalized Amplitude');
xlabel('Range');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mix=reshape(Mix,[Nr,Nd]);
sig_fft2=fft2(Mix,Nr,Nd);
sig_fft2=sig_fft2(1:Nr/2,1:Nd);
sig_fft2=fftshift(sig_fft2);

rdm=abs(sig_fft2);
rdm=10*log10(rdm);

dopller=linspace(-100,100,Nd);
range_axis=linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','Range from FFT2')
surf(dopller,range_axis,rdm);
title('Amplitude and Range FFT2');
xlabel('Speed');
ylabel('Range');
zlabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tr=10;Td=8;Gr=4;Gd=4;offset=1.4;

rdm=rdm/max(max(rdm));
for i=Tr+Gr+1:(Nr/2)-(Gr+Tr)
    for j=Td+Gd+1:Nd-(Gd+Td)
        noise=zeros(1,1);
        for p=i-(Tr+Gr):i+(Tr+Gr)
            for q=j-(Td+Gd):j+(Td+Gd)
                if(abs(i-p)>Gr||abs(j-q)>Gd)
                    noise=noise+db2pow(rdm(p,q));
                end
            end
        end
        
        threshold=pow2db(noise/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
        threshold=threshold+offset;
        CUT=rdm(i,j);
        if(CUT<threshold)
            rdm(i,j)=0;
        else
            rdm(i,j)=1;
        end
    end
end
rdm(1:(Tr+Gr), :) = 0;
rdm(end-Tr-Gr:end, :) = 0;
rdm(:, 1:(Td+Gd)) = 0;
rdm(:, end-Td-Gd:end) = 0;
figure('Name','CA-CFAR RDM');
surf(dopller,range_axis,rdm);
colorbar;
title( 'CA-CFAR RDM');
xlabel('Speed');
ylabel('Range');
zlabel('Normalized Amplitude');
view(315, 45);





























