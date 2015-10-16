clear;
s='C:\Users\PeterHu-T420\Desktop\ODE45_DensityMatrix_Gauss_withoutchirp_xi=4';
sp=strcat(s,'\relative_parameters.dat');
collection_effect=input('No.1 considers collection effect, No.0 not:\n');

fid=fopen(sp,'r');
omega0=fread(fid,1,'double');
omega1=fread(fid,1,'double');
Omega1=fread(fid,1,'double');
mu=fread(fid,1,'double');
len=fread(fid,1,'unsigned long');
dt=fread(fid,1,'double');
T=fread(fid,1,'double');
fclose(fid);
wmg=2*pi/dt;
   

s0=strcat(s,'\source.dat');
fid0=fopen(s0,'r');
source=fread(fid0,[2 len],'double');
t=source(1,:);
source=source(2,:);
fclose(fid0);
figure
plot(t/T,source,'b-','linewidth',2);
title('Laser field');
xlabel('Time(In Optical Period )','fontsize',14);
ylabel('\it{E(t)}(a.u.)','fontsize',14);

s1=strcat(s,'\dipole.dat');
fid1=fopen(s1,'r');
dipole=fread(fid1,[1,len], 'double');
if collection_effect==1
    dipole=dipole*7.5E24*(5.29E-11)^3;
end
dipole=dipole(1:len-1);
len2=length(dipole);
t1=t(1:len2);
figure
plot(t1/T,dipole,'r-','linewidth',2)
title('dipole');
xlabel('Time(In Optical Period)','fontsize',14);
ylabel('Dipole Intensity (a.u.)','fontsize',14);

fre=(0:round(len2/2)-1)/len2*wmg/omega1;
FFA=fft(dipole);
hff=FFA;
FFA=abs(FFA(1:round(len2/2))/length(dipole));%
FFA=2*log10(abs(FFA)+eps);
figure
plot(fre,FFA,'r-','linewidth',2);
title('HHG');
xlabel('Harmonic Order(\omega/\omega_L)','fontsize',14);
ylabel('Harmonic Intensity(arb.unit)','fontsize',14);
fclose(fid1);

% s3=strcat(s,'\population_inversion.dat');
% fid3=fopen(s3,'r');
% pop_inv=fread(fid3,[1, len],'double');
% figure
% plot(t/T,pop_inv,'r-.','linewidth',2);
% xlabel('Time(In Optical Period)','fontsize',14);
% ylabel('Population Inversion','fontsize',14);
% fclose(fid3);
% 
% s4=strcat(s,'\dress_down.dat');
% s5=strcat(s,'\dress_up.dat');
% s6=strcat(s,'\dress_inv.dat');
% fid4=fopen(s4,'r');
% fid5=fopen(s5,'r');
% fid6=fopen(s6,'r');
% dress_down=fread(fid4,[1,len],'double');
% dress_up=fread(fid5,[1,len],'double');
% dress_inv=fread(fid6,[1,len],'double');
% figure
% subplot(3,1,1)
% plot(t/T,dress_down,'r-.','linewidth',2);
% % xlabel('Time(In Optical Period)','fontsize',14);
% ylabel('|C_1^A(t)|^2','fontsize',14);
% subplot(3,1,2)
% plot(t/T,dress_up,'r-.','linewidth',2);
% % xlabel('Time(In Optical Period)','fontsize',14);
% ylabel('|C_2^A(t)|^2','fontsize',14);
% subplot(3,1,3)
% plot(t/T,dress_inv,'r-.','linewidth',2);
% xlabel('Time(In Optical Period)','fontsize',14);
% ylabel('dress_-inv','fontsize',14);

s7=strcat(s,'\dipoleawt.txt');
cwtd=textread(s7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%this parameters should be same to the file wavelet.cpp
s8=strcat(s,'\wtplot.dat');
fid8=fopen(s8,'r');
nt=fread(fid8,1,'int');
mw=fread(fid8,1,'int');
tcenter=fread(fid8,[1 nt],'double');
omg=fread(fid8,[1 mw],'double');
cwtd=cwtd(:,1:mw);
flipud(cwtd);
fliplr(tcenter);
figure
[tcenter, omg]=meshgrid(tcenter,omg);
%  surf(tcenter,omg,log10(abs(cwtd)'+eps));shading('interp');
surf(tcenter,omg,abs(cwtd'));shading('interp');
%  imagesc(abs(cwtd));shading('interp');

view(0,90)

%合成脉冲

min_order=input('选择谐波最小级次：');   %合成脉冲选择谐波级次范围，根据实际情况确定
max_order=input('选择谐波最大级次：');
fre=(0:len2-1)/len2*wmg/omega1;
order_selected=find(fre<min_order|fre>max_order);
hff(order_selected)=0;
hpulse=abs(ifft(hff));
hpulse=hpulse.*hpulse;
pulse_intensity=hpulse;
figure
plot(t(1:length(pulse_intensity))/T,pulse_intensity,'-.b','linewidth',2);%
xlabel('Time(In Optical Period)','fontsize',14);
ylabel('Pulse Intensity (arb.units)','fontsize',14);
%计算合成脉冲的FWHM
t2=t(1:length(pulse_intensity));
Max_pulseintensity=max(pulse_intensity);
HalfMax_order=find(abs(pulse_intensity-Max_pulseintensity/2.0)<1e-3&abs(pulse_intensity-Max_pulseintensity/2.0)/Max_pulseintensity<1e-3);
FWHM=abs(2*t2(HalfMax_order(1))*2.41888e-17); %单位，秒
display('合成脉冲的宽度：(阿秒as)');
FWHM=FWHM*1e18  %单位。阿秒
