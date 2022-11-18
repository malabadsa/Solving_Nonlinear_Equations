clear all; close all; clc;

% specify input parameter
distance = input('Enter fiber length (in units of LD) = ');
beta2 = input('dispersion: 1 for normal, -1 for anomalous ');
N = input('Nonlinear parameter N = ');
mshape = input('m = 0 for sech, m > 0 for super-Gaussian = ');
chirp0 = 0;

% set simulation parameter
nt = 1024; Tmax = 32;  %fft points and window size
step_num = round(20*distance*N^2); %no. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2*Tmax)/nt; %step size in tau

%tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau;
omega = (pi/Tmax)*[(0:nt/2-1) (-nt/2:-1)];  % frequency grid

% input field profile
if mshape ==0
    uu = sech(tau).*exp(-0.5i*chirp0*tau.^2);
else
    uu = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));
end

%plot input pulse shape and spectrum
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); %spectrum
figure;
subplot(211), plot(tau, abs(uu).^2,'--k');
hold on;
axis([-5 5 0 inf]);
xlabel('Normalized time');
ylabel('Normalized power');
title('Input and output pulse shape and spectrum');
subplot(212), plot(fftshift(omega)/(2*pi), abs(temp).^2, '--k'); hold on;
axis([-.5 .5 0 inf]);
xlabel('Normalized Frequency');
ylabel('Spectral power');

%store dispersive phase shifts to speedup code
dispersion = exp(i*0.5*beta2*omega.^2*deltaz); %phase factor
hhz = 1i*N^2*deltaz;    %nonlinear phase factor

%main loop
temp = uu.*exp(abs(uu).^2.*hhz/2);
for n = 1:step_num
    f_temp = ifft(temp).*dispersion;
    uu = fft(f_temp);
    temp = uu.*exp(abs(uu).^2.*hhz);
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); %final field
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);    %final spectrum


%plot output pulse shape and spectrum
subplot(211), plot(tau, abs(uu).^2,'-k');
subplot(212), plot(fftshift(omega)/(2*pi), abs(temp).^2, '-k');
