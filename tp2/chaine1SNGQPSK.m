% Transmission QPSK sur frequence porteuse

clear all
close all
% clc

% Initialisation
Nbits = 1000 ; % Nombre de bits emis
M = 4 ; % Taille de la constellation
n = log2(M) ; % Nombre de bits / symbole
Fe = 10000 ; % Frequence d'echantillonnage
Rs = 1000 ; % Debit symbole
Ns = Fe/Rs ; % Nombre d'echantillons / periode symbole
alphaSRRCF = 0.35 ; % Coefficient roll off
long_filtre_Ts = 6 ; % Longueur du filtre en SRRCF en nombre de periode Ts
h=rcosdesign(alphaSRRCF,long_filtre_Ts,Ns);
nfft=2^nextpow2(Nbits/n*Ns) ; % Nombre de points de la FFT
axe_f=linspace(-0.5,0.5,nfft) ; % Axe des frequences
fp=2000;
fpN=fp/Fe;

% Emission
bits=randi([0 1],1,Nbits);  % Generation des bits
symbI=1-2*bits(1:2:Nbits); % Codage symboles sur la voie I
symbQ=1-2*bits(2:2:Nbits); % Codage symboles sur la voie Q
peigne=zeros(1,Nbits/n*Ns);  
peigne(1:Ns:Nbits/n*Ns)=symbI+1i*symbQ; % Symboles du peigne de Dirac
xe=filter(h,1,peigne) ; % Generation du signal passe-bas
x=real(exp(1i*2*pi*fpN*[1:Nbits/n*Ns]).*xe) ; % Generation du signal passe-bande

figure
semilogy(axe_f,fftshift(abs(fft(xe,nfft)).^2))
hold on
semilogy(axe_f,fftshift(abs(fft(x,nfft)).^2),'r')
xlabel('frequence')
title('Densite spectrale de puissance')
legend('S_xe(f)','S_x(f)')
grid

