% QPSK en bande transposee

clear all
close all
clc

% Initialisation
Eb_N0_dB=[1:6]; % vecteur Eb/N0 en dB
Eb_N0_lin=10.^(Eb_N0_dB./10);

Nerr=zeros(1,length(Eb_N0_dB))
Nerrlimite = 100; % precision de 10%

Nbits=10000;
M=4;
offsetQPSK=pi/4;

Fe=10000; % frequence d'echantillonnage
Rs=1000; % debit symbole
Ns=Fe/Rs; % Nombre d'echantillons par periode porteuse
retard_Ts=6;
alphaSRRCF=0.35;
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);

fp=2000; % frequence porteuse
fpN=fp/Fe    % frequence porteuse normalisee


% for k=1:length(Eb_N0_dB)
    
  %  while (Nerr(k)<Nerrlimite)

        % Emission
        bits=randi([0 1],1,Nbits)
        test=reshape(bits,2,Nbits/2)
        bits_dec=bi2de(test')
        symboles=pskmod(bits_dec,M,offsetQPSK,'gray')/(sqrt(2)/2)
        peigne=zeros(1,Nbits/log2(M)*Ns);
        peigne(1:Ns:Nbits/log2(M)*Ns)=symboles;
        signal_passe_bas=filter(h,1,peigne);
        signal_passe_bande=real(signal_passe_bas).*cos(2*pi*fpN.*[1:Nbits/log2(M)*Ns]) ...
            - imag(signal_passe_bas).*sin(2*pi*fpN.*[1:Nbits/log2(M)*Ns]);

        % Canal


        % Reception
        
   % end
% end

nfft=2^nextpow2(Nbits/log2(M)*Ns);
semilogy(linspace(-0.5,0.5,nfft),fftshift(abs(fft(signal_passe_bas,nfft)).^2),'r')
hold on
semilogy(linspace(-0.5,0.5,nfft),fftshift(abs(fft(signal_passe_bande,nfft)).^2))
xlabel('frequence normalisee')
grid
legend('PSD signal passe-bas','PSD signal passe-bande')
