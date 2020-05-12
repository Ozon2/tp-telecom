%% Comparaison de modulatins sur fréquence proteuse
clear all
close all

%% Initialisation
Eb_sur_N0_dB = [1:6];                       % vecteur Eb/N0 en dB
Eb_sur_N0 = 10.^(Eb_sur_N0_dB./10);         % vecteur Eb/N0 linéaire

Nerrlimite = 100;                           % nombre d'erreur limite pour une precision de 10%

Nbits = 3000;                               % Nombre de bits emis

Fe=144000;                                   % frequence d'echantillonnage
Rb=48000;                                   % debit binaire
retard_Ts=6;                                % Longueur du filtre en SRRCF en nombre de periode Ts
alphaSRRCF=0.5;                             % Coefficient roll off
Nbrep = 20;


%% Chaine 4-ASK

% Initialisation
TEB1s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
M = 4;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

for i=1:Nbrep
    % Emission
    bits=randi([0 1],1,Nbits);      % Generation des bits
    symboles = A_PSKmod(bits, Nbits);
    peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
    peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
    xe1=filter(h,1,peigne_allonge); % Generation du signal passe-bas
    
    % Canal
    signal_canal = xe1; % Signal du canal
    
    % Reception
    z=filter(hr,1,signal_canal);                               % Signal canal filtré
    % Echantilonage
    ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
    
    % Decision
    symbole_estimes = 3*(real(ze)>2)-3*(real(ze)<-2)+(real(ze)>0 & real(ze)<2)-(real(ze)<0 & real(ze)>-2);
    
    % Demapping
    bits_estimes = A_PSKdemod(symbole_estimes, Nbits);
    
    TEB1 = sum(bits ~= bits_estimes); % Comparaison des bits
    if (i==1)
        fprintf("Le TEB sans bruit vaut : %d \n", TEB1);
    end
    
    % Implantation avec bruit
    Pre = mean(abs(xe1).^2);                       % Calcul de la puissance du signal recu
    Sigma2 = Pre*Ns./(2*log2(M)*Eb_sur_N0);       % Calcul de la variance du bruit
    
    Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
    nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
    for k=1:length(Eb_sur_N0_dB)
        
        while (Nerr(k) < Nerrlimite)
            
            
            
            % Canal
            signal_canal = xe1 + sqrt(Sigma2(k))*(randn(1,length(xe1)) + 1i*randn(1,length(xe1))); % Signal du canal
            
            % Reception
            z=filter(hr,1,signal_canal);                               % Signal canal filtré
            % Echantilonage
            ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
            
            % Decision
            symbole_estimes = 3*(real(ze)>2)-3*(real(ze)<-2)+(real(ze)>0 & real(ze)<2)-(real(ze)<0 & real(ze)>-2);
            
            % Demapping
            bits_estimes = A_PSKdemod(symbole_estimes, Nbits);
            
            
            NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
            Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
            
            nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
        end
        TEB1s(k) = TEB1s(k) + Nerr(k)/(nbEssais(k)*Nbits);
    end
end
TEB1s = TEB1s/Nbrep;

TES_theo = 2*(1-1/M)*qfunc(sqrt(6*log2(M)/(M^2-1)*Eb_sur_N0));
TEB_theo = TES_theo/log2(M);
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB1s,'blue *'); hold on
semilogy(Eb_sur_N0_dB,TEB_theo,'magenta');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé 4-ASK","TEB théorique 4-ASK");

figure;
plot(real(ze) + eps*1i,'b+'); hold on
plot(symboles + eps*1i,'r*', 'LineWidth',3);
plot([-0.5,0.5], [0,0], 'b-');
plot([0,0],[-0.5,0.5],'g-');
plot([2,2],[-0.5,0.5],'g-');
plot([-2,-2],[-0.5,0.5],'g-');
title("Constellation 4-ASK");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles sortie échantillonneur","symboles sortie mapping");

%% Chaine QPSK

% Initialisation
TEB2s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
bits_estimes = zeros(1,Nbits);              % Initialisation de bits_estime
M = 4;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

for i=1:Nbrep
    % Emission
    bits=randi([0 1],1,Nbits);      % Generation des bits
    bits_reshape = reshape(bits, 2, Nbits/2);
    bits_dec = bi2de(bits_reshape');
    symboles = pskmod(bits_dec',M,pi/M,'gray');
    peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
    peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
    xe2=filter(h,1,peigne_allonge); % Generation du signal passe-bas
    
    % Canal
    signal_canal = xe2; % Signal du canal
    
    % Reception
    z=filter(hr,1,signal_canal);                               % Signal canal filtré
    % Echantilonage
    ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
    
    % Demapping
    bits_estimes_dec = pskdemod(ze,M,pi/M,'gray');
    bits_estimes_reshape = de2bi(bits_estimes_dec)';
    bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
    
    TEB2 = sum(bits ~= bits_estimes); % Comparaison des bits
    if (i==1)
        fprintf("Le TEB sans bruit vaut : %d \n", TEB2);
    end
    
    % Implantation avec bruit
    Pre = mean(abs(xe2).^2);                       % Calcul de la puissance du signal recu
    Sigma2 = Pre*Ns./(2*log2(M)*Eb_sur_N0);       % Calcul de la variance du bruit
    
    Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
    nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
    for k=1:length(Eb_sur_N0_dB)
        
        while (Nerr(k) < Nerrlimite)
            
            
            
            % Canal
            signal_canal = xe2 + sqrt(Sigma2(k))*(randn(1,length(xe2)) + 1i*randn(1,length(xe2))); % Signal du canal
            
            % Reception
            z=filter(hr,1,signal_canal);                               % Signal canal filtré
            % Echantilonage
            ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
            
            % Demapping
            bits_estimes_dec = pskdemod(ze,M,pi/M,'gray');
            bits_estimes_reshape = de2bi(bits_estimes_dec)';
            bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
            
            
            NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
            Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
            
            nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
        end
        TEB2s(k) = TEB2s(k) + Nerr(k)/(nbEssais(k)*Nbits);
    end
end
TEB2s = TEB2s/Nbrep;

TEB_theo = qfunc(sqrt(2*Eb_sur_N0));
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB2s,'black *'); hold on
semilogy(Eb_sur_N0_dB,TEB_theo,'magenta');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé QPSK","TEB théorique QPSK");

figure;
plot(ze,'b+'); hold on
plot(symboles,'r*','LineWidth',3);
plot([-2,2], [0,0], 'g-');
plot([0,0],[-2,2],'g-');
title("Constellation QPSK");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles sortie échantillonneur","symboles sortie mapping");

%% chaine 8-PSK

% Initialisation
TEB3s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
M = 8;                                      % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

for i=1:Nbrep
    % Emission
    bits=randi([0 1],1,Nbits);      % Generation des bits
    bits_reshape = reshape(bits, 3, Nbits/3);
    bits_dec = bi2de(bits_reshape');
    symboles = pskmod(bits_dec',M,pi/M,'gray');
    peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
    peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
    xe3=filter(h,1,peigne_allonge); % Generation du signal passe-bas

    % Canal
    signal_canal = xe3; % Signal du canal

    % Reception
    z=filter(hr,1,signal_canal);                               % Signal canal filtré
    % Echantilonage
    ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);

    % Demapping
    bits_estimes_dec = pskdemod(ze,M,pi/M,'gray');
    bits_estimes_reshape = de2bi(bits_estimes_dec)';
    bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
        
    TEB3 = sum(bits ~= bits_estimes); % Comparaison des bits
    if (i==1)
        fprintf("Le TEB sans bruit vaut : %d \n", TEB3);
    end
    
    % Implantation avec bruit
    Pre = mean(abs(xe3).^2);                       % Calcul de la puissance du signal recu
    Sigma2 = Pre*Ns./(2*log2(M)*Eb_sur_N0);       % Calcul de la variance du bruit

    Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
    nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
    for k=1:length(Eb_sur_N0_dB)
    
        while (Nerr(k) < Nerrlimite)
        
        

            % Canal
            signal_canal = xe3 + sqrt(Sigma2(k))*(randn(1,length(xe3)) + 1i*randn(1,length(xe3))); % Signal du canal

            % Reception
            z=filter(hr,1,signal_canal);                               % Signal canal filtré
            % Echantilonage
            ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
        
            % Demapping
            bits_estimes_dec = pskdemod(ze,M,pi/M,'gray');
            bits_estimes_reshape = de2bi(bits_estimes_dec)';
            bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
        
        
            NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
            Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
        
            nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
        end
        TEB3s(k) = TEB3s(k) + Nerr(k)/(nbEssais(k)*Nbits);
    end
end
TEB3s = TEB3s/Nbrep;

TES_theo = 2*qfunc(sqrt(2*log2(M)*Eb_sur_N0)*sin(pi/M));
TEB_theo = TES_theo/log2(M);
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB3s,'green *'); hold on
semilogy(Eb_sur_N0_dB,TEB_theo,'magenta');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé 8-PSK","TEB théorique 8-PSK");

figure;
plot(ze,'b+'); hold on
plot(symboles,'r*','LineWidth',3);
plot([-2,2], [0,0], 'g-');
plot([0,0],[-2,2],'g-');
plot([-2,2],[-2,2],'g-');
plot([-2,2],[2,-2],'g-');
title("Constellation 8-PSK");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles sortie échantillonneur","symboles sortie mapping");

%% chaine 16-QAM

% Initialisation
TEB4s = zeros(1,length(Eb_sur_N0));         % Initialisation du TEB
M = 16;                                     % Taille de la constellation
Rs = Rb/log2(M);                            % debit symbole
Ns=Fe/Rs;                                   % Nombre d'echantillons par periode porteuse
h=rcosdesign(alphaSRRCF,retard_Ts,Ns);      % Filtre d'emission
hr = fliplr(h);                             % Filtre de réception adapte
t0 = retard_Ts*Ns+1;                        % Instant de prise de décision pour d0 emis a t=0

for i=1:Nbrep
    % Emission
    bits=randi([0 1],1,Nbits);      % Generation des bits
    bits_reshape = reshape(bits, 4, Nbits/4);
    bits_dec = bi2de(bits_reshape');
    symboles = qammod(bits_dec',M,'gray');
    peigne= kron(symboles, [1, zeros(1,Ns-1)]); % Symboles du peigne de Dirac
    peigne_allonge=[peigne zeros(1,retard_Ts*Ns)];
    xe4=filter(h,1,peigne_allonge); % Generation du signal passe-bas
    
    % Canal
    signal_canal = xe4; % Signal du canal
    
    % Reception
    z=filter(hr,1,signal_canal);                               % Signal canal filtré
    % Echantilonage
    ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
    
    % Demapping
    bits_estimes_dec = qamdemod(ze,M,'gray');
    bits_estimes_reshape = de2bi(bits_estimes_dec)';
    bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
    
    TEB4 = sum(bits ~= bits_estimes); % Comparaison des bits
    if (i == 1)
        fprintf("Le TEB sans bruit vaut : %d \n", TEB4);
    end
    
    % Implantation avec bruit
    Pre = mean(abs(xe4).^2);                       % Calcul de la puissance du signal recu
    Sigma2 = Pre*Ns./(2*log2(M)*Eb_sur_N0);       % Calcul de la variance du bruit
    
    Nerr=zeros(1,length(Eb_sur_N0_dB));         % Initialisation du nombre d'erreur
    nbEssais=zeros(1,length(Sigma2));           % Initialisation du nombre d'essais
    for k=1:length(Eb_sur_N0_dB)
        
        while (Nerr(k) < Nerrlimite)
            
            % Canal
            signal_canal = xe4 + sqrt(Sigma2(k))*(randn(1,length(xe4)) + 1i*randn(1,length(xe4))); % Signal du canal
            
            % Reception
            z=filter(hr,1,signal_canal);                               % Signal canal filtré
            % Echantilonage
            ze = z(t0:Ns:(Nbits/log2(M)+retard_Ts)*Ns);
            
            % Demapping
            bits_estimes_dec = qamdemod(ze,M,'gray');
            bits_estimes_reshape = de2bi(bits_estimes_dec)';
            bits_estimes = reshape(bits_estimes_reshape, 1, Nbits);
            
            
            NerrActuel = sum(bits ~= bits_estimes); % Comparaison des bits
            Nerr(k) = Nerr(k) + NerrActuel; % Calcul du nombre d'erreur
            
            nbEssais(k) = nbEssais(k) + 1; % Calcul du nombre d'essais
        end
        TEB4s(k) = TEB4s(k) + Nerr(k)/(nbEssais(k)*Nbits);
    end
end
TEB4s = TEB4s/Nbrep;

TES_theo = 4*(1-1/sqrt(M))*qfunc(sqrt(3*log2(M)/(M-1)*Eb_sur_N0));
TEB_theo = TES_theo/log2(M);
Pb = TEB_theo; % Puissance du bruit
variance_simu=Pb.*(1-Pb)./(nbEssais.*Nbits); % Variance de la simu

figure;
semilogy(Eb_sur_N0_dB,TEB4s,'red *'); hold on
semilogy(Eb_sur_N0_dB,TEB_theo,'magenta');
semilogy(Eb_sur_N0_dB,Pb+sqrt(variance_simu),'c')
semilogy(Eb_sur_N0_dB,Pb-sqrt(variance_simu),'c')
title("TEB en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB simulé 16-QAM","TEB théorique 16-QAM");

figure;
plot(ze,'b+'); hold on
plot(symboles,'r*','LineWidth',3);
plot([-5,5], [0,0], 'g-');
plot([-5,5], [2,2], 'g-');
plot([-5,5], [-2,-2], 'g-');
plot([0,0],[-5,5],'g-');
plot([2,2],[-5,5],'g-');
plot([-2,-2],[-5,5],'g-');
title("Constellation 16-QAM");
xlabel("Partie réelle");
ylabel("Partie imaginaire");
legend("symboles sortie échantillonneur","symboles sortie mapping");

%% Figures

figure;
semilogy(Eb_sur_N0_dB,TEB1s,'blue-+'); hold on
semilogy(Eb_sur_N0_dB,TEB2s,'black-o'); hold on
semilogy(Eb_sur_N0_dB,TEB3s,'green-*'); hold on
semilogy(Eb_sur_N0_dB,TEB4s,'red-d'); hold on
title("TEB des différentes chaines en fonction de (Eb/N0) (dB)");
xlabel("(Eb/N0) (dB)");
ylabel("TEB");
legend("TEB chaine 4-ASK","TEB chaine QPSK","TEB chaine 8-PSK","TEB chaine 16-QAM");

nfft1=2^nextpow2(Nbits/log2(4)*Ns);      % Nombre de points de la FFT
axe_f1=linspace(-0.5,0.5,nfft1) ;          % Axe des frequences
nfft2=2^nextpow2(Nbits/log2(4)*Ns);      % Nombre de points de la FFT
axe_f2=linspace(-0.5,0.5,nfft2) ;          % Axe des frequences
nfft3=2^nextpow2(Nbits/log2(8)*Ns);      % Nombre de points de la FFT
axe_f3=linspace(-0.5,0.5,nfft3) ;          % Axe des frequences
nfft4=2^nextpow2(Nbits/log2(16)*Ns);     % Nombre de points de la FFT
axe_f4=linspace(-0.5,0.5,nfft4) ;          % Axe des frequences

figure
subplot(2,2,1); semilogy(axe_f1,fftshift(abs(fft(xe1,nfft1)).^2), 'blue'); hold on
title("DSP chaine 4-ASK");
xlabel('frequence normalisee')
grid
subplot(2,2,2); semilogy(axe_f2,fftshift(abs(fft(xe2,nfft2)).^2), 'black');
title("DSP chaine QPSK");
xlabel('frequence normalisee')
grid
subplot(2,2,3); semilogy(axe_f3,fftshift(abs(fft(xe3,nfft3)).^2), 'green');
title("DSP chaine 8-PSK");
xlabel('frequence normalisee')
grid
subplot(2,2,4); semilogy(axe_f4,fftshift(abs(fft(xe4,nfft4)).^2), 'red');
title("DSP chaine 16-QAM");
xlabel('frequence normalisee')
grid

