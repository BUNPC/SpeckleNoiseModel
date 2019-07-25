function [SNRdB,PD,Noise,NEP,Noisecomponents,Noisefraction,psi,beta]=noisecalcfunc(rho,source,medium,photodetector,electronics)
%Calculates the noise for NIRS as a function of soruce detector separation,
%for a specific light source, medium, photodetector and electronics.
%Returns the signal to noise ratio (dB), the signal (PD), Noise, how the
%noise is distributed in three categories (electronic, shot and speckle),
%the NEP, the fraction of variance contributed by each noise component, the
%"sensitivity of the system" at a given gain psy, and the speckle contrast
%constan beta

%Coded by Antonio Ortega, B.O.A.S. lab Boston University, auxiliary
%functions from Xiaojun Cheng (BU) and Steven Jacques (https://omlc.org/news/apr08/skinspectra/index.html)

%assumes CV of laser source is 0 (no source drift(

%all units are SI. So that means that extinction coefficients are in 1/m
%rho is source detector separation 

%source is a structure that defines the characteristics of the source such
%as wavelength, coherence length, and power. The calculations assume a
%monochromatic source despite the coherence length

%medium is a structure containing the optical characteristics of the medium

%electronics is a structure containing the characteristics of the
%electronics post photodetector; they are assumed to include an electronic
%gain stage after photodetection and potentially an ADC that converts the
%photocurrent to DL

%photodetector is a structure containing the characteristics of the
%photodetector: photodetector.radius is used to calculate the area of
%sensing and thus the amount of power calculated at a distance rho;
%noisefactor is the excess noise factor; for a PIN diode it should be 1;
%photodetector.gain is the internal gain; should be 1 unless it's an APD.
%photodetector.darkcurrent is the mean dark current specified in the
%datasheet; do note this is a function of the APD gain


%%


%Physical constants
c=299792458;   %m/s
qe=1.60217662e-19; %C
h=6.626e-34;  %   Planck's constant Js
kB=1.38064852e-23; %Boltzmann constant J/K


T=273+electronics.T;

a=photodetector.radius;
P0=source.Power;

lambda=source.wavelength;
L_c=source.coherencelength;


G=electronics.gain;
tau=electronics.integrationtime;
BW=electronics.bandwidth;
Rsh=electronics.Rsh;


if BW>(1/2/tau)
    BW=1/2/tau;
end

Z=electronics.conversionfactor;

if isfield(photodetector,'noisefactor')
    F=photodetector.noisefactor;
else
    F=1;
end

if isfield(photodetector,'gain')
    gamma=photodetector.gain;
else
    gamma=1;
end
    

idark=photodetector.darkcurrent;

alpha=medium.alpha;
BFi=medium.BFi;


n=medium.n;
mus=medium.musp;
mua=medium.mua;
Ep=h*c/lambda; %photon energy

if isfield(photodetector,'R0')
    R0=photodetector.R0;  %detector responsivity (internal gain=1)
else
    %if R0 is not present, calculate from quantum efficiency, assuming it
    %is there
    R0=photodetector.quanteff*qe/Ep;
end


NA=photodetector.NA;

if isfield(photodetector,'M')
    M=photodetector.M;
else
    M=(2*pi*a*NA/lambda).^2/2; %This approximation is actually not good unless M is big
    if M<3
        M=1;
    end
end

Nidark=sqrt((qe/tau*(gamma^2*F))*idark+2*kB*T/tau/Rsh);   %RMS of dark count; it's supposed to be Poisson distributed, A

k0=2*pi/lambda;

T = getRr(mua, mus, rho, n); %photons/m^2 using diffuse reflectance model
PD=P0*T*pi.*a.^2;
R=gamma*R0;  % Photodetector responsivity A/W
psi=Z*G*R;  %sensitivity of system a.u./W

tau_c=1./(6*BFi*mus.^2*k0.^2*alpha*rho.^2);

%% calculation of beta
beta=calcBeta1(L_c*1e2,rho*1e2,lambda*1e2,mua/1e2,mus/1e2,BFi*1e4,n);
beta=beta/2;  %because of polarization. Maybe.


%% calculation of noise
Nspeck=psi*sqrt(beta*tau_c/tau./M).*PD*sqrt(BW*tau*2);
%Nshot=sqrt(psi*Z*G*qe*gamma*F)*sqrt(PD)*sqrt(BW*tau*2);
Nshot=psi*sqrt(qe*F/tau/R0)*sqrt(PD)*sqrt(BW*tau*2);
Nelec=psi*Nidark/R*sqrt(BW*tau*2)*ones(size(Nshot));      
Noise=sqrt(Nspeck.^2+Nshot.^2+Nelec.^2); %model is not taking into account amplifier noise
Noisecomponents=[Nelec',Nshot',Nspeck'];
Noisefraction=Noisecomponents.^2./Noise'.^2;

%% Calculation of signal and SNR
S=psi*PD;
SNRdB=20*log10(S./Noise);
[~,indi]=near(SNRdB,0);
NEP=PD(indi)./sqrt(BW);

A_shotwatts=sqrt(qe*F/tau/R0)*sqrt(BW*tau*2);
N_elecwatts=Nelec*G./psi;



