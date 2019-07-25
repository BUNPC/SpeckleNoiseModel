% Example of calculation of noise for a breast phantom, OPT101, using VHG

clear

rho=3e-2;  % source-detector separation

c=299792458;   %m/s
qe=1.60217662e-19; %C
h=6.626e-34;  %   Planck's constant Js


%Source 1. VHG
SM_VHG.wavelength=785e-9;
SM_VHG.coherencelength=10;  %m
SM_VHG.Power=50e-3;
SM_VHG.Power=logspace(-11,0,1000);

%medium. Breast phantom
breastphantom.BFi=40e-15;
breastphantom.alpha=1;
breastphantom.n=1.4;
breastphantom.musp=900;  %per meter; if scalar, assumes spectrally flat
breastphantom.mua=6.5;  %per meter



%OPT101, photodiode only
OPT101.NA=0.22;
OPT101.radius=100e-6;  %sensing area assuming fully illuminated (collected with a fiber of this core diameter)
OPT101.darkcurrent=2.5e-12; %should be specified for internal gain=1. Mean value.
OPT101.R0=0.57;  %A/W, if scalar, assumes responsivity is spectrally flat
OPT101.capacitance=1.2e-9;  %produces additional noise
OPT101.label='Photodiode';


%dummy electronics
electronics.integrationtime=0.5;  %should be 1/fs
electronics.gain=100; %electronic gain
electronics.bandwidth=1; %should be <= fs/2, considered single sided; code enforces this
% dummyel.conversionfactor=7.5e12*40e-3; 
electronics.conversionfactor=6e12; %this is the conversion factor from photocurrent to DL for Gain=1
electronics.Rsh=1e9;  %load for photodiode
electronics.T=25;  %load for photodiode

photodetector=OPT101;
medium=breastphantom;
source=SM_VHG;

[SNRdB,PD,Noise,NEP,Noisecomponents,Noisefraction,psi]=noisecalcfunc(rho,source,medium,photodetector,electronics);


%%
figure(1)
semilogx(PD,SNRdB)
xlabel('Power incident in detector [W]')
ylabel('SNR [dB]')
ylim([0 110])
grid on
xlim([NEP,max(PD)])




figure(2)
loglog(PD,Noise/psi)
xlim([min(PD),max(PD)])
xlim([NEP*sqrt(1),max(PD)])
xlabel('Power incident in detector [W]')
ylabel('Equivalent optical power noise [W]')
grid on

