function beta=calcBeta1(L_c,rho,lambda,mu_a,mu_sp,BFi,n)

%%

%% Calculation of speckle contrast beta. Code modified by Antonio Ortega (BU) 
%% based on code provided by Xiaojun Cheng (BU)

%             L_c=.05;% cm  An ideal sech^2 0.5 ns FWHM pulse has a linewidth of 0.315 / 0.5 ns = 0.63 GHz
%             %The coherence length of such a transform-limited pulse is c / pi / 0.63  GHz = 15 cm
%
%             rho = .9; % cm distance between source and detector
%             lambda = 830e-7; %cm wavelength
%
%             mu_a = 0.1; % 1/cm
%             mu_sp = 10; % 1/cm
%
%             BFi = 400e-12;   % cm^2/s
%
%             n = 1.4;

v=299792458e2/n; % cm/s

%% TPSF
D=1/(3*(mu_a+mu_sp));  %difussion constant for photon flow
alb=3*mu_sp.*D; %albedo

mu_eff=sqrt(mu_a./D);

sigma=0.1;
mu = 5*sigma;

t_irf=-mu:mu/100:mu*4;

t=t_irf(t_irf>0)*1e-9; % s  times of flights in seconds

Slt = t * v;  %distance of photon flight?

z0=1/mu_sp;
Reff=.493;   %the function of refraction index I guess
A=(1+Reff)/(1-Reff);
zb=2*A*D;  %extrapolated boundary location
zminus=sqrt(z0^2+rho^2);   %virtual source location -
zplus=sqrt((z0+2*zb)^2+rho^2);  %virtual source location +
phi=v./(4*pi*D*Slt).^(3/2).*(exp(-zminus.^2./(4*D*Slt))-exp(-zplus.^2./(4*D*Slt))).*exp(-mu_a*Slt); %TOF distribution
T = getRr(mu_a, mu_sp, rho, n);
% Plt = phi / max(phi);   %normalized to amplitude=1, but not area=1
%             t = t * 1e9; %ns   times of flight in nanoseconds
Plt=phi;
%%
% get the path length distribution for the given sample time delay

% flip the IRF around t_sample


Plt = Plt / sum(Plt);

k0 = 2*pi/lambda; % Free-space wavenumber
k = k0*n; % Wavenumber in tissue

% tau = [0, logspace(-7,-2,300)]'; %time delays for autocorrelation (?)

tau=0;

%% Optimized code to obtain beta using Bellini

[Plt1,Plt2]=meshgrid(Plt,Plt);
Pltp=Plt1.*Plt2;

[Slt1,Slt2]=meshgrid(Slt,Slt);

%                        BFi=1/(6*tau_c*mu_sp^2*k0^2*rho^2);

CorrMatrix=zeros(size(Plt1,1),size(Plt1,2),length(tau));
%             CorrMatrix0=CorrMatrix;
for ii=1:length(tau)
    %                 CorrMatrix0(:,:,ii)=Pltp.*exp( -2*mu_sp*BFi*k^2*Slt1*tau(ii)).*exp( -2*mu_sp*BFi*k^2*Slt2*tau(ii));
    CorrMatrix(:,:,ii)=Pltp.* exp( -2*mu_sp*BFi*k^2*Slt1*tau(ii)).*exp( -2*mu_sp*BFi*k^2*Slt2*tau(ii)).*exp(-2*(Slt1-Slt2).^2/L_c.^2);
end
%             CorrMatrix0=CorrMatrix0/sum(Plt).^2;
CorrMatrix = CorrMatrix/sum(Plt).^2;

%             g20=1+squeeze(sum(sum(CorrMatrix0)));
g2=1+squeeze(sum(sum(CorrMatrix)));
%   figure(99)
%             semilogx(tau,g2)
%% Vectorized previous optimization but no performance gain

%             [Plt1,Plt2]=meshgrid(Plt,Plt);
%             Pltp=Plt1.*Plt2;
%             [Slt1,Slt2,Tau]=meshgrid(Slt,Slt,tau);
%
%
%             CorrMatrix0=Pltp.*exp( -2*mu_sp*BFi*k^2*Slt1.*Tau).*exp( -2*mu_sp*BFi*k^2*Slt2.*Tau);
%             CorrMatrix= Pltp.*exp( -2*mu_sp*BFi*k^2*Slt1.*Tau).*exp( -2*mu_sp*BFi*k^2*Slt2.*Tau).*exp(-2*(Slt1-Slt2).^2/L_coherence.^2);
%
%             CorrMatrix0=CorrMatrix0/sum(Plt).^2;
%             CorrMatrix = CorrMatrix/sum(Plt).^2;
%
%             g20=1+squeeze(sum(sum(CorrMatrix0),2));
%             g2=1+squeeze(sum(sum(CorrMatrix),2));

%%

beta=g2(1)-1;

