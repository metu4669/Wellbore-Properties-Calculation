clear
SGm = 0.71 % Specific Gravity of Mixture
N2m = 0.02 % Mole Fraction of N2
CO2m = 0.01 % Mole Fraction of CO2
H2Sm = 0.07 % Mole Fraction of H2S

Pwh = 110 % Well Head Pressure(psi)
Twh = 100 % Well Head Temp(F)
Pbh = 1200 % Bottomhole Pressure(psi)
Tbh = 170 % Bottomhole Temperature(F)
 %_______________________________________%
 
 %{
 To calculate Ppc and Tpc, we will use Sutton Corralation.
 %}
  
 SPpc = 756.8-131*SGm-3.6*SGm*SGm
 STpc = 169.2+349.5*SGm-74*SGm*SGm
 
 
 %{
We will correct our result with Wichert-Aziz Correlation and this is the
final corrected results
 %}
 A= H2Sm + CO2m
 B= H2Sm
 Eps = 120*(power(A,0.9)-power(A,1.6))+15*(power(B,0.5)-power(B,4))
 Tpc =STpc-Eps %Final Corrected Tpc
 Ppc = SPpc*Tpc/(STpc+H2Sm*(1-H2Sm)*Eps) %Final Corrected Ppc
 
%  To calculate ZFactor, I used Redlich Kwong Z Factor formula.
k = 150 %Division #
Ma = SGm * 29 % Apparent Molecular Weight
R = 10.7316

Temperature = transpose(linspace(Twh,Tbh,k))+460
Pressure = transpose(linspace(Pwh,Pbh,k))

TemperatureR = Temperature/Tpc
PressureR = Pressure/Ppc


TmPr = zeros(k,2) % Temp and press relatively depth
TPr = zeros(k,2) %Relative to depth, 1 st col is temp and 2nd col is pressure.
 

 Density = zeros(k,1) %Density Matrix
 ZFactor = zeros(k,1) %Z Factor Matrix

 
 for i = 1:k               
%-----------------------------------------------------------------------%
    L = 0.4278*PressureR(i,1)/TemperatureR(i,1)^2.5;
	M = 0.0867*PressureR(i,1)/TemperatureR(i,1);
	ZfactorDeployed = [1 -1 (L-M+M^2) -L*M];
	Zx1 = roots(ZfactorDeployed);
	Zx1 = Zx1(imag(Zx1)==0); % Save only the real roots
	Zx1 = sort(Zx1);
	Z = Zx1(1);
    ZFactor(i,1) = Z(1)
 end
 
 % We will use Lee-Gonzales and Eakin Correlation for calculation of visc

 LGECX = zeros(k,1)
 LGECY = zeros(k,1)
 LGECK = zeros(k,1)
 Viscosity = zeros(k,1)
 
 for i = 1:k
    
     
     Density(i,1) = Pressure(i,1)*Ma /(ZFactor(i,1)*R*Temperature(i,1))
     
     LGECX(i,1) = 3.5+(986/Temperature(i,1))+0.01*Ma
     
     LGECY(i,1) = 2.4 - 0.2*LGECX(i,1)
     
     LGECK(i,1) = (9.4+0.02*Ma)*(Temperature(i,1)^1.5)*(10^(-4))/(209+19*Ma+Temperature(i,1))
     Viscosity(i,1) = (LGECK(i,1))* exp(LGECX(i,1)*((Density(i,1)/62.4)^LGECY(i,1)))
 
 end
 
 % Pseudo Pressure
 l = 10
 
 P = zeros(l,1)
 u = zeros(l,1)
 z = zeros(l,1)
 Den = zeros(l,1)
 Del = zeros(l,1)

 Pseudo = zeros(k,1)
 for i = 1:k
    for t = 1:l+1     
        Press = Pressure(i,1)
        P(t,1) = (t-1)*Press/l
%-----------------------------------------------
    L2 = 0.4278*P(t,1)/Ppc/(((Temperature(i,1))/Tpc)^2.5);
	M2 = 0.0867*P(t,1)/Ppc/((Temperature(i,1))/Tpc);
	ZfactorDeployed2 = [1 -1 (L2-M2+M2^2) -L2*M2];
	Zx2 = roots(ZfactorDeployed2);
	Zx2 = Zx2(imag(Zx2)==0); % Save only the real roots
	Zx2 = sort(Zx2);
	Z2 = Zx2(1);
    z(t,1) = Z2(1)
%--------------------------------------------------
        Den(t,1) = P(t,1)*Ma/(z(t,1)*R*Temperature(i,1))
        LX = 3.5+(986/Temperature(i,1))+0.01*Ma
        LY = 2.4 - 0.2*LX
        LK = (9.4+0.02*Ma)*(Temperature(i,1)^1.5)*(10^(-4))/(209+19*Ma+Temperature(i,1))
        u(t,1) = LK* exp(LX*((Den(t,1)/62.4)^LY))
        Del(t,1) = 2 * P(t,1)/(u(t,1)*z(t,1)) 
        Del(isnan(Del)) = 0
       S = sum(Del)
       Pseudo(i,1) = 2*S(1,1) 
    end
    
    Pseudo(i,1) = P(t,1)*(Pseudo(i,1)-Del(1,1)-Del(l,1))/(2*l)
    
 end
 
 
 %Comparing Part
%-----------------------------------------------------------------------
%I will use PVTProp ZFactor Calculation Tool

%We will choose some random Temperature values from the well and compare
%them with our values
sample = 110
Random = randi([1 k],sample,1)
Random =sort(Random)
 ZFRandom = zeros(sample,2) % 1st Column Ours, 2nd Column PVT
 
 for i = 1:sample
    ZFRandom(i,1) = ZFactor(Random(i,1),1) 
    ZFRandom(i,2) = ZFactorPVTProps(TemperatureR(Random(i,1),1),PressureR(Random(i,1),1))
 end
 
 figure
 plot(Pressure(Random(1:sample),1),ZFRandom(1:sample,1),Pressure(Random(1:sample),1),ZFRandom(1:sample,2))
 ylabel('ZFactor')
 xlabel('Pressure')
 title('ZFactor vs Pressure')
%-----------------------------------------------------------------------
 %Drawing
 x = 1:k
 y = Density(1:k)
 z = Viscosity(1:k)
 c = Pseudo(1:k)

 figure
 scatter(x,y) 
 title('Density(lb/ft^3) vs Segment Number')
 xlabel('Segment Number')
 ylabel('Density(lb/ft^3)')
 
 figure
 scatter(x,z)
 title('Viscosity(cp) vs Segment Number')
 xlabel('Segment Number')
 ylabel('Viscosity(cp)')
 
 figure
 scatter(x,c)
 title('Pseudo Pressure(psi^2/cp) vs Segment Number')
 xlabel('Segment Number')
 ylabel('Pseudo Pressure(psi^2/cp)')
 
 figure
 scatter(x,ZFactor(1:k,1))
 title('ZFactor vs Segment Number')
 xlabel('Segment Number')
 ylabel('ZFactor')
 clc
                
 
