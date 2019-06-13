
function [BIR4s1, BIR4s2decHinder, BIR4s2obs, BIR4s3Hinder, BIR4180] = createIRRUPT(zeta,kappa,Tp,deltaomegamax,res,delta)
%this function generates a complex description of a the pulses for an IRRUPT sequence 
%BIR4s1 is the first segment played in I and S  
%BIR4s2decHinder is the second segment played in S 
%BIR4s2obs is the second segment played in I 
%BIR4s3Hinder is the third segment played in I and S
%BIR4180 is the adiabatic inversion played on both I and S
%delta is the phasejump that controls the amount of transfered
%res is the spacing in seconds between two timepoints of the output
%Tp is the duration in seconds for one BIR4 pulse
%If delta = 90 this generates pulsesfor an BINEPT sequence (Merkle 1992) 
%zeta, kappa are parameters as described in Garwood 1991
%recommended zeta = 10, kappa = atan(20)
%deltaomegamax is the sweepwidth, recommended in Garwood 1991:  2 pi 45/Tp
deltaphidec = pi + (delta*pi/180)/2 ;
deltaphiobs = pi ;
deltaphi180 = pi + pi/2 ;

phimax = -deltaomegamax.*Tp.*log(cos(kappa))/(kappa.*tan(kappa)) ;
tpcal = 0:res:(Tp/4-res);           % in s.
dT = tpcal(2)-tpcal(1);        
revAHP = tanh(zeta.*(1-4.*tpcal/Tp)).*exp(1i.*(phimax/4 -(deltaomegamax.*Tp/(4.*kappa.*tan(kappa))).*log(cos(4.*kappa.*tpcal/Tp)/cos(kappa))));
AHP =fliplr(revAHP);

BIR4s1 = revAHP;
BIR4s2obs = horzcat(AHP.*exp(1i.*deltaphiobs),revAHP.*exp(1i.*deltaphiobs));
BIR4s2decHinder = horzcat(AHP.*exp(1i.*deltaphidec),revAHP.*exp(1i.*deltaphidec));
BIR4s3Hinder = AHP.*exp(1i.*pi);
BIR4180 = horzcat(revAHP,AHP.*exp(1i.*deltaphi180),revAHP.*exp(1i.*deltaphi180),AHP);

end