function mod = doTripleInterp(a, b, c, param, f)
%Author: Diletta Balta
%Department of Electronics and Telecommunication
%Politecnico di Torino 
%diletta.balta@polito.it

%This function provides LM, LE and GT positions obtained by weighting each template contribution by a nonlinear sinusoid weight function based on the percentage of the gait cycle

%inputs 
%a = anatomical landmark positions obtained by the loading template
%b = anatomical landmark positions obtained by the static template
%c = anatomical landmark positions obtained by the swingg template
%param = modulation parameters
%f = stance frame

%outputs 
%mod = final anatomical landmark positions obtained by weighting each template contribution by a nonlinear sinusoid weight function based on the percentage of the gait cycle

beta = param.beta';
b1 = param.b1';
b2 = param.b2';
b3 = param.b3';
indB1 = param.indB1';
indB2 = param.indB2';
indB3 = param.indB3';
modStart = b(1,1:f)+0.5*(a(1,1:f)-b(1,1:f)).*(1-sin(pi*(beta(1:f)-b2)/(b1-b2)+0.5*pi));
modEnd = c(1,f+1:end)+0.5*(b(1,f+1:end)-c(1,f+1:end)).*(1-sin(pi*(beta(f+1:end)-b3)/(b2-b3)+0.5*pi));
mod = [modStart,modEnd];
mod(indB1) = a(1,indB1);
mod(indB2) = b(1,indB2);
mod(indB3) = c(1,indB3);
end

