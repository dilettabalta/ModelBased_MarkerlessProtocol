function results = kinematic_evaluationV2(d_fol,side_type,events,results)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides the sagittal lower-limb joint angles following the methods described in the methodological paper (paragraph
%JOINT KINEMATICS ESTIMATION).
% For each joint, three kinematic curves were obtained based on the three sets of body segment templates.
% These curves were then combined into a single curve by a nonlinear sinusoid weight function based on the percentage of the gait cycle

%inputs
%d_fol = folder containing the dynamic acquisition
%side_type = R (right) or L (left)
%events = vector containing [first frame, last frame, loading frame, stance
%frame, swing frame]
%results = a MATLAB structure containing LM, LE and GT positions for each template (static, loading and
%swing models)

%outputs
%results = a MATLAB structure containing ankle, knee and hip kinematics

set(0,'DefaultFigureVisible','on')
time = 0:100;

switch side_type
    case 'R'
        a = 1;
        b = 0;
    case 'L'
        a = -1; %on left trials, angles are specular than right ones
        b = 180;%on left trials, angles are specular than right ones
end


%% Computing foot angle from markerless analysis

events_t = events;
start = events_t(1) + 1;
frameFlex_t = events_t(5)+2-start;
frameLoad_t = events_t(3)+2-start;
frameSta_t = events_t(4)+2-start;
foot_mksSta = results.static.foot_angles';
foot_mksFlex = results.flex.foot_angles';
foot_mksLoad = results.load.foot_angles';

foot_mksSta = foot_mksSta*a+b;%left trial need the specular angle
foot_mksFlex = foot_mksFlex*a+b;%left trial need the specular angle
foot_mksLoad = foot_mksLoad*a+b;%left trial need the specular angle

ind_drops = 0;
fc = 7;
fs = 30;
[p1,p2] = butter(4,fc*2/fs);

ind_drop_rit = ind_drops - start + 1;
ind2cut = ind_drop_rit(ind_drop_rit>=1 & ind_drop_rit < length(foot_mksFlex));%frame sbagliati
x_mks = 1:length(foot_mksFlex)+length(ind2cut);
xq_mks = 1:length(foot_mksFlex)+length(ind2cut);
n_drop = ind2cut(ind2cut>1&ind2cut<length(foot_mksFlex));
new_smp = nan(1,length(n_drop));
for cc=1:length(n_drop)
    new_smp(cc) = ind2cut(cc)+length(ind2cut(1:cc))-1;
end

frameFlex = frameFlex_t + sum(new_smp<=frameFlex_t);
frameLoad = frameLoad_t + sum(new_smp<=frameLoad_t);
frameSta = frameSta_t + sum(new_smp<=frameSta_t);
results.frameFlex = frameFlex;
results.frameLoad = frameLoad;
results.frameSta = frameSta;
x_mks(new_smp) = [];
foot_mksLoad_clean = interp1(x_mks,foot_mksLoad,xq_mks,'spline')';
foot_mksFlex_clean = interp1(x_mks,foot_mksFlex,xq_mks,'spline')';
foot_mksSta_clean = interp1(x_mks,foot_mksSta,xq_mks,'spline')';

% computing shank angles
shank_mksFlex = results.flex.shank_angles';
shank_mksLoad = results.load.shank_angles';
shank_mksSta = results.static.shank_angles';
shank_mksFlex = shank_mksFlex*a+b; %left trial need the specular angle
shank_mksLoad = shank_mksLoad*a+b; %left trial need the specular angle
shank_mksSta = shank_mksSta*a+b; %left trial need the specular angle

shank_mksFlex_clean = interp1(x_mks,shank_mksFlex,xq_mks,'spline')';
shank_mksLoad_clean = interp1(x_mks,shank_mksLoad,xq_mks,'spline')';
shank_mksSta_clean = interp1(x_mks,shank_mksSta,xq_mks,'spline')';

% computing thigh angle
thigh_mksFlex = results.flex.thigh_angles';
thigh_mksLoad = results.load.thigh_angles';
thigh_mksSta = results.static.thigh_angles';

thigh_mksFlex = thigh_mksFlex*a+b;
thigh_mksLoad = thigh_mksLoad*a+b;
thigh_mksSta = thigh_mksSta*a+b;

thigh_mksFlex_clean = interp1(x_mks,thigh_mksFlex,xq_mks,'spline')';
thigh_mksLoad_clean = interp1(x_mks,thigh_mksLoad,xq_mks,'spline')';
thigh_mksSta_clean = interp1(x_mks,thigh_mksSta,xq_mks,'spline')';

%% computing knee angles
knee_mksFlex = thigh_mksFlex_clean-shank_mksFlex_clean;
knee_mksLoad = thigh_mksLoad_clean-shank_mksLoad_clean;
knee_mksSta = thigh_mksSta_clean-shank_mksSta_clean;

beta = mean([knee_mksLoad,knee_mksSta,knee_mksFlex],2);
b1 = beta(frameLoad);
b2 = beta(frameSta);
b3 = beta(frameFlex);
indB1 = find(beta(1:frameSta)>=b1);
indB2 = find(beta<=b2);
indB3 = frameSta + find(beta(frameSta:end)>=b3) - 1;

results.modulationParam.beta = beta;
results.modulationParam.b1 = b1;
results.modulationParam.b2 = b2;
results.modulationParam.b3 = b3;
results.modulationParam.indB1 = indB1;
results.modulationParam.indB2 = indB2;
results.modulationParam.indB3 = indB3;

%% non-linear sinusoid weight function for combining the three curves (one for each set of body set templates) into a single one
knee_mks_modStart = knee_mksSta(1:frameSta)+0.5*(knee_mksLoad(1:frameSta)-knee_mksSta(1:frameSta)).*(1-sin(pi*(beta(1:frameSta)-b2)/(b1-b2)+0.5*pi));
knee_mks_modEnd = knee_mksFlex(frameSta+1:end)+0.5*(knee_mksSta(frameSta+1:end)-knee_mksFlex(frameSta+1:end)).*(1-sin(pi*(beta(frameSta+1:end)-b3)/(b2-b3)+0.5*pi));
knee_mks_mod = [knee_mks_modStart;knee_mks_modEnd];
knee_mks_mod(indB1) = knee_mksLoad(indB1);
knee_mks_mod(indB2) = knee_mksSta(indB2);
knee_mks_mod(indB3) = knee_mksFlex(indB3);

knee_mks_fin = filtfilt(p1,p2,knee_mks_mod);
knee_mks_fin100 = interp1(0:length(knee_mks_fin)-1,knee_mks_fin,linspace(0,length(knee_mks_fin)-1,101),'spline')';
knee_mks_fin100 = roundd(knee_mks_fin100,1);
angle.knee = knee_mks_fin100;
knee_mksRaw = knee_mks_fin;
knee_mksRaw(new_smp) = [];
angle.kneeRaw = knee_mksRaw;
K1 = knee_mks_fin100(1);
K2 = max(knee_mks_fin100(1:41));
[K3,indK3_t] = min(knee_mks_fin100(26:76));
indK3 = indK3_t + 25;
K4 = mean(knee_mks_fin100(indK3-5:indK3+5));
K5 = max(knee_mks_fin100(51:end));
K6 = max(knee_mks_fin100) - min(knee_mks_fin100);
parameters = round([K1 K2 K3 K4 K5 K6]);
results.parameters = parameters;

%% computing hip angle (thigh angle with respect to the y-axis)
hip_mksFlex = thigh_mksFlex_clean - 90;
hip_mksSta = thigh_mksSta_clean - 90;
hip_mksLoad = thigh_mksLoad_clean - 90;

%% non-linear sinusoid weight function for combining the three curves (one for each set of body set templates) into a single one
hip_mks_modStart = hip_mksSta(1:frameSta)+0.5*(hip_mksLoad(1:frameSta)-hip_mksSta(1:frameSta)).*(1-sin(pi*(beta(1:frameSta)-b2)/(b1-b2)+0.5*pi));
hip_mks_modEnd = hip_mksFlex(frameSta+1:end)+0.5*(hip_mksSta(frameSta+1:end)-hip_mksFlex(frameSta+1:end)).*(1-sin(pi*(beta(frameSta+1:end)-b3)/(b2-b3)+0.5*pi));
hip_mks_mod = [hip_mks_modStart;hip_mks_modEnd];
hip_mks_mod(indB1) = hip_mksLoad(indB1);
hip_mks_mod(indB2) = hip_mksSta(indB2);
hip_mks_mod(indB3) = hip_mksFlex(indB3);
hip_mks_fin = filtfilt(p1,p2,hip_mks_mod);
hip_mks_fin100 = interp1(0:length(hip_mks_fin)-1,hip_mks_fin,linspace(0,length(hip_mks_fin)-1,101),'spline')';
hip_mks_fin100 = roundd(hip_mks_fin100,1);
angle.hip = hip_mks_fin100;
hip_mksRaw = hip_mks_fin;
hip_mksRaw(new_smp) = [];
angle.hipRaw = hip_mksRaw;

%% computing ankle angles
ankle_mksFlex = 90 - shank_mksFlex_clean + foot_mksFlex_clean;
ankle_mksSta = 90 - shank_mksSta_clean + foot_mksSta_clean;
ankle_mksLoad = 90 - shank_mksLoad_clean + foot_mksLoad_clean;

%% non-linear sinusoid weight function for combining the three curves (one for each set of body set templates) into a single one
ankle_mks_modStart = ankle_mksSta(1:frameSta)+0.5*(ankle_mksLoad(1:frameSta)-ankle_mksSta(1:frameSta)).*(1-sin(pi*(beta(1:frameSta)-b2)/(b1-b2)+0.5*pi));
ankle_mks_modEnd = ankle_mksFlex(frameSta+1:end)+0.5*(ankle_mksSta(frameSta+1:end)-ankle_mksFlex(frameSta+1:end)).*(1-sin(pi*(beta(frameSta+1:end)-b3)/(b2-b3)+0.5*pi));
ankle_mks_mod = [ankle_mks_modStart;ankle_mks_modEnd];
ankle_mks_mod(indB1) = ankle_mksLoad(indB1);
ankle_mks_mod(indB2) = ankle_mksSta(indB2);
ankle_mks_mod(indB3) = ankle_mksFlex(indB3);
ankle_mks_fin = filtfilt(p1,p2,ankle_mks_mod);
ankle_mks_fin100 = interp1(0:length(ankle_mks_fin)-1,ankle_mks_fin,linspace(0,length(ankle_mks_fin)-1,101),'spline')';
ankle_mks_fin100 = roundd(ankle_mks_fin100,1);
angle.ankle = ankle_mks_fin100;
ankle_mksRaw = ankle_mks_mod;
ankle_mksRaw(new_smp) = [];
angle.ankleRaw = ankle_mksRaw;
results.angle = angle;

%% figures
drawnow
h = figure;
set(gcf,'Position',[1 1 1366 685]) %fullscreen plot
plot(time,hip_mks_fin100,'b')
legend('mkl')
title('hip angle')
xlabel('%')
ylabel('degree')
savename1 = ['hip.tif'];
saveas(h,fullfile(d_fol,savename1));
%
% ankle figure
h = figure;
set(gcf,'Position',[1 1 1366 685]) %fullscreen plot
plot(time,ankle_mks_fin100,'b')
title('ankle angle')
legend('mkl')
xlabel('%')
ylabel('degree')
savename1 = ['ankle.tif'];
saveas(h,fullfile(d_fol,savename1));
%
% knee figure
h = figure;
set(gcf,'Position',[1 1 1366 685]) %fullscreen plot
plot(time,knee_mks_fin100,'b')
legend('mkl')
title('knee angle')
xlabel('%')
ylabel('degree')
savename1 = ['knee.tif'];
saveas(h,fullfile(d_fol,savename1));
delete([d_fol 'kinematics.xlsx']);
xlswrite([d_fol,'kinematics.xlsx'],ankle_mks_fin100','Ankle','A3')
xlswrite([d_fol,'kinematics.xlsx'],knee_mks_fin100','Knee','A3')
xlswrite([d_fol,'kinematics.xlsx'],hip_mks_fin100','Hip','A3')

end
