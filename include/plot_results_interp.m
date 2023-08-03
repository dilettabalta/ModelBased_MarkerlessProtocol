function plot_results_interp(d_fol,~,results)
%Author: Diletta Balta
%Department of Electronics and Telecommunication
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides all the frames of the gait cycle with LM, LE and GT
%positions superimposed on the images ('Result_interp' folder)

%inputs
%d_fol = folder containing the dynamic acquisition
%results = a MATLAB structure containing LM, LE and GT positions for each template (static, loading and
%swing models)


close all;
set(0, 'DefaultFigureVisible', 'off')
warning off

inizio = 1;
dim1_f = 720; % image height
dim2_f = 1280; % width

savefolder = [d_fol '\Results_interp\']; % destination folder

fc = 7;
[p1,p2] = butter(4,fc*2/30);
status = mkdir(savefolder);

if status
    mkdir(savefolder);
end

listu = ls([d_fol '\Video']);
lis = zeros(1,size(listu,1)-3);
cont = 1;

%list containing the number of frames of the gait cycle
for jj=3:size(listu,1)
    for jjj=1:size(listu(jj,:),2)
        if strcmp(listu(jj,jjj),'.')
            lis(1,cont)=str2num(listu(jj,1:(jjj-1)));
            cont=cont+1;
        end
    end
end

lis=sort(lis,'ascend');

%% interpolate models's contributions based on the percentage of the gait cycle by using a nonlinear sinusoid weight function
frameLoad = results.frameLoad;
frameFlex = results.frameFlex;
frameSta = results.frameSta;

ML_load_pl = results.load.ML_pl;
ML_flex_pl = results.flex.ML_pl;
ML_sta_pl = results.static.ML_pl;

ML_load_imm = results.load.ML_imm;
ML_flex_imm = results.flex.ML_imm;
ML_sta_imm = results.static.ML_imm;

KN_load_pl = results.load.KN_pl;
KN_flex_pl = results.flex.KN_pl;
KN_sta_pl = results.static.KN_pl;

KN_load_imm = results.load.KN_imm;
KN_flex_imm = results.flex.KN_imm;
KN_sta_imm = results.static.KN_imm;

TR_load_pl = results.load.TR_pl;
TR_flex_pl = results.flex.TR_pl;
TR_sta_pl = results.static.TR_pl;

TR_load_imm = results.load.TR_imm;
TR_flex_imm = results.flex.TR_imm;
TR_sta_imm = results.static.TR_imm;

ankle = results.angle.ankle;
knee = results.angle.knee;
hip = results.angle.hip;

%% LATERAL MALLEOLUS
%% plot coordinates

beta = results.modulationParam.beta';
b1 = results.modulationParam.b1';
b2 = results.modulationParam.b2';
b3 = results.modulationParam.b3';
indB1 = results.modulationParam.indB1';
indB2 = results.modulationParam.indB2';
indB3 = results.modulationParam.indB3';
[mod] = doTripleInterp(ML_load_pl(1,:),ML_sta_pl(1,:),ML_flex_pl(1,:),results.modulationParam,frameSta);
ML_mod_pl_x = mod;
[mod] = doTripleInterp(ML_load_pl(2,:),ML_sta_pl(2,:),ML_flex_pl(2,:),results.modulationParam,frameSta);
ML_mod_pl_y = mod;
ML_mod_pl = [ML_mod_pl_x;ML_mod_pl_y];

%% LATERAL MALLEOLUS
%% image coordinates

[mod] = doTripleInterp(ML_load_imm(1,:),ML_sta_imm(1,:),ML_flex_imm(1,:),results.modulationParam,frameSta);
ML_mod_imm_x = filtfilt(p1,p2,mod);
[mod] = doTripleInterp(ML_load_imm(2,:),ML_sta_imm(2,:),ML_flex_imm(2,:),results.modulationParam,frameSta);
ML_mod_imm_y = filtfilt(p1,p2,mod);
ML_mod_imm = [ML_mod_imm_x;ML_mod_imm_y];

%% FOOT ANGLES
foot_ang_load = results.load.foot_angles;
foot_ang_static = results.static.foot_angles;
foot_ang_flex = results.flex.foot_angles;

foot_ang = doTripleInterp(foot_ang_load,foot_ang_static,foot_ang_flex,results.modulationParam,frameSta);
foot_ang = filtfilt(p1,p2,foot_ang);

%% LATERAL EPICONDYLE
%% plot coordinates
[mod] = doTripleInterp(KN_load_pl(1,:),KN_sta_pl(1,:),KN_flex_pl(1,:),results.modulationParam,frameSta);
KN_mod_pl_x = filtfilt(p1,p2,mod);
[mod] = doTripleInterp(KN_load_pl(2,:),KN_sta_pl(2,:),KN_flex_pl(2,:),results.modulationParam,frameSta);
KN_mod_pl_y = filtfilt(p1,p2,mod);
KN_mod_pl = [KN_mod_pl_x;KN_mod_pl_y];

%% image coordinates
[mod] = doTripleInterp(KN_load_imm(1,:),KN_sta_imm(1,:),KN_flex_imm(1,:),results.modulationParam,frameSta);
KN_mod_imm_x = filtfilt(p1,p2,mod);
[mod] = doTripleInterp(KN_load_imm(2,:),KN_sta_imm(2,:),KN_flex_imm(2,:),results.modulationParam,frameSta);
KN_mod_imm_y = filtfilt(p1,p2,mod);
KN_mod_imm = [KN_mod_imm_x;KN_mod_imm_y];

%% GREAT TROCANTHER
%% plot coordinates
[mod] = doTripleInterp(TR_load_pl(1,:),TR_sta_pl(1,:),TR_flex_pl(1,:),results.modulationParam,frameSta);
TR_mod_pl_x = filtfilt(p1,p2,mod);
[mod] = doTripleInterp(TR_load_pl(2,:),TR_sta_pl(2,:),TR_flex_pl(2,:),results.modulationParam,frameSta);
TR_mod_pl_y = filtfilt(p1,p2,mod);
TR_mod_pl = [TR_mod_pl_x;TR_mod_pl_y];

%% image coordinates
[mod] = doTripleInterp(TR_load_imm(1,:),TR_sta_imm(1,:),TR_flex_imm(1,:),results.modulationParam,frameSta);
TR_mod_imm_x = filtfilt(p1,p2,mod);
[mod] = doTripleInterp(TR_load_imm(2,:),TR_sta_imm(2,:),TR_flex_imm(2,:),results.modulationParam,frameSta);
TR_mod_imm_y = filtfilt(p1,p2,mod);
TR_mod_imm = [TR_mod_imm_x;TR_mod_imm_y];

for lu = 1:length(lis)
    fr=lis(lu);
    if fr<10
        rgb = imread([d_fol '\Video\000' num2str(fr) '.bmp']);
    elseif fr>=10 && fr<100
        rgb = imread([d_fol '\Video\00' num2str(fr) '.bmp']);
    elseif fr>=100
        rgb = imread([d_fol '\Video\0' num2str(fr) '.bmp']);
    end
    figure
    imshow(rgb)
    hold on
    plot(ML_mod_imm(2,lu),ML_mod_imm(1,lu),'or','Linewidth',2)
    plot(KN_mod_pl(2,lu),abs(KN_mod_pl(1,lu)-dim1_f)+1,'or','Linewidth',2)
    line([KN_mod_pl(2,lu) ML_mod_imm(2,lu)], [abs(KN_mod_pl(1,lu)-dim1_f)+1 ML_mod_imm(1,lu)],'Color','k','Linewidth',2)
    plot(TR_mod_pl(2,lu),abs(TR_mod_pl(1,lu)-dim1_f)+1,'or','Linewidth',2)
    line([TR_mod_pl(2,lu) KN_mod_imm(2,lu)], [abs(TR_mod_pl(1,lu)-dim1_f)+1 KN_mod_imm(1,lu)],'Color','k','Linewidth',2)
    foot_imm = results.load.foot1_imm{lu};
    [r_foot, c_foot] = find(foot_imm);
    m_foot = tan(deg2rad(foot_ang(lu)));
    x_rect_foot = min(c_foot):max(c_foot);
    y_rect_foot = m_foot*(x_rect_foot-ML_mod_pl(2,lu))+ML_mod_pl(1,lu);
    y_rect_foot_imm = abs(y_rect_foot-dim1_f)+1;
    ind = find(min(r_foot)-50<y_rect_foot_imm & y_rect_foot_imm<max(r_foot));
    y_rect_foot_imm = y_rect_foot_imm(min(r_foot)-50<y_rect_foot_imm & y_rect_foot_imm<max(r_foot));
    plot(x_rect_foot(ind),y_rect_foot_imm,'k','Linewidth',2)
    print('-dpng','-r600',[savefolder 'res_' num2str(fr)]);
    close all
end
