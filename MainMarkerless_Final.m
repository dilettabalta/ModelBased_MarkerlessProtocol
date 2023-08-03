%% MARKERLESS DATA PROCESSING
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino 
%diletta.balta@polito.it

clear all
close all
clc

warning off
addpath('D:\Codici_MLM_final\include')
pathname = 'D:\Codici_MLM_final'; %path containing the MATLAB codes

% FOLDERS with the dynamic and the static acquisitions
disp('Select the dynamic acquisition folder:')
d_fol = [uigetdir() '\'];
disp('Select the static acquisition folder:')
d_fol1 = [uigetdir() '\'];

sub = str2double(d_fol(end-33));
side = d_fol(end-18);

%% Creation of three folders that will contain RGB and Depth and 3D point clouds 
cd(d_fol)
mkdir color_stream
try
    movefile 'color*.png' color_stream
catch
    disp('Color folder already exists')
end

mkdir depth_stream
try
    movefile 'depth*.bin' depth_stream
catch
    disp('Depth folder already exists')
end

mkdir cloud
try
    movefile 'cloud*.bin' cloud
catch
    disp('Cloud folder already exists or there are no cloud files')
end
cd ..

set(0, 'DefaultFigureVisible', 'off') 

try
    load ([d_fol '\start.mat']);
    load ([d_fol '\stop.mat']);

catch
    disp('Select the first frame (both feet must be fully visible):')
    cd ([d_fol 'color_stream']);
    start = uigetfile('*.png');
    start = str2double(regexp(start,'\d*','Match'));

    disp('Select the last frame (both feet must be fully visible):')
    stop = uigetfile('*.png');
    stop = str2double(regexp(stop,'\d*','Match'));
    cd ..
    save start start
    save stop stop
end

%% Gait cycle identification
conv_factor = 0.002; %see appendix A of the methodological paper for more details on how to set this parameter
[start1, stop1, start2, pos_IC1, pos_IC2, pos_IC1_c, stride_length, step_length, gait_speed] = gaitcycle(d_fol, start, stop, side, conv_factor);
cd ..
save start1 start1
save stop1 stop1
B = imread([d_fol 'color_stream\color0.png']); %Background frame (B)
crop_seq_txt_Azure(d_fol,start1+1,stop1+1); %creation of 2 folders ('VIDEO' e 'MATCH_d_raw') containing RGB and Depth images of the analyzed gait cycle
mkdir([d_fol 'Misc\Segm']); %folder containing the segmentation masks for each gait cycle
lis = ls([d_fol 'Video']);
frames = size(lis,1);
frames = frames-2;

%% Segmentation - Dynamic trial
subject_segmentation(d_fol,frames,B,side,sub); %function for subject segmentation in dynamic frames
          
%% Segmentation - Static trial
cd(d_fol1)
mkdir color_stream
try
    movefile 'color*.png' color_stream
catch
    disp('Color folder already exists')
end

mkdir depth_stream
try
    movefile 'depth*.bin' depth_stream
catch
    disp('Depth folder already exists')
end

mkdir cloud
try
    movefile 'cloud*.bin' cloud
catch
    disp('Cloud folder already exists or there are no cloud files')
end
cd ..

addpath(d_fol1)
crop_seq_txt_Azure(d_fol1,1,1);
mkdir([d_fol1 '\Misc\Segm'])
lis1 = ls([d_fol1 '\Video']);
frames = size(lis1,1);
frames = frames-2;

%% Segmentation - Static trial
set(0, 'DefaultFigureVisible', 'off')
subject_segmentation(d_fol1,frames,B,side,sub); %function for subject segmentation 
%% Subject-specific models definition
cd(pathname)
disp('Static Model')
template.static = dynamic_calib_ok_provaAzure_v1(d_fol1,0,side,'static',sub); %creazione del modello lower-limb 2D in statica

try
    load ([d_fol '\ref_flex.mat']);
    load ([d_fol '\ref_load.mat']);
    load ([d_fol '\stance.mat']);

catch
    cd ([d_fol '\Video'])
    disp('Select the LOAD frame of the gait cycle:')
    ref_load = uigetfile('*.bmp');
    ref_load = str2double(regexp(ref_load,'\d*','Match'));

    disp('Select the FLEX frame of the gait cycle:')
    ref_flex = uigetfile('*.bmp');
    ref_flex = str2double(regexp(ref_flex,'\d*','Match'));

    disp('Select the STANCE frame of the gait cycle:')
    stance = uigetfile('*.bmp');
    stance = start1 + str2double(regexp(stance,'\d*','Match'));

    cd ..
    save ref_load ref_load
    save ref_flex ref_flex
    save stance stance
    cd ..
end

disp('Load Model')
template.load = dynamic_calib_ok_provaAzure_v1(d_fol,ref_load,side,'load',sub); %creazione del modello lower-limb 2D in load
disp('Flex Model')
template.flex = dynamic_calib_ok_provaAzure_v1(d_fol,ref_flex,side,'flex',sub); %creazione del modello lower-limb 2D in swing

%%	JOINT CENTERS TRAJECTORIES ESTIMATION
listu = lis;
lis2=zeros(1,size(listu,1)-2);
cont=1;
%list containing the number of frames of the gait cycle 
for jj=3:size(listu,1)
    for jjj=1:size(listu(jj,:),2)
        if strcmp(listu(jj,jjj),'.')
            lis2(1,cont)=str2num(listu(jj,1:(jjj-1))); %#ok<ST2NM> 
            cont=cont+1;
        end
    end
end
lis2=sort(lis2,'ascend');
addpath('D:\Codici_MLM_final\include')
disp('Kinematic estimations in:')
res_static_f = foot_matching_sc(d_fol,side,lis2,template.static,sub); %foot fitting for estimating the lateral malleoulus through the static model
res_static_f_s = shank_matching_Azure_nuova_v1(d_fol,lis2,template.static,res_static_f);%shank fitting for estimating the lateral epicondyle through the static  model
res_static_f_s_t = thigh_matching_Azure_v1(d_fol,lis2,template.static,res_static_f_s);%thigh fitting for estimating the great throcanter through the static model
results.static = res_static_f_s_t;
disp ('static')

res_load_f = foot_matching_sc(d_fol,side,lis2,template.load,sub);%foot fitting for estimating the lateral malleoulus through the loading model
res_load_f_s = shank_matching_Azure_nuova_v1(d_fol,lis2,template.load,res_load_f);%shank fitting for estimating the lateral epicondyle through the loading model
res_load_f_s_t = thigh_matching_Azure_v1(d_fol,lis2,template.load,res_load_f_s);%thigh fitting for estimating the great throcanter through the loading model
results.load = res_load_f_s_t;
disp ('load')

res_flex_f = foot_matching_sc(d_fol,side,lis2,template.flex,sub);%foot fitting for estimating the lateral malleoulus through the swing model
res_flex_f_s = shank_matching_Azure_nuova_v1(d_fol,lis2,template.flex,res_flex_f);%shank fitting for estimating the lateral epicondyle through the swing model
res_flex_f_s_t = thigh_matching_Azure_v1(d_fol,lis2,template.flex,res_flex_f_s); %thigh fitting for estimating the great throcanter through the swing model
results.flex = res_flex_f_s_t;
disp ('flex')

events = [start1, stop1, ref_load+start1,stance,ref_flex+start1];
results = kinematic_evaluationV2(d_fol,side,events,results); %sagittal lower-limb joint angle estimation
plot_results_interp(d_fol,lis2,results);


