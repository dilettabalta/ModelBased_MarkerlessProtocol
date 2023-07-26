%% MARKERLESS DATA PROCESSING


%d_fol1 = input(prompt,'s');
sub = str2double(d_fol(end-33));
side = d_fol(end-18);

%x = inputdlg({'Dynamic Acquisitions Folder path','Static Acquisitions Folder path',},'Markerless Data Processing',[1 100; 1 100]);
%d_fol = x(1);  %[cd '\Alba_Rago_Diletta_Balta_trial_1_D_L_21_03_2022_17_16\'];
%d_fol1 = x(2); %[cd '\Alba_Rago_Diletta_Balta_trial_0_S_L_21_03_2022_17_13\'];
%side = questdlg('Select side of the foreground leg','Markerless Data Processing','L','R','');

%% Creazione cartelle color, depth, cloud in cui saranno spostati i file
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


set(0, 'DefaultFigureVisible', 'off') %disabilito la visualizzazione delle figure (quando faccio girare i codici)
%lo start e lo stop devono essere definiti dall'algoritmo
%dell'identificazione del tallone e della punta

try
    load ([d_fol '\start.mat']);
    load ([d_fol '\stop.mat']);
catch
    disp('Select the first frame of the gait cycle:')
    cd ([d_fol 'color_stream']);
    start = uigetfile('*.png');
    start = str2double(regexp(start,'\d*','Match'));

    disp('Select the last frame of the gait cycle:')
    stop = uigetfile('*.png');
    stop = str2double(regexp(stop,'\d*','Match'));
    cd ..
    save start start
    save stop stop
end

B = imread([d_fol 'color_stream\color0.png']); %frame rappresentante il background
crop_seq_txt_Azure(d_fol,start+1,stop+1); %allineamento spaziale tra RGB e Depth in modo da creare 2 cartelle ('VIDEO' e 'MATCH_d_raw')
mkdir([d_fol 'Misc\Segm']); %cartella che conterrà le segmentazioni
lis = ls([d_fol 'Video']);
frames = size(lis,1);
frames = frames-2;
%% Segmentation - Dynamic trial
 subject_segmentation(d_fol,frames,B,side,sub);
          
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
 subject_segmentation(d_fol1,frames,B,side,sub);

%% Models definition
if d_fol(10) == 'r'
    cd('C:\Users\rgalb\Dropbox (Politecnico Di Torino Studenti)\Tesi - Alba Rago - Azure Kinect\Acquisitions\Codici')
else
    cd('C:\Users\balta\Dropbox\Tesi - Alba Rago - Azure Kinect (1)\Acquisitions\Codici')
end
disp('Static Model')
template.static = dynamic_calib_ok_provaAzure(d_fol1,0,side,'static',sub);

try
    load ([d_fol '\ref_flex.mat']);
    load ([d_fol '\ref_load.mat']);
    load ([d_fol '\stance.mat']);
%     disp([ref_load ref_flex]);
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
    stance = start + str2double(regexp(stance,'\d*','Match'));

    cd ..
    save ref_load ref_load
    save ref_flex ref_flex
    save stance stance
    cd ..
end
disp('Load Model')
template.load = dynamic_calib_ok_provaAzure(d_fol,ref_load,side,'load',sub);
disp('Flex Model')
template.flex = dynamic_calib_ok_provaAzure(d_fol,ref_flex,side,'flex',sub);

%% Kinematic estimation
listu = lis;
lis2=zeros(1,size(listu,1)-2);
cont=1;
%riempe lis con il numero di ogni foto
for jj=3:size(listu,1)
    for jjj=1:size(listu(jj,:),2)
        if strcmp(listu(jj,jjj),'.')
            lis2(1,cont)=str2num(listu(jj,1:(jjj-1)));
            cont=cont+1;
        end
    end
end
lis2=sort(lis2,'ascend');
addpath('C:\Users\balta\Dropbox\Tesi - Alba Rago - Azure Kinect(1)\Acquisitions\Codici')
display('Kinematic estimations in:')
res_static_f = foot_matching_sc(d_fol,side,lis2,template.static,sub);
res_static_f_s = shank_matching_Azure_nuova(d_fol,lis2,template.static,res_static_f);
res_static_f_s_t = thigh_matching_Azure(d_fol,lis2,template.static,res_static_f_s);
results.static = res_static_f_s_t;
disp ('static')

res_load_f = foot_matching_sc(d_fol,side,lis2,template.load,sub);
res_load_f_s = shank_matching_Azure_nuova(d_fol,lis2,template.load,res_load_f);
res_load_f_s_t = thigh_matching_Azure(d_fol,lis2,template.load,res_load_f_s);
results.load = res_load_f_s_t;
disp ('load')
res_flex_f = foot_matching_sc(d_fol,side,lis2,template.flex,sub);
res_flex_f_s = shank_matching_Azure_nuova(d_fol,lis2,template.flex,res_flex_f);
res_flex_f_s_t = thigh_matching_Azure(d_fol,lis2,template.flex,res_flex_f_s);
results.flex = res_flex_f_s_t;
disp ('flex')

events = [start, stop, ref_load+start,stance,ref_flex+start];
[results,beta,b1,b2,new_smp] = kinematic_evaluationV2(d_fol,side,events,results);
plot_results_interp(d_fol,lis2,results);

%% SDK BODY TRACKING PROCESSING
cd('C:\Users\balta\Dropbox\Tesi - Alba Rago - Azure Kinect (1)\Codici\SDK')
results_sdk = json_bodytracking(d_fol, side, start, stop);
close all

%% generazione di un video

writer = VideoWriter([cd '\AzureKinect_MLM.mp4']);
writer.FrameRate = 10;
open(writer)

for i = 1:43
   %img = imread([ d_fol, '\SDK Results\leg\leg' num2str(i-1) '.jpg']);
   img = imread([ d_fol, '\Results_interp\res_' num2str(i-1) '.png']);
   writeVideo(writer,img);
end
close(writer)
