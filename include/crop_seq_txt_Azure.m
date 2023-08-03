function crop_seq_txt_Azure(d_fol,start,stop,frame)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function creates 2 folders ('VIDEO' e 'MATCH_d_raw') containing RGB and Depth images of the analyzed gait cycle

%inputs
%d_fol = folder containing the dynamic acquisition
%start = first frame of the gait cycle
%stop = last frame of the gait cycle

o_fol_rect = [d_fol 'Video\'];
o_fol_match_d1 = [d_fol 'Misc\MATCH_d_raw\'];

status_rect = mkdir(o_fol_rect);
status_match_d1 = mkdir(o_fol_match_d1);

if status_rect
    mkdir(o_fol_rect);
end
if status_match_d1
    mkdir(o_fol_match_d1);
end

list_temp = ls([d_fol 'color_stream']);
list_temp_d = ls([d_fol 'depth_stream']);
type = 'g';
len = size(list_temp,2);

for k=len:-1:1
    ind_l = find(list_temp(:,k)==type);
    if isempty(ind_l)
        break
    end
    if k==len
        list = list_temp(ind_l,:);
        list_d = list_temp_d(ind_l,:);
    else
        list = [list_temp(ind_l,:);list];
        list_d = [list_temp_d(ind_l,:);list_d];
    end
end

%% Saving of the RGB and depth images of the selected gait cycle
if ~exist('frame','var')
    fr = 0;
else
    fr = frame;
    start = frame+start;
    stop = start;
end

for l1 = start:stop
    drawnow
    rgb_dist = imread([d_fol 'color_stream\' list(l1,:)]);
    depth_dist = readAzureBin([d_fol 'depth_stream\' list_d(l1,:)]);

    if fr<=9
        p_name1_rect = [o_fol_rect '000' num2str(fr) '.bmp'];
        p_name1_match_d1 = [o_fol_match_d1 '000' num2str(fr) '.txt'];

    elseif fr>=10 && fr<100
        p_name1_rect = [o_fol_rect '00' num2str(fr) '.bmp'];
        p_name1_match_d1=[o_fol_match_d1 '00' num2str(fr) '.txt'];

    elseif fr>=100
        p_name1_rect=[o_fol_rect '0' num2str(fr) '.bmp'];
        p_name1_match_d1=[o_fol_match_d1 '0' num2str(fr) '.txt'];
    end

    %Video
    imwrite(rgb_dist,p_name1_rect,'bmp');

    %MATCH_d_raw
    im_match_d1 = depth_dist;
    file_match_d1 = fopen(p_name1_match_d1,'wt');
    for ii = 1:size(im_match_d1,1)
        fprintf(file_match_d1,'%d\t',im_match_d1(ii,:));
        fprintf(file_match_d1,'\n');
    end
    fclose(file_match_d1);
    fr=fr+1;
end

close all

end

