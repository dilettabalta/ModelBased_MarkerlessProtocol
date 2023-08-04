function template = dynamic_calib_ok_provaAzure_v1(d_fol,frame,side,Type,~)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides a MATLAB structure containing a set of body segments
%templates (foot template, shank template and thigh template) following the methods described in the methodological paper (paragraph
%ANATOMICAL CALIBRATION AND BODY SEGMENT TEMPLATES DEFINITION)

%inputs
%d_fol = folder containing the dynamic acquisition
%frame = static image/ loading image/ swing image
%side = R (right) or L (left)
%Type = 'static', 'load', 'flex'

%outputs
%template is MATLAB structure containing a set of body segments templates (foot template, shank template and thigh template)

show_save_all = true;
dim1_f = 720;
dim2_f = 1280;
thr_cut = 0.9;
thr_sole_start = 0.2;
thr_sole_end = 0.5;

if frame<10
    ref = ['000' num2str(frame)];
elseif frame<100
    ref = ['00' num2str(frame)];
elseif frame<1000
    ref = ['0' num2str(frame)];
end

save_fol = [d_fol 'Calib\'];
status_graphs=mkdir(save_fol);

if status_graphs
    mkdir(save_fol);
end

rgb_rect = imread([d_fol 'Video\' ref '.bmp']);
h = figure;
set(h, 'MenuBar', 'none');
set(h, 'ToolBar', 'none');
imshow(rgb_rect);
hold on

%% ANATOMICAL CALIBRATION AND BODY SEGMENT TEMPLATES DEFINITION
[x_joints, y_joints] = ginputc(4,'Color', 'r','ShowPoints',true);
x_test = x_joints(4);
y_test = y_joints(4);
x_mal = x_joints(1);
y_mal = y_joints(1);
x_gin = x_joints(2);
y_gin = y_joints(2);
x_tr = x_joints(3);
y_tr = y_joints(3);
MLRef=[y_mal, x_mal];
KNRef=[y_gin, x_gin];
TRRef=[y_tr, x_tr];
shank_length= pdist([x_mal,y_mal;x_gin, y_gin],'euclidean');
m_shank = -(y_mal-y_gin)/(x_mal-x_gin);
alpha_shank = atan(m_shank);
if alpha_shank<0
    alpha_shank = alpha_shank + deg2rad(180);
end

thigh_length=pdist([x_tr, y_tr;x_gin, y_gin],'euclidean');

m_thigh= -(y_tr-y_gin)/(x_tr-x_gin); %inverto il segno perchè sono coordinate immagine
alpha_thigh= atan(m_thigh);
if alpha_thigh<0
    alpha_thigh= alpha_thigh+ deg2rad(180);
end

up_lim = y_mal-(shank_length+thigh_length)*1.1;
template.up_lim = round(up_lim);
segm_rect = logical(imread([d_fol 'Misc\Segm\' num2str(ref) '.bmp']));

%estraggo subito il piede per avere distanza di riferimento
segm_rect_t = segm_rect;

[r_seg_t, ~] = find(segm_rect);
half_seg_r = (max(r_seg_t)+min(r_seg_t))/2;
segm_rect_t(1:half_seg_r,:) = false;
hsv_or = rgb2hsv(rgb_rect);
D_imm_r_t = false(dim1_f,dim2_f);
D_imm_l_t = false(dim1_f,dim2_f);
I = rgb_rect;

[D_imm_l_t, D_imm_r_t] = Feet_segmentation_dyn1 (I);
figure,imshowpair(D_imm_r_t,D_imm_l_t,'montage')
% figure,imshow(D_imm_r_t|D_imm_l_t)

% Feet segmentation refinement through morphological operators

if side == 'R'
    dim1_f = size(I,1);
    dim2_f = size(I,2);
    D1_imm_t = D_imm_r_t;
    D1_imm_t = bwmorph(D1_imm_t,'bridge','inf');
    D1_imm_t = imfill(D1_imm_t,'holes');
    D1_imm_t = bwmorph(D1_imm_t,'majority','inf');
    D1_el = bwconncomp(D1_imm_t,8);
    arsr_D1 = cellfun(@numel,D1_el.PixelIdxList);
    [~, X_D1] = max(arsr_D1);
    D1_imm = false(size(I,1),size(I,2));
    D1_imm(D1_el.PixelIdxList{X_D1})=1;
    foot1_rect = D1_imm;
    D2_imm_t = D_imm_l_t;
    D2_imm_t = bwmorph(D2_imm_t,'bridge','inf');
    D2_imm_t = imfill(D2_imm_t,'holes');
    D2_imm_t = bwmorph(D2_imm_t,'majority','inf');
    D2_el = bwconncomp(D2_imm_t,8);
    if ~isempty(D2_el.PixelIdxList)
        arsr_D2 = cellfun(@numel,D2_el.PixelIdxList);
        [~, X_D2] = max(arsr_D2);
        D2_imm = false(size(I,1),size(I,2));
        D2_imm(D2_el.PixelIdxList{X_D2})=1;
        foot2_rect = D2_imm;
    else
        foot2_rect = D2_imm_t;
    end

else
    dim1_f = size(I,1);
    dim2_f = size(I,2);
    D1_imm_t = D_imm_l_t;
    D1_imm_t = bwmorph(D1_imm_t,'bridge','inf');
    D1_imm_t = imfill(D1_imm_t,'holes');
    D1_imm_t = bwmorph(D1_imm_t,'majority','inf');
    D1_el = bwconncomp(D1_imm_t,8);
    arsr_D1 = cellfun(@numel,D1_el.PixelIdxList);
    [~, X_D1] = max(arsr_D1);
    D1_imm = false(size(I,1),size(I,2));
    D1_imm(D1_el.PixelIdxList{X_D1})=1;
    foot1_rect = D1_imm;
    D2_imm_t = D_imm_r_t;
    D2_imm_t = bwmorph(D2_imm_t,'bridge','inf');
    D2_imm_t = imfill(D2_imm_t,'holes');
    D2_imm_t = bwmorph(D2_imm_t,'majority','inf');
    D2_el = bwconncomp(D2_imm_t,8);
    if ~isempty(D2_el.PixelIdxList)
        arsr_D2 = cellfun(@numel,D2_el.PixelIdxList);
        [~, X_D2] = max(arsr_D2);
        D2_imm = false(size(I,1),size(I,2));
        D2_imm(D2_el.PixelIdxList{X_D2})=1;
        foot2_rect = D2_imm;
    else
        foot2_rect = D2_imm_t;
    end
end

D1_imm_t = bwmorph(D1_imm_t,'bridge','inf');
D1_imm_t = imfill(D1_imm_t,'holes');
D1_imm_t = bwmorph(D1_imm_t,'majority','inf');

D1_el = bwconncomp(D1_imm_t,8);
arsr_D1 = cellfun(@numel,D1_el.PixelIdxList);
[~, X_D1] = max(arsr_D1);

D1_imm = false(dim1_f,dim2_f);
D1_imm(D1_el.PixelIdxList{X_D1})=1;

foot1_rect = D1_imm;
D2_imm_t = bwmorph(D2_imm_t,'bridge','inf');
D2_imm_t = imfill(D2_imm_t,'holes');
D2_imm_t = bwmorph(D2_imm_t,'majority','inf');

try
    D2_el = bwconncomp(D2_imm_t,8);
    arsr_D2 = cellfun(@numel,D2_el.PixelIdxList);
    [~, X_D2] = max(arsr_D2);

    D2_imm = false(dim1_f,dim2_f);
    D2_imm(D2_el.PixelIdxList{X_D2}) = 1;
catch
    D2_imm = false(dim1_f,dim2_f);
end

foot2_rect = D2_imm;

segm_rect(1:y_tr,:) = false;
depth_match = readmatrix([d_fol '\Misc\MATCH_d_raw\'  num2str(ref) '.txt']);

%%contour
[e_r_in_T,e_c_in_T] = find(edge(foot1_rect,'canny'));

%inverting and scaling y-axis
e_rr_in_T = abs(e_r_in_T-dim1_f)+1;
e_cc_in_T = e_c_in_T;

e_rr_T = e_rr_in_T;
e_cc_T = e_cc_in_T;

rr_low_T = 0;
for i=min(e_cc_in_T):max(e_cc_in_T)
    rr_low_T(i-min(e_cc_in_T)+1) = min(e_rr_in_T(e_cc_in_T==i));
end
cc_low_T = min(e_cc_in_T):max(e_cc_in_T);

clearvars rr_in_low_T cc_in_low_T
for i = min(e_cc_in_T):max(e_cc_in_T)
    cand_T = e_rr_in_T(e_cc_in_T==i)';
    cand_cut_T = find(diff(cand_T)<-16,1,'last');
    if isempty(cand_cut_T)
        ok_cand_T = cand_T;
    else
        ok_cand_T = cand_T(cand_cut_T+1:end);
    end

    if exist('rr_in_low_T','var')
        rr_in_low_T = [rr_in_low_T ok_cand_T];
    else
        rr_in_low_T = ok_cand_T;
    end

    if exist('cc_in_low_T','var')
        cc_in_low_T = [cc_in_low_T ones(1,length(ok_cand_T))*i];
    else
        cc_in_low_T = ones(1,length(ok_cand_T))*i;
    end
end

if side == 'L'
    cc_low_T = fliplr(cc_low_T);
    rr_low_T = fliplr(rr_low_T);
end

%cutting the template computing the slope of the sole as the line
%best fitting the point between heel and metatarsus

ind_start_sole_T = round(thr_sole_start*length(cc_low_T));
ind_end_sole_T = round(thr_sole_end*length(cc_low_T));
rr_sole_T = rr_low_T(ind_start_sole_T:ind_end_sole_T);
cc_sole_T = cc_low_T(ind_start_sole_T:ind_end_sole_T);
fun_T = fit(cc_sole_T',rr_sole_T','poly1');
coeff_temp_T = coeffvalues(fun_T);
m_sole_T = coeff_temp_T(1);
q_sole_T = coeff_temp_T(2);

alpha_foot = atan(m_sole_T);
R0_1_T = [cos(alpha_foot) -sin(alpha_foot);sin(alpha_foot) cos(alpha_foot)];
foot_perim = [e_rr_T';e_cc_T'];
O1_1_T = -R0_1_T*[rr_low_T(10);cc_low_T(10)];
T0_1_T = zeros(3);
T0_1_T(1:2,1:2) = R0_1_T;
T0_1_T(:,3) = [O1_1_T;1];
data_1_T = T0_1_T*[foot_perim;ones(1,length(e_rr_T))];
data_1_T = data_1_T(1:2,:);

thr = 0.90;
check_heel_T = 0;
heel_T = 0;
while check_heel_T<=thr*length(e_cc_T)
    heel_T=heel_T-1;
    check_heel_T = sum(data_1_T(2,:)>heel_T);
end

check_toe_T = 0;
toe_T = 0;
while check_toe_T<=thr*length(e_cc_T)
    toe_T=toe_T+1;
    check_toe_T = sum(data_1_T(2,:)<toe_T);
end

if side == 'L'
    temp = toe_T;
    toe_T = heel_T;
    heel_T = temp;
end

foot_length = abs(toe_T-heel_T);
cut_length = foot_length*thr_cut;

switch side
    case 'R'
        cut_T = heel_T + cut_length;
        ind_to_cut_T = find(data_1_T(2,:)>cut_T);
    case 'L'
        cut_T = heel_T - cut_length;
        ind_to_cut_T = find(data_1_T(2,:)<cut_T);
end

foot_perim_cut = foot_perim;
foot_perim_cut(:,ind_to_cut_T) = [];

if show_save_all
    figure;plot(foot_perim_cut(2,:),foot_perim_cut(1,:),'or')
    axis equal
end

foot_perim_cut(1,:) = abs(foot_perim_cut(1,:)-dim1_f)+1;
foot_cut = false(size(foot1_rect));
ind_foot_cut_dist = mat2cell(foot_perim_cut', size(foot_perim_cut', 1), ones(1, size(foot_perim_cut', 2)));
foot_cut(sub2ind(size(foot_cut), ind_foot_cut_dist{:})) = true;
if show_save_all
    figure;imshow(foot_cut);
end
%foot_cut is the edge of the foot
[r_foot_rect, c_foot_rect] = find(foot_cut);

%% Shank identification
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = dim2_f;
imageSizeY = dim1_f;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = MLRef(2);
centerY = MLRef(1);
radius = shank_length;
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
segm_mask = segm_rect & circlePixels;
if show_save_all
    figure;imshow(segm_mask)
end

shank_low_cut_par = 1/5;
shank_high_cut_par = 4/5;
radius_ext = shank_high_cut_par*shank_length;
circlePixels_ext = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_ext.^2;
radius_int = shank_low_cut_par*shank_length;
circlePixels_int = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_int.^2;
annulus = (circlePixels_ext & not(circlePixels_int));

mask_comp = bwconncomp(segm_mask,8);
segm_mask = segm_mask & not(foot1_rect | foot2_rect);
segm_mask(MLRef(1):end,:) = false;
depth_segm_el = immultiply(depth_match,segm_mask);
depth_segm_el(depth_segm_el==0) = NaN;
depth_segm_el = round(depth_segm_el);
depth_match(depth_match == 0) = nan;
depth_segm_el(depth_segm_el == 0) = nan;
most_freq = nanmedian(depth_match(foot1_rect)); %#ok<NANMEDIAN> % median not mean to exclude outliers
lim_left = 2200;
lim_right = 2850;
depth_segm_el(depth_segm_el>lim_right) = NaN;
depth_segm_el(depth_segm_el<lim_left) = NaN;

h = histogram(depth_segm_el);

%% Otsu technique to identify the foregroung shank
level = multithresh(depth_segm_el,1);
binLocations=h.BinEdges';
binLocations(end)=[];
counts=h.Values';
first_clust=binLocations<=level(1); % First cluster
second_clust=binLocations>level(1) & binLocations<=lim_right; %Second cluster
[~, ind_firstpeak] = max(counts(first_clust)); % First cluster peak
[~, ind_secondpeak] = max(counts(second_clust)); % Second cluster peak
vettore_ind = [binLocations(ind_firstpeak) binLocations(find(second_clust,1)+ind_secondpeak-1)];
d_ind = diff(vettore_ind);
ind_merg = (d_ind<=30);
if ind_merg == 1
    shank_start = lim_left;
    shank_stop = lim_right;
else
    %figure()
    %histogram(depth_segm_el,'BinEdges',lim_left:lim_right)
    %xline(level(1),'--g','LineWidth',0.9)

    shank_start = lim_left;
    shank_stop = level;
end

shank_depth = depth_segm_el;
shank_depth(shank_depth<shank_start) = NaN;
shank_depth(shank_depth>shank_stop) = NaN;
shank_imm_t = segm_mask;
shank_imm_t(isnan(depth_segm_el)) = false;
shank_imm_t(depth_segm_el>shank_stop) = false;
shank_imm_t(depth_segm_el<shank_start) = false;
shank_depth = immultiply(shank_depth,annulus);
shank_imm_t = immultiply(shank_imm_t,annulus);
shank_imm_t = bwmorph(shank_imm_t,'bridge','inf');
shank_imm_t = imfill(shank_imm_t,'holes'); %Foreground shank

if show_save_all
    figure;imshow(shank_imm_t);
end

shank_comp = bwconncomp(shank_imm_t,8);
arsr_shank = cellfun(@numel,shank_comp.PixelIdxList);

[~, X_shank] = max(arsr_shank);
shank_imm_temp = false(size(shank_depth));
shank_imm_temp(shank_comp.PixelIdxList{X_shank})=1;

[r_shank_t, c_shank_t] = find(shank_imm_temp);
r_shank_t = abs(r_shank_t-dim1_f)+1;

shank_cdm_r = mean(r_shank_t);
shank_cdm_c = mean(c_shank_t);
m_ax_shank = (shank_cdm_r-abs(MLRef(1)-dim1_f))/(shank_cdm_c-MLRef(2));
q_ax_shank = abs(MLRef(1)-dim1_f)-m_ax_shank*MLRef(2);
m_ax_ort_shank = -1/m_ax_shank;
[cut_low_shank_X, cut_low_shank_Y] = linecirc(m_ax_shank,q_ax_shank,MLRef(2),abs(MLRef(1)-dim1_f)+1,shank_length*shank_low_cut_par);

[cut_low_shank(2), temp] = max(cut_low_shank_Y);
cut_low_shank(1) = cut_low_shank_X(temp);
q_ax_ort_shank = cut_low_shank(2)-m_ax_ort_shank*(cut_low_shank(1));

shank_low_control = r_shank_t > c_shank_t.*m_ax_ort_shank+q_ax_ort_shank;
r_shank_rect = abs(r_shank_t(shank_low_control)-dim1_f)+1;
c_shank_rect = c_shank_t(shank_low_control);

shank_rect = false(size(shank_imm_temp));
ind_shank = mat2cell([r_shank_rect c_shank_rect], length(r_shank_rect), ones(1, 2));
shank_rect(sub2ind(size(shank_rect), ind_shank{:})) = true; %shank template

%% creating a circle centered on forward thigh, in order to
% reduce depth data to elaborate
% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = dim2_f;
imageSizeY = dim1_f;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the circle in the image.
centerX = KNRef(2);
centerY = KNRef(1);
radius = thigh_length*(1);
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
segm_mask = segm_rect & circlePixels & not(foot1_rect | foot2_rect);
segm_mask(KNRef(1):end,:) = false;
depth_segm_el = immultiply(depth_match,segm_mask);

depth_segm_el(depth_segm_el==0) = NaN;
depth_segm_el = round(depth_segm_el);

shank_depth = immultiply(shank_rect,depth_match);
most_freq = nanmedian(nanmedian(shank_depth(shank_depth~=0))); %#ok<NANMEDIAN>

lim_left = 2300;
lim_right = 2850;
depth_segm_el(depth_segm_el>lim_right) = NaN;
depth_segm_el(depth_segm_el<lim_left) = NaN;

h = histogram(depth_segm_el,lim_left:lim_right);

%% Otsu technique to identify the foregroung thigh
level = multithresh(depth_segm_el,2);
% figure
% histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('First histogram with 2 thresholds'), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
% set(gca,'Color',[0.50, 0.50, 0.50]) % Colore sfondo
% hold on
% xline(level(1),'--g','LineWidth',0.9)
% xline(level(2),'--r','LineWidth',0.9)

binLocations=h.BinEdges';
binLocations(end)=[];
counts=h.Values';
first_clust=binLocations<=level(1); % First cluster
second_clust=binLocations>level(1) & binLocations<=level(2); % Second cluster
third_cluster = binLocations>level(2) & binLocations<=lim_right; % Third cluster
[First_peak, ind_firstpeak] = max(counts(first_clust)); % First cluster peak
[Second_peak, ind_secondpeak] = max(counts(second_clust)); % Second cluster peak
[Third_peak, ind_thirdpeak] = max(counts(third_cluster)); %Third cluster peak
vettore_ind = [binLocations(ind_firstpeak) binLocations(find(second_clust,1)+ind_secondpeak) binLocations(find(third_cluster,1)+ind_thirdpeak-1)];
d_ind = diff(vettore_ind);
ind_merg = d_ind<=100;
vettore_peak = [First_peak Second_peak Third_peak];
d_peak = diff(vettore_peak);
peak_merg = abs(d_peak)<=180;
inters = ind_merg & peak_merg;

if nnz(inters)==0 %in the circle there are the foreground hand, the foreground and background thighs
    %     figure()
    %     histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('Depth histogram'), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
    %     set(gca,'Color',[0.50, 0.50, 0.50]) % Colore sfondo
    %     hold on
    %     xline(level(1),'--g','LineWidth',0.9)
    %     xline(level(2),'--r','LineWidth',0.9)
    %     legend('','First Otsu thereshold','Second Otsu thereshold')
    idx = depth_segm_el>level(1) & depth_segm_el<level(2);
    thigh_hand_start = lim_left;
    thigh_hand_stop = level(2);
    thigh_start = level(1);
    thigh_stop = level(2);
else
    ind_min = find(inters);
    level = multithresh(depth_segm_el,1);
    first_clust=binLocations<=level(1);
    second_clust=binLocations>level(1);
    [First_peak, ~] = max(counts(first_clust)); % First cluster peak
    [Second_peak, ~] = max(counts(second_clust)); % Second cluster peak
    if First_peak>Second_peak
        %         figure()
        %         histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('Depth histogram '), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
        %         set(gca,'Color',[0.50, 0.50, 0.50])
        %         hold on
        %         xline(level,'--r','LineWidth',0.9)
        %         legend('','Otsu thereshold')
        if sum(counts(first_clust))>sum(counts(second_clust)) %the first cluster is the foreground thigh
            idx=depth_segm_el<level;
            thigh_hand_start = lim_left;
            thigh_hand_stop = lim_right;
            thigh_start = lim_left;
            thigh_stop = level;
        else %the second cluster is the foreground thigh (the first one is the foreground hand)
            thigh_hand_start = lim_left;
            thigh_hand_stop = lim_right;
            thigh_start = level;
            thigh_stop = lim_right;
        end
    else
        if sum(counts(first_clust))>sum(counts(second_clust)) %the first cluster is the foreground thigh
            idx=depth_segm_el<level;
            thigh_hand_start = lim_left;
            thigh_hand_stop = lim_right;
            thigh_start = lim_left;
            thigh_stop = level;
        else %the second cluster is the foreground thigh
            thigh_hand_start = lim_left;
            thigh_hand_stop = lim_right;
            thigh_start = level;
            thigh_stop = lim_right;
        end
    end
end

if strcmp(Type,'static')
    level = multithresh(depth_segm_el,1);
    binLocations=h.BinEdges';
    binLocations(end)=[];
    counts=h.Values';
    first_clust=binLocations<=level(1); % First cluster
    second_clust=binLocations>level(1) & binLocations<=lim_right; %Second cluster
    [~, ind_firstpeak] = max(counts(first_clust)); % First cluster peak
    [~, ind_secondpeak] = max(counts(second_clust)); % Second cluster peak
    vettore_ind = [binLocations(ind_firstpeak) binLocations(find(second_clust,1)+ind_secondpeak-1)];
    d_ind = diff(vettore_ind);
    ind_merg = d_ind<=100;
    if ind_merg == 1
        thigh_hand_start = lim_left;
        thigh_hand_stop = lim_right;
        thigh_start = lim_left;
        thigh_stop = lim_right;
    else
        %         figure()
        %         histogram(depth_segm_el,'BinEdges',lim_left:lim_right)
        %         xline(level(1),'--g','LineWidth',0.9)
        thigh_hand_start = lim_left;
        thigh_hand_stop = level;
        thigh_start = lim_left;
        thigh_stop = level;
    end
end


thigh_low_cut_par = 1/6;
thigh_high_cut_par = 2/3;

radius_ext = thigh_high_cut_par*thigh_length;
circlePixels_ext = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_ext.^2;
radius_int = thigh_low_cut_par*thigh_length;
circlePixels_int = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_int.^2;
annulus = (circlePixels_ext & not(circlePixels_int));

thigh_depth = depth_segm_el;
thigh_depth(thigh_depth<thigh_start) = NaN;
thigh_depth(thigh_depth>thigh_stop) = NaN;
thigh_imm_t = segm_mask;
thigh_imm_t(isnan(depth_segm_el)) = false;
thigh_imm_t(depth_segm_el>thigh_stop) = false;
thigh_imm_t(depth_segm_el<thigh_start) = false;

thigh_depth = immultiply(thigh_depth,annulus);
thigh_imm_t = immultiply(thigh_imm_t,annulus);
thigh_imm_t = bwmorph(thigh_imm_t,'bridge','inf');
thigh_imm_t= imfill(thigh_imm_t,'holes');

thigh_comp = bwconncomp(thigh_imm_t,8);
arsr_thigh = cellfun(@numel,thigh_comp.PixelIdxList);

[~, X_thigh] = max(arsr_thigh);
thigh_imm_temp = false(size(thigh_depth));
thigh_imm_temp(thigh_comp.PixelIdxList{X_thigh})=1;

[r_thigh_t, c_thigh_t] = find(thigh_imm_temp);
r_thigh_t = abs(r_thigh_t-dim1_f)+1;

thigh_cdm_r = mean(r_thigh_t);
thigh_cdm_c = mean(c_thigh_t);
m_ax_thigh = (thigh_cdm_r-abs(KNRef(1)-dim1_f))/(thigh_cdm_c-KNRef(2));
q_ax_thigh = abs(KNRef(1)-dim1_f)-m_ax_thigh*KNRef(2);
m_ax_ort_thigh = -1/m_ax_thigh;
[cut_low_thigh_X, cut_low_thigh_Y] = linecirc(m_ax_thigh,q_ax_thigh,KNRef(2),abs(KNRef(1)-dim1_f)+1,thigh_length*thigh_low_cut_par);
[cut_low_thigh(2), temp] = max(cut_low_thigh_Y);
cut_low_thigh(1) = cut_low_thigh_X(temp);
q_ax_ort_thigh = cut_low_thigh(2)-m_ax_ort_thigh*(cut_low_thigh(1));
thigh_low_control = r_thigh_t > c_thigh_t.*m_ax_ort_thigh+q_ax_ort_thigh;
r_thigh_rect = abs(r_thigh_t(thigh_low_control)-dim1_f)+1;
c_thigh_rect = c_thigh_t(thigh_low_control);
thigh_rect = false(size(thigh_imm_temp));
ind_thigh = mat2cell([r_thigh_rect c_thigh_rect], length(r_thigh_rect), ones(1, 2));
thigh_rect(sub2ind(size(thigh_rect), ind_thigh{:})) = true; %thigh template
template_rect = (foot_cut | shank_rect | thigh_rect); %set of body segments template

if show_save_all
    figure;imshow(template_rect);
    hold on
    plot(MLRef(2),MLRef(1),'or','Linewidth',2);
    hold on
    plot(KNRef(2),KNRef(1),'or','Linewidth',2);
    hold on
    plot(TRRef(2),TRRef(1),'or','Linewidth',2);
end

template.frame = frame+1;
template.mal = MLRef;
template.knee = KNRef;
template.troc = TRRef;
template.foot = foot_cut;
template.shank = shank_rect;
template.thigh = thigh_rect;
template.shank_length = shank_length;
template.thigh_length = thigh_length;
template.cut = cut_length;
template.m_sole = m_sole_T;
template.q_sole = q_sole_T;
switch side
    case 'R'
        template.alpha_foot = alpha_foot;
        template.angle_foot = rad2deg(template.alpha_foot);
        template.alpha_shank = alpha_shank;
        template.angle_shank = rad2deg(template.alpha_shank);
        template.alpha_thigh = alpha_thigh;
        template.angle_thigh = rad2deg(template.alpha_thigh);
        %         save ([d_fol 'template'], 'template', '-v7.3');
    case 'L'
        template.alpha_foot = alpha_foot + deg2rad(180);
        template.angle_foot = rad2deg(template.alpha_foot);
        template.alpha_shank = alpha_shank;
        template.angle_shank = rad2deg(template.alpha_shank);
        template.alpha_thigh = alpha_thigh;
        template.angle_thigh = rad2deg(template.alpha_thigh);
        %         save ([d_fol 'template'], 'template', '-v7.3');
end

set(0, 'DefaultFigureVisible', 'on')
h = figure;
set(h, 'MenuBar', 'none');
set(h, 'ToolBar', 'none');
imshow(rgb_rect)
hold on
plot(c_foot_rect,r_foot_rect,'.b')
hold on
plot(c_shank_rect,r_shank_rect,'.b')
hold on
plot(c_thigh_rect,r_thigh_rect,'.b')
hold on
plot(MLRef(2),MLRef(1),'or','Linewidth',2);
hold on
plot(KNRef(2),KNRef(1),'or','Linewidth',2);
hold on
plot(TRRef(2),TRRef(1),'or','Linewidth',2);
hold on
line([MLRef(2) KNRef(2)], [MLRef(1) KNRef(1)],'Color','g','Linewidth',2)
hold on
line([KNRef(2) TRRef(2)], [KNRef(1) TRRef(1)],'Color','g','Linewidth',2)
set(h,'units','normalized','outerposition',[0 0 1 1])
% print('-dpng','-r600',[d_fol 'Calib\template_' type]);
saveas(h,[d_fol 'Calib\template_' Type '.bmp']);
close all;
