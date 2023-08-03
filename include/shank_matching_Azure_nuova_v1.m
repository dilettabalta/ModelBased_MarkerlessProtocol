function results = shank_matching_Azure_nuova_v1(d_fol,lis,template,results)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides LE positions in each frame by matching the shank template to each shank mask by using ICP following the methods described in the methodological paper
% paragraph JOINT CENTERS TRAJECTORIES ESTIMATION - knee joint center (KJC) estimation).

%inputs
%d_fol = folder containing the dynamic acquisition
%lis = list containing the number of frames of the gait cycle
%template = static/loading/swing shank template

%outputs
%results = a MATLAB structure containing LE position on each frame

close all;
set(0, 'DefaultFigureVisible', 'off')
dim1_f = 720;
dim2_f = 1280;
inizio = 1;
show_save_all = false;
save_imp = false;

if save_imp
    mkdir([d_fol 'matching_shank\']);
    if exist([d_fol 'matching_shank\visual\'],'dir')==7
        rmdir([d_fol 'matching_shank\visual\'],'s');
    end
    mkdir([d_fol 'matching_shank\visual\']);
end

ML_imm = results.ML_imm;
ML_pl = results.ML_pl;
foot1_rect = results.foot1_imm;
foot2_rect = results.foot2_imm;

shank_imm_T = template.shank;

[r_shank_T, c_shank_T] = find(shank_imm_T);
template_pl = [abs(r_shank_T-dim1_f)+1, c_shank_T]';

list_res = zeros(1,length(lis));
list_KN_pl = zeros(2,length(lis));
list_KN_imm = zeros(2,length(lis));
list_angles = zeros(1,length(lis));
clearvars shank_match shank_data

for lu =  1:length(lis)

    drawnow
    scaling = 1;
    KN_imm_T = template.knee';
    KN_pl_T = scaling*[abs(KN_imm_T(1)-dim1_f)+1; KN_imm_T(2)];
    ML_imm_T = template.mal';
    ML_pl_T = scaling*[abs(ML_imm_T(1)-dim1_f)+1; ML_imm_T(2)];

    rigid_source_t = scaling*[template_pl(2,:); template_pl(1,:)];
    rigid_KN_t = [KN_pl_T(2);KN_pl_T(1)];
    cdm_rigid_source_t = mean(rigid_source_t,2);
    m_ax1_rigid_source_t = (cdm_rigid_source_t(2)-ML_pl_T(1))/(cdm_rigid_source_t(1)-ML_pl_T(2));
    rigid_source_t = [rigid_source_t rigid_KN_t];

    fr = lis(lu);
    disp(fr)
    if fr<10
        depth_match = load([d_fol 'Misc\MATCH_d_raw\000' num2str(fr) '.txt']);
    elseif fr>=10 && fr<100
        depth_match = load([d_fol 'Misc\MATCH_d_raw\00' num2str(fr) '.txt']);
    elseif fr>=100
        depth_match = load([d_fol 'Misc\MATCH_d_raw\0' num2str(fr) '.txt']);
    end

    if fr<10
        segm_rect = logical(imread([d_fol 'Misc\Segm\000' num2str(fr) '.bmp']));
    elseif fr>=10 && fr<100
        segm_rect = logical(imread([d_fol 'Misc\Segm\00' num2str(fr) '.bmp']));
    elseif fr>=100
        segm_rect = logical(imread([d_fol 'Misc\Segm\0' num2str(fr) '.bmp']));
    end

    depth_match_large = depth_match;

    %distance between feet
    foot_dist_1 = nanmedian(depth_match_large(foot1_rect{lu}));
    foot_dist_2 = nanmedian(depth_match_large(foot2_rect{lu}));
    list_foot_dist_1(lu) = foot_dist_1;
    list_foot_dist_2(lu) = foot_dist_2;

    depth_segm = immultiply(depth_match_large,segm_rect);
    if show_save_all
        h5 = figure;imshow(depth_segm,[]);
        hold on
        plot(ML_imm(2,lu),ML_imm(1,lu),'or','Linewidth',2)
    end

    shank_length = scaling*template.shank_length;

    % Create a logical image of a circle with specified
    % diameter, center, and image size.
    % First create the image.
    imageSizeX = dim2_f;
    imageSizeY = dim1_f;
    [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    % Next create the circle in the image.
    centerX = ML_imm(2,lu);
    centerY = ML_imm(1,lu);
    radius = shank_length*(1);
    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    segm_mask = segm_rect & circlePixels;
    
    if show_save_all
        figure;imshow(segm_mask)
    end

    low_cut_par = 1/5;
    high_cut_par = 4/5;
    radius_ext = high_cut_par*shank_length;
    circlePixels_ext = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_ext.^2;
    radius_int = low_cut_par*shank_length;
    circlePixels_int = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_int.^2;
    annulus = (circlePixels_ext & not(circlePixels_int));

    segm_mask_temp = bwmorph(segm_mask,'spur');
    segm_mask_temp = bwmorph(segm_mask_temp,'clean');
    mask_comp = bwconncomp(segm_mask_temp,8);

    segm_mask = segm_mask & not(foot1_rect{lu} | foot2_rect{lu});
    segm_mask(ML_imm(1,lu):end,:) = false;
    depth_segm_el = immultiply(depth_match_large,segm_mask);
    depth_segm_el(depth_segm_el==0) = NaN;
    depth_segm_el = round(depth_segm_el);

    if show_save_all
        figure;imshow(depth_segm_el,[]);
    end

    depth_match(depth_match == 0) = nan;
    depth_segm_el(depth_segm_el == 0) = nan;
    most_freq = nanmedian(depth_match(foot1_rect{lu})); %#ok<NANMEDIAN> % median not mean to exclude outliers
    lim_left = 2000;
    lim_right = 2850;
    depth_segm_el(depth_segm_el>lim_right) = NaN;
    depth_segm_el(depth_segm_el<lim_left) = NaN;

    %% Identification of foregroung shank through Otsu technique
    h = histogram(depth_segm_el);
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
        figure()
        histogram(depth_segm_el,'BinEdges',lim_left:lim_right)
        xline(level(1),'--g','LineWidth',0.9)

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
    shank_imm_t=bwmorph(shank_imm_t,'bridge','inf');
    shank_imm_t= imfill(shank_imm_t,'holes');

    if show_save_all
        figure;imshow(shank_imm_t);
    end

    shank_comp = bwconncomp(shank_imm_t,8);
    arsr_shank = cellfun(@numel,shank_comp.PixelIdxList);

    [~, X_shank] = max(arsr_shank);
    shank_imm_temp = false(size(shank_depth));
    shank_imm_temp(shank_comp.PixelIdxList{X_shank})=1;

    if show_save_all
        saveas(h1,[d_fol 'matching_shank\visual\' template.type '_hist_' num2str(fr)],'bmp');
    end

    %h=figure;imshow(depth_match,[]), saveas(h,['depth_' num2str(lu-1) '.jpg']);
    %h=figure;imshow(shank_imm_temp), saveas(h,['shank static' num2str(lu-1) '.jpg']);
    [r_shank_t, c_shank_t] = find(shank_imm_temp);
    r_shank_t = abs(r_shank_t-dim1_f)+1; 
    shank_cdm_r = mean(r_shank_t);
    shank_cdm_c = mean(c_shank_t);
    m_ax_shank = (shank_cdm_r-ML_pl(1,lu))/(shank_cdm_c-ML_pl(2,lu));
    q_ax_shank = ML_pl(1,lu)-m_ax_shank*ML_pl(2,lu);
    m_ax_ort_shank = -1/m_ax_shank;
    [cut_low_shank_X, cut_low_shank_Y] = linecirc(m_ax_shank,q_ax_shank,ML_pl(2,lu),ML_pl(1,lu),shank_length*low_cut_par);

    [cut_low_shank(2), temp] = max(cut_low_shank_Y);
    cut_low_shank(1) = cut_low_shank_X(temp);
    q_ax_ort_shank = cut_low_shank(2)-m_ax_ort_shank*(cut_low_shank(1));

    shank_low_control = r_shank_t > c_shank_t.*m_ax_ort_shank+q_ax_ort_shank;
    r_shank_imm = abs(r_shank_t(shank_low_control)-dim1_f)+1;
    c_shank_imm = c_shank_t(shank_low_control);

    shank_imm = false(size(shank_imm_temp));
    ind_shank = mat2cell([r_shank_imm c_shank_imm], length(r_shank_imm), ones(1, 2));
    shank_imm(sub2ind(size(shank_imm), ind_shank{:})) = true;
    list_shank_imm{lu} = shank_imm;

    [r_shank, c_shank] = find(shank_imm);

    if show_save_all
        h_test = figure;imshow(segm_rect);
        hold on
        plot(c_shank,r_shank,'.r')
        hold on
        plot(ML_imm(2,lu),ML_imm(1,lu),'*r','Linewidth',2)
        %         saveas(h_test,[d_fol 'matching_shank\visual\' template.type '_plot_' num2str(fr)],'bmp');
    end

    %% full icp
    shank_pl = [abs(r_shank-dim1_f)+1, c_shank]';
    rigid_target_t = [shank_pl(2,:); shank_pl(1,:)];
    cdm_rigid_target_t = mean(rigid_target_t,2); %centroid computation
    m_ax1_rigid_target_t = (cdm_rigid_target_t(2)-ML_pl(1,lu))/(cdm_rigid_target_t(1)-ML_pl(2,lu));

    alpha_source_t = atan(m_ax1_rigid_source_t);
    if alpha_source_t<0
        alpha_source_t = alpha_source_t + deg2rad(180);
    end

    alpha_target_t = atan(m_ax1_rigid_target_t);
    if alpha_target_t<0
        alpha_target_t = alpha_target_t + deg2rad(180);
    end

    alpha = alpha_target_t - alpha_source_t;
    rotaz = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
    rigid_source_rot = rotaz*(rigid_source_t-repmat(cdm_rigid_source_t,1,size(rigid_source_t,2)))+repmat(cdm_rigid_source_t,1,size(rigid_source_t,2));
    rigid_KN_rot = rotaz*(rigid_KN_t-cdm_rigid_source_t)+cdm_rigid_source_t;
    ML_pl_rot_T = rotaz*([ML_pl_T(2);ML_pl_T(1)]-cdm_rigid_source_t)+cdm_rigid_source_t;

    rigid_source_rot_tr = rigid_source_rot - repmat(ML_pl_rot_T,1,size(rigid_source_rot,2));
    rigid_KN_rot_tr = rigid_KN_rot - ML_pl_rot_T;
    rigid_target_tr = rigid_target_t - repmat([ML_pl(2,lu);ML_pl(1,lu)],1,size(rigid_target_t,2));

    if show_save_all
        figure;plot(rigid_target_tr(1,:),rigid_target_tr(2,:),'or')
        hold on;plot(rigid_source_rot_tr(1,:),rigid_source_rot_tr(2,:),'ob')
        hold on;plot(rigid_KN_rot_tr(1),rigid_KN_rot_tr(2),'og')
        hold on;plot(0,0,'og')
        hold on;line([0 rigid_KN_rot_tr(1)], [0 rigid_KN_rot_tr(2)], 'Color','g','Linewidth',2)
        axis equal
    end

    [~,~,template_fin_t,shank_fin_t,res] = icp_2_rot(rigid_target_tr,rigid_source_rot_tr);

    template_fin = template_fin_t + repmat([ML_pl(2,lu);ML_pl(1,lu)],1,size(template_fin_t,2));
    KN_fin = template_fin(:,end);
    template_fin(:,end) = [];
    shank_match{lu} = template_fin;
    shank_fin = shank_fin_t + repmat([ML_pl(2,lu);ML_pl(1,lu)],1,size(shank_fin_t,2));
    shank_data{lu} = shank_fin;

    cdm_template_fin = mean(template_fin,2);
    m_ax1_template_fin = (cdm_template_fin(2)-ML_pl(1,lu))/(cdm_template_fin(1)-ML_pl(2,lu));

    alpha_template_fin = atan(m_ax1_template_fin);
    if alpha_template_fin<0
        alpha_template_fin = alpha_template_fin + deg2rad(180);
    end

    list_KN_pl(:,lu) = flipud(KN_fin)';
    list_KN_imm(:,lu) = [abs(KN_fin(2)-dim1_f)+1; KN_fin(1)];
    list_angles(lu) = rad2deg(atan2(KN_fin(2)-ML_pl(1,lu),KN_fin(1)-ML_pl(2,lu)));

    list_res(lu) = res;
    if show_save_all
        set(0, 'DefaultFigureVisible', 'off')
        h_icp1 = figure;
        plot(template_fin(1,:),template_fin(2,:),'ob');
        hold on;
        plot(shank_fin(1,:),shank_fin(2,:),'or');
        hold on
        plot(ML_pl(2,lu),ML_pl(1,lu),'og','Linewidth',2)
        hold on
        plot(KN_fin(1),KN_fin(2),'og','Linewidth',2)
        hold on
        line([KN_fin(1) ML_pl(2,lu)], [KN_fin(2) ML_pl(1,lu)],'Color','g','Linewidth',2)
        axis equal
    end

    edge_template_fin = edge_plot(template_fin);
    edge_shank_fin = edge_plot(shank_fin);
    if save_imp
        h_icp2 = figure;
        plot(edge_template_fin(1,:),edge_template_fin(2,:),'ob');
        hold on;
        plot(edge_shank_fin(1,:),edge_shank_fin(2,:),'or');
        hold on
        plot(ML_pl(2,lu),ML_pl(1,lu),'og','Linewidth',2)
        hold on
        plot(KN_fin(1),KN_fin(2),'og','Linewidth',2)
        hold on
        line([KN_fin(1) ML_pl(2,lu)], [KN_fin(2) ML_pl(1,lu)],'Color','g','Linewidth',2)
        axis equal
        %         saveas(h_icp2,[d_fol 'matching_shank\' template.type '_match_' num2str(fr)],'bmp');
    end

    if show_save_all
        h3 = figure;
        imshow(shank_imm);
    end

    close all
end

results.KN_pl = list_KN_pl;
results.KN_imm = list_KN_imm;
results.shank_pl = shank_match;
results.shank_imm = list_shank_imm;
results.shank_angles = list_angles;
results.shank_data_pl = shank_data;
results.shank_res = list_res;
results.foot_dist_1 = list_foot_dist_1;
results.foot_dist_2 = list_foot_dist_2;
save ([d_fol 'results'], 'results', '-v7.3');

if save_imp
    h_res = figure;
    plot(list_res,'b')
    hold on
    plot(list_res,'*b')
    %     saveas(h_res,[d_fol 'matching_shank\' template.type '_res'],'bmp');

    h_ang = figure;
    plot(list_angles,'b')
    hold on
    plot(list_angles,'*b')
    %     saveas(h_ang,[d_fol 'matching_shank\' template.type '_angles'],'bmp');
end

close all

end