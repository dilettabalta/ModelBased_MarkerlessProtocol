function results = thigh_matching_Azure_v1(d_fol,lis,template,results)
%Author: Diletta Balta
%Department of Electronics and Telecommunication
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides GT positions in each frame by matching the shank template to each shank mask by using ICP following the methods described in the methodological paper
% paragraph JOINT CENTERS TRAJECTORIES ESTIMATION - hip joint center (HJC) estimation).

%inputs
%d_fol = folder containing the dynamic acquisition
%lis = list containing the number of frames of the gait cycle
%template = static/loading/swing thigh template

%outputs
%results = a MATLAB structure containing GT position on each frame

close all;
set(0, 'DefaultFigureVisible', 'off')


dim1_f = 720;
dim2_f = 1280;

inizio = 1;
show_save_all = false;
save_imp = false;

if save_imp
    mkdir([d_fol 'matching_thigh\']);

    if exist([d_fol 'matching_thigh\visual\'],'dir')==7
        rmdir([d_fol 'matching_thigh\visual\'],'s');
    end
    mkdir([d_fol 'matching_thigh\visual\']);
end

KN_imm = results.KN_imm;
KN_pl = results.KN_pl;
foot1_rect = results.foot1_imm;
foot2_rect = results.foot2_imm;

thigh_imm_T = template.thigh;
[r_thigh_T, c_thigh_T] = find(thigh_imm_T);
template_pl = [abs(r_thigh_T-dim1_f)+1, c_thigh_T]';

list_res = zeros(1,length(lis));
list_TR_pl = zeros(2,length(lis));
list_TR_imm = zeros(2,length(lis));
list_angles = zeros(1,length(lis));
clearvars thigh_match thigh_data

for lu = 1:length(lis)
    disp(lu)
    scaling = 1;
    TR_imm_T = template.troc';
    TR_pl_T = scaling*[abs(TR_imm_T(1)-dim1_f)+1; TR_imm_T(2)];
    KN_imm_T = template.knee';
    KN_pl_T = scaling*[abs(KN_imm_T(1)-dim1_f)+1; KN_imm_T(2)];
    rigid_source_t = scaling*[template_pl(2,:); template_pl(1,:)];
    rigid_TR_t = [TR_pl_T(2);TR_pl_T(1)];
    cdm_rigid_source_t = mean(rigid_source_t,2);
    m_ax1_rigid_source_t = (cdm_rigid_source_t(2)-KN_pl_T(1))/(cdm_rigid_source_t(1)-KN_pl_T(2));
    rigid_source_t = [rigid_source_t rigid_TR_t];

    fr = lis(lu);
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
    depth_segm = immultiply(depth_match_large,segm_rect);
    if show_save_all
        h5 = figure;imshow(depth_segm,[]);
        hold on
        plot(KN_imm(2,lu),KN_imm(1,lu),'or','Linewidth',2)
    end

    thigh_length = scaling*template.thigh_length;

    %% creating a circle centered on forward thigh, in order to
    % reduce depth data to elaborate
    % Create a logical image of a circle with specified
    % diameter, center, and image size.
    % First create the image.
    imageSizeX = dim2_f;
    imageSizeY = dim1_f;
    [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    % Next create the circle in the image.
    centerX = KN_imm(2,lu);
    centerY = KN_imm(1,lu);
    radius = thigh_length*(1);
    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
    segm_mask = segm_rect & circlePixels & not(foot1_rect{lu} | foot2_rect{lu});
    segm_mask(KN_imm(1,lu):end,:) = false;
    depth_segm_el = immultiply(depth_match_large,segm_mask);
    depth_segm_el(depth_segm_el==0) = NaN;
    depth_segm_el = round(depth_segm_el);

    if lu==inizio
        shank_depth = immultiply(results.shank_imm{lu},depth_match_large);
        % median not mean to exclude outliers
        most_freq = nanmedian(nanmedian(shank_depth(shank_depth~=0)));
    else % taking as most frequently distance the one resulting from previous frame
        most_freq = nanmedian(nanmedian(shank_depth(shank_depth~=0)));
    end

    lim_left = 2300;
    lim_right = 2750;
    depth_segm_el(depth_segm_el>lim_right) = NaN;
    depth_segm_el(depth_segm_el<lim_left) = NaN;
    lim_right = max(max(depth_segm_el));
    lim_left = min(min(depth_segm_el));

    %% Identification of foreground thigh through Otsu technique
    h = histogram(depth_segm_el,lim_left:lim_right);
    level = multithresh(depth_segm_el,2);

    %     figure()
    %     histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('First histogram with 2 thresholds'), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
    %     set(gca,'Color',[0.50, 0.50, 0.50])
    %     hold on
    %     xline(level(1),'--g','LineWidth',0.9)
    %     xline(level(2),'--r','LineWidth',0.9)

    binLocations=h.BinEdges';
    binLocations(end)=[];
    counts=h.Values';
    first_clust=binLocations<=level(1); % First cluster
    second_clust=binLocations>level(1) & binLocations<=level(2); % Second cluster
    third_cluster = binLocations>level(2) & binLocations<=lim_right; %
    [First_peak, ind_firstpeak] = max(counts(first_clust)); % First cluster peak
    [Second_peak, ind_secondpeak] = max(counts(second_clust)); % Second cluster peak
    [Third_peak, ind_thirdpeak] = max(counts(third_cluster)); % Third cluster peak
    vettore_ind = [binLocations(ind_firstpeak) binLocations(find(second_clust,1)+ind_secondpeak) binLocations(find(third_cluster,1)+ind_thirdpeak-1)];
    d_ind = diff(vettore_ind);
    ind_merg = (d_ind<=80);
    vettore_peak = [First_peak Second_peak Third_peak];
    d_peak = diff(vettore_peak);
    peak_merg = (abs(d_peak)<=180);
    inters = ind_merg & peak_merg;

    if nnz(inters)==0 %in the circle there are the foreground hand, the foreground and background thighs
        %         figure()
        %         histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('Depth histogram'), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
        %         set(gca,'Color',[0.50, 0.50, 0.50]) % Colore sfondo
        %         hold on
        %         xline(level(1),'--g','LineWidth',0.9)
        %         xline(level(2),'--r','LineWidth',0.9)
        %         legend('','First Otsu thereshold','Second Otsu thereshold')
        idx = depth_segm_el>level(1) & depth_segm_el<level(2); % Ricerca dei valori che sono all'interno dell'intervello tra le soglie
        thigh_hand_start = lim_left;
        thigh_hand_stop = level(2);
        thigh_start = level(1);
        thigh_stop = level(2);
    else
        level = multithresh(depth_segm_el,1);
        first_clust=binLocations<=level(1);
        second_clust=binLocations>level(1);
        [First_peak, ind_firstpeak] = max(counts(first_clust)); % First cluster peak
        [Second_peak, ind_secondpeak] = max(counts(second_clust)); % Second cluster peak
        if First_peak>Second_peak
            %             figure()
            %             histogram(depth_segm_el,'BinEdges',lim_left:lim_right),title('Depth histogram '), xlabel('Distance from camera (mm)'),ylabel('Pixels occurences')
            %             set(gca,'Color',[0.50, 0.50, 0.50])
            %             hold on
            %             xline(level,'--r','LineWidth',0.9)
            %             legend('','Otsu thereshold')
            if sum(counts(first_clust))>sum(counts(second_clust))  %the first cluster is the foreground thigh
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


    if show_save_all
        mkdir ([d_fol '\Thigh_matching'])
        hold on;temp = xline(thigh_start,'--c');
        set(temp,'Linewidth',2)
        hold on;temp = xline(thigh_stop,'--k');
        set(temp,'Linewidth',2)
    end
    %     saveas(gcf,[ d_fol 'Thigh_matching\Thigh' num2str(fr) '.jpg']);

    low_cut_par = 1/6;
    high_cut_par = 2/3;

    radius_ext = high_cut_par*thigh_length;
    circlePixels_ext = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_ext.^2;
    radius_int = low_cut_par*thigh_length;
    circlePixels_int = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius_int.^2;
    annulus = (circlePixels_ext & not(circlePixels_int));

    thigh_depth = depth_segm_el;
    thigh_depth(thigh_depth<thigh_start) = NaN;
    thigh_depth(thigh_depth>thigh_stop) = NaN;

    thigh_imm_t = segm_mask;
    thigh_imm_t(isnan(depth_segm_el)) = false;
    thigh_imm_t(depth_segm_el>thigh_stop) = false;
    thigh_imm_t(depth_segm_el<thigh_start) = false;
    thigh_hand_depth = depth_segm_el;
    thigh_hand_depth(thigh_hand_depth<thigh_hand_start) = NaN;
    thigh_hand_depth(thigh_hand_depth>thigh_hand_stop) = NaN;
    thigh_hand_imm_t = segm_mask;
    thigh_hand_imm_t(isnan(depth_segm_el)) = false;
    thigh_hand_imm_t(depth_segm_el>thigh_hand_stop) = false;
    thigh_hand_imm_t(depth_segm_el<thigh_hand_start) = false;
    thigh_depth = immultiply(thigh_depth,annulus);
    thigh_imm_t = immultiply(thigh_imm_t,annulus);
    thigh_imm_t=bwmorph(thigh_imm_t,'bridge','inf');
    thigh_imm_t= imfill(thigh_imm_t,'holes');
    thigh_hand_depth = immultiply(thigh_hand_depth,annulus);
    thigh_hand_imm_t = immultiply(thigh_hand_imm_t,annulus);
    thigh_hand_imm_t=bwmorph(thigh_hand_imm_t,'bridge','inf');
    thigh_hand_imm_t= imfill(thigh_hand_imm_t,'holes');
    thigh_comp = bwconncomp(thigh_imm_t,8);
    arsr_thigh = cellfun(@numel,thigh_comp.PixelIdxList);
    [area_thigh, X_thigh] = max(arsr_thigh);
    thigh_imm_temp = false(size(thigh_depth));
    thigh_imm_temp(thigh_comp.PixelIdxList{X_thigh})=1;

    thigh_comp_T = bwconncomp(thigh_imm_T,8);
    area_thigh_T = cellfun(@numel,thigh_comp_T.PixelIdxList);
    if area_thigh < area_thigh_T*0.6
        arsr_thigh(X_thigh) = 0;
        [~, X_thigh] = max(arsr_thigh);
        thigh_imm_temp(thigh_comp.PixelIdxList{X_thigh})=1;
    end

    thigh_hand_comp = bwconncomp(thigh_hand_imm_t,8);
    arsr_thigh_hand = cellfun(@numel,thigh_hand_comp.PixelIdxList);
    [area_thigh_hand, X_thigh_hand] = max(arsr_thigh_hand);
    thigh_hand_imm_temp = false(size(thigh_hand_depth));
    thigh_hand_imm_temp(thigh_hand_comp.PixelIdxList{X_thigh_hand})=1;

    if area_thigh_hand < area_thigh_T*0.6
        arsr_thigh_hand(X_thigh_hand) = 0;
        [~, X_thigh_hand] = max(arsr_thigh_hand);
        thigh_hand_imm_temp(thigh_hand_comp.PixelIdxList{X_thigh_hand})=1;
    end

    [r_thigh_t, c_thigh_t] = find(thigh_imm_temp);
    %     h=figure;imshow(thigh_imm_temp), saveas(h,['thigh flex' num2str(lu-1) '.jpg']);
    r_thigh_t = abs(r_thigh_t-dim1_f)+1;
    [r_thigh_hand_t, c_thigh_hand_t] = find(thigh_hand_imm_temp);
    r_thigh_hand_t = abs(r_thigh_hand_t-dim1_f)+1;

    thigh_cdm_r = mean(r_thigh_t);
    thigh_cdm_c = mean(c_thigh_t);
    m_ax_thigh = (thigh_cdm_r-KN_pl(1,lu))/(thigh_cdm_c-KN_pl(2,lu));
    q_ax_thigh = KN_pl(1,lu)-m_ax_thigh*KN_pl(2,lu);
    m_ax_ort_thigh = -1/m_ax_thigh;
    [cut_low_thigh_X, cut_low_thigh_Y] = linecirc(m_ax_thigh,q_ax_thigh,KN_pl(2,lu),KN_pl(1,lu)+1,thigh_length*low_cut_par);

    [cut_low_thigh(2), temp] = max(cut_low_thigh_Y);
    cut_low_thigh(1) = cut_low_thigh_X(temp);
    q_ax_ort_thigh = cut_low_thigh(2)-m_ax_ort_thigh*(cut_low_thigh(1));

    thigh_low_control = r_thigh_t > c_thigh_t.*m_ax_ort_thigh+q_ax_ort_thigh;
    r_thigh_imm = abs(r_thigh_t(thigh_low_control)-dim1_f)+1;
    c_thigh_imm = c_thigh_t(thigh_low_control);

    thigh_imm = false(size(thigh_imm_temp));
    ind_thigh = mat2cell([r_thigh_imm c_thigh_imm], length(r_thigh_imm), ones(1, 2));
    thigh_imm(sub2ind(size(thigh_imm), ind_thigh{:})) = true;
    [r_thigh, c_thigh] = find(thigh_imm);

    list_thigh_imm{lu} = thigh_imm;

    thigh_hand_low_control = r_thigh_hand_t > c_thigh_hand_t.*m_ax_ort_thigh+q_ax_ort_thigh;
    r_thigh_hand_imm = abs(r_thigh_hand_t(thigh_hand_low_control)-dim1_f)+1;
    c_thigh_hand_imm = c_thigh_hand_t(thigh_hand_low_control);

    thigh_hand_imm = false(size(thigh_hand_imm_temp));
    ind_thigh_hand = mat2cell([r_thigh_hand_imm c_thigh_hand_imm], length(r_thigh_hand_imm), ones(1, 2));
    thigh_hand_imm(sub2ind(size(thigh_hand_imm), ind_thigh_hand{:})) = true;
    [r_thigh_hand, c_thigh_hand] = find(thigh_hand_imm);

    %% full icp
    thigh_pl = [abs(r_thigh-dim1_f)+1, c_thigh]';
    rigid_target_t = [thigh_pl(2,:); thigh_pl(1,:)];
    cdm_rigid_target_t = mean(rigid_target_t,2);
    m_ax1_rigid_target_t = (cdm_rigid_target_t(2)-KN_pl(1,lu))/(cdm_rigid_target_t(1)-KN_pl(2,lu));

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
    rigid_TR_rot = rotaz*(rigid_TR_t-cdm_rigid_source_t)+cdm_rigid_source_t;
    KN_pl_rot_T = rotaz*([KN_pl_T(2);KN_pl_T(1)]-cdm_rigid_source_t)+cdm_rigid_source_t;

    rigid_source_rot_tr = rigid_source_rot - repmat(KN_pl_rot_T,1,size(rigid_source_rot,2));
    rigid_TR_rot_tr = rigid_TR_rot - KN_pl_rot_T;
    rigid_target_tr = rigid_target_t - repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(rigid_target_t,2));

    if show_save_all
        figure;plot(rigid_target_tr(1,:),rigid_target_tr(2,:),'or')
        hold on;plot(rigid_source_rot_tr(1,:),rigid_source_rot_tr(2,:),'ob')
        hold on;plot(rigid_TR_rot_tr(1),rigid_TR_rot_tr(2),'og')
        hold on;plot(0,0,'og')
        hold on;line([0 rigid_TR_rot_tr(1)], [0 rigid_TR_rot_tr(2)], 'Color','g','Linewidth',2)
        axis equal
    end

    [~,~,template_fin_t,thigh_fin_t,res] = icp_2_rot(rigid_target_tr,rigid_source_rot_tr);

    template_fin = template_fin_t + repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(template_fin_t,2));
    TR_fin = template_fin(:,end);
    template_fin(:,end) = [];
    thigh_match{lu} = template_fin;
    thigh_fin = thigh_fin_t + repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(thigh_fin_t,2));
    thigh_data{lu} = thigh_fin;
    thigh_hand_fin = [c_thigh_hand abs(r_thigh_hand-dim1_f)+1]';

    if show_save_all
        figure;
        plot(template_fin(1,:),template_fin(2,:),'ob');
        hold on;
        plot(thigh_fin(1,:),thigh_fin(2,:),'or');
        axis equal
    end

    k_thigh = boundary(thigh_fin(1,:)',thigh_fin(2,:)',1);
    k_thigh_hand = boundary(thigh_hand_fin(1,:)',thigh_hand_fin(2,:)',1);
    k_template = boundary(template_fin(1,:)',template_fin(2,:)',1);
    x_thigh = thigh_fin(1,:)';
    y_thigh = thigh_fin(2,:)';
    x_thigh_hand = thigh_hand_fin(1,:)';
    y_thigh_hand = thigh_hand_fin(2,:)';
    x_template = template_fin(1,:)';
    y_template = template_fin(2,:)';
    bound_thigh = [x_thigh(k_thigh)';y_thigh(k_thigh)'];
    bound_thigh_hand = [x_thigh_hand(k_thigh_hand)';y_thigh_hand(k_thigh_hand)'];
    bound_template = [x_template(k_template)';y_template(k_template)'];
    BW_thigh = poly2mask(bound_thigh(1,:),abs(bound_thigh(2,:)-dim1_f)+1,dim1_f,dim2_f);
    %     figure;imshow(BW_thigh)
    BW_thigh_hand = poly2mask(bound_thigh_hand(1,:),abs(bound_thigh_hand(2,:)-dim1_f)+1,dim1_f,dim2_f);
    %     figure;imshow(BW_thigh_hand)
    BW_template = poly2mask(bound_template(1,:),abs(bound_template(2,:)-dim1_f)+1,dim1_f,dim2_f);
    %     figure;imshow(BW_template)


    new_thigh_imm = BW_thigh | (BW_thigh_hand & BW_template);
    %     figure;imshow(new_thigh_imm)
    thigh_diff = new_thigh_imm - BW_thigh;
    %     figure;imshow(thigh_diff);
    thigh_diff_pix = sum(sum(thigh_diff));
    old_thigh_diff_pix = 999999;

    while (old_thigh_diff_pix-thigh_diff_pix)>5
        old_thigh_diff_pix = thigh_diff_pix;
        [r_thigh, c_thigh] = find(new_thigh_imm);
        thigh_pl = [abs(r_thigh-dim1_f)+1, c_thigh]';

        rigid_target_t = [thigh_pl(2,:); thigh_pl(1,:)];
        cdm_rigid_target_t = mean(rigid_target_t,2);
        m_ax1_rigid_target_t = (cdm_rigid_target_t(2)-KN_pl(1,lu))/(cdm_rigid_target_t(1)-KN_pl(2,lu));

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
        rigid_TR_rot = rotaz*(rigid_TR_t-cdm_rigid_source_t)+cdm_rigid_source_t;

        KN_pl_rot_T = rotaz*([KN_pl_T(2);KN_pl_T(1)]-cdm_rigid_source_t)+cdm_rigid_source_t;

        rigid_source_rot_tr = rigid_source_rot - repmat(KN_pl_rot_T,1,size(rigid_source_rot,2));
        rigid_TR_rot_tr = rigid_TR_rot - KN_pl_rot_T;
        rigid_target_tr = rigid_target_t - repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(rigid_target_t,2));

        [~,~,template_fin_t,thigh_fin_t,res] = icp_2_rot(rigid_target_tr,rigid_source_rot_tr);

        template_fin = template_fin_t + repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(template_fin_t,2));
        TR_fin = template_fin(:,end);
        template_fin(:,end) = [];
        thigh_match{lu} = template_fin;
        thigh_fin = thigh_fin_t + repmat([KN_pl(2,lu);KN_pl(1,lu)],1,size(thigh_fin_t,2));
        thigh_data{lu} = thigh_fin;

        if show_save_all
            figure;
            plot(template_fin(1,:),template_fin(2,:),'ob');
            hold on;
            plot(thigh_fin(1,:),thigh_fin(2,:),'or');
            axis equal
        end

        k_thigh = boundary(thigh_fin(1,:)',thigh_fin(2,:)',1);
        k_thigh_hand = boundary(thigh_hand_fin(1,:)',thigh_hand_fin(2,:)',1);
        k_template = boundary(template_fin(1,:)',template_fin(2,:)',1);
        x_thigh = thigh_fin(1,:)';
        y_thigh = thigh_fin(2,:)';
        x_thigh_hand = thigh_hand_fin(1,:)';
        y_thigh_hand = thigh_hand_fin(2,:)';
        x_template = template_fin(1,:)';
        y_template = template_fin(2,:)';
        bound_thigh = [x_thigh(k_thigh)';y_thigh(k_thigh)'];
        bound_thigh_hand = [x_thigh_hand(k_thigh_hand)';y_thigh_hand(k_thigh_hand)'];
        bound_template = [x_template(k_template)';y_template(k_template)'];
        BW_thigh = poly2mask(bound_thigh(1,:),abs(bound_thigh(2,:)-dim1_f)+1,dim1_f,dim2_f);
        BW_thigh_hand = poly2mask(bound_thigh_hand(1,:),abs(bound_thigh_hand(2,:)-dim1_f)+1,dim1_f,dim2_f);
        BW_template = poly2mask(bound_template(1,:),abs(bound_template(2,:)-dim1_f)+1,dim1_f,dim2_f);
        new_thigh_imm = BW_thigh | (BW_thigh_hand & BW_template);
        %         figure;imshow(new_thigh_imm)
        thigh_diff = new_thigh_imm - BW_thigh;
        %         figure;imshow(thigh_diff);
        thigh_diff_pix = sum(sum(thigh_diff));
    end

    cdm_template_fin = mean(template_fin,2);
    m_ax1_template_fin = (cdm_template_fin(2)-KN_pl(1,lu))/(cdm_template_fin(1)-KN_pl(2,lu));

    alpha_template_fin = atan(m_ax1_template_fin);
    if alpha_template_fin<0
        alpha_template_fin = alpha_template_fin + deg2rad(180);
    end

    list_TR_pl(:,lu) = flipud(TR_fin)';
    list_TR_imm(:,lu) = [abs(TR_fin(2)-dim1_f)+1; TR_fin(1)];
    list_angles(lu) = rad2deg(atan2(TR_fin(2)-KN_pl(1,lu),TR_fin(1)-KN_pl(2,lu)));
    list_res(lu) = res;
    area = cell2mat(struct2cell(regionprops(logical(thigh_imm),'Area')));
    area1(lu) = area(1);
    centroid = cell2mat(struct2cell(regionprops(logical(thigh_imm),'Centroid')));
    centroid1(lu,:) = centroid(:,1:2);


    if show_save_all
        h_icp1 = figure;
        plot(template_fin(1,:),template_fin(2,:),'ob');
        hold on;
        plot(thigh_fin(1,:),thigh_fin(2,:),'or');
        hold on
        plot(KN_pl(2,lu),KN_pl(1,lu),'og','Linewidth',2)
        hold on
        plot(TR_fin(1),TR_fin(2),'og','Linewidth',2)
        hold on
        line([TR_fin(1) KN_pl(2,lu)], [TR_fin(2) KN_pl(1,lu)],'Color','g','Linewidth',2)
        axis equal
    end

    edge_template_fin = edge_plot(template_fin);
    edge_thigh_fin = edge_plot(thigh_fin);
    if show_save_all
        h_icp2 = figure;
        plot(edge_template_fin(1,:),edge_template_fin(2,:),'ob');
        hold on;
        plot(edge_thigh_fin(1,:),edge_thigh_fin(2,:),'or');
        hold on
        plot(KN_pl(2,lu),KN_pl(1,lu),'og','Linewidth',2)
        hold on
        plot(TR_fin(1),TR_fin(2),'og','Linewidth',2)
        hold on
        line([TR_fin(1) KN_pl(2,lu)], [TR_fin(2) KN_pl(1,lu)],'Color','g','Linewidth',2)
        axis equal
        %         saveas(h_icp2,[d_fol 'matching_thigh\' template.type '_match2_' num2str(fr)],'bmp');
    end

    if show_save_all
        h2 = figure;
        imshow(thigh_depth,[]);
        h3 = figure;
        imshow(thigh_imm);
        %         saveas(h1,[d_fol 'matching_thigh\visual\' template.type '_hist_' num2str(fr)],'bmp');
    end
    close all


end
cd (d_fol)
cd ..

results.TR_pl = list_TR_pl;
results.TR_imm = list_TR_imm;
results.thigh_pl = thigh_match;
results.thigh_imm = list_thigh_imm;
results.thigh_data_pl = thigh_data;
results.thigh_angles = list_angles;
results.thigh_res = list_res;

if show_save_all
    h_res = figure;
    plot(list_res,'b')
    hold on
    plot(list_res,'*b')
    h_ang = figure;
    plot(list_angles,'b')
    hold on
    plot(list_angles,'*b')
end
set(0, 'DefaultFigureVisible', 'off')
end