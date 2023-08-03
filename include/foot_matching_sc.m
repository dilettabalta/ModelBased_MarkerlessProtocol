function results = foot_matching_sc(d_fol,side_type,lis,template,~)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides LM positions in each frame by matching the foot template to each foot segmentation mask by using ICP following the methods described in the methodological paper
% paragraph JOINT CENTERS TRAJECTORIES ESTIMATION - Ankle joint center (AJC) estimation).

%inputs
%d_fol = folder containing the dynamic acquisition
%lis = list containing the number of frames of the gait cycle
%template = static/loading/swing foot template

%outputs
%results = a MATLAB structure containing LM position on each frame

close all;
set(0, 'DefaultFigureVisible', 'off')
inizio = 1;
dim1_f = 720;
dim2_f = 1280;
show_save_all = true;
thr = 0.90;
thr_sole_start = 0.2;
thr_sole_end = 0.5; %threshold for the slope of the foot
fl_sc = 0;

%% template
template_rect = template.foot;
[e_r_T,e_c_T] = find(template_rect);
e_r_pl_T = abs(e_r_T-dim1_f)+1;
e_c_pl_T = e_c_T;
data_pl_T = [e_r_pl_T';e_c_pl_T'];
cdm_pl_T = mean(data_pl_T,2);

if show_save_all
    mkdir([d_fol 'matching_foot\']);
end

%creating variables
list_ML_pl = zeros(2,length(lis));
list_ML_imm = zeros(2,length(lis));

list_alpha_Reallignedrect = zeros(1,length(lis));
list_angles_D = zeros(1,length(lis));
scaling_icp = zeros(1,length(lis));
scaling_edge_cut  = zeros(1,length(lis));
clearvars list_foot_pl

for lu = inizio:length(lis)
    drawnow
    fr = lis(lu);
    disp(fr)
    if fr<10
        rgb_or = imread([d_fol 'Video\000' num2str(fr) '.bmp']);
    elseif fr>=10 && fr<100
        rgb_or = imread([d_fol 'Video\00' num2str(fr) '.bmp']);
    elseif fr>=100
        rgb_or = imread([d_fol 'Video\0' num2str(fr) '.bmp']);
    end

    if fr<10
        segm_rect = logical(imread([d_fol 'Misc\Segm\000' num2str(fr) '.bmp']));
    elseif fr>=10 && fr<100
        segm_rect = logical(imread([d_fol 'Misc\Segm\00' num2str(fr) '.bmp']));
    elseif fr>=100
        segm_rect = logical(imread([d_fol 'Misc\Segm\0' num2str(fr) '.bmp']));
    end

    %removing useless data for the feet
    segm_rect_t = segm_rect;
    segm_rect_t(1:template.knee(1),:) = false;

    %%feet segmentation
    I = rgb_or;
    [D_imm_l_t, D_imm_r_t] = Feet_segmentation_dyn1(I);

    %figure,imshow(D_imm_r_t|D_imm_l_t)
    %h = figure; imshow(I), hold on,visboundaries(D_imm_r_t|D_imm_l_t)

    % Feet refinement through morphological operators
    if side_type == 'L'
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
    else
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
        results.foot1_imm{lu} = foot1_rect;
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
    end

    D1_imm_t = bwmorph(D1_imm_t,'bridge','inf');
    D1_imm_t = imfill(D1_imm_t,'holes');
    D1_imm_t = bwmorph(D1_imm_t,'majority','inf');

    D1_el = bwconncomp(D1_imm_t,8);
    arsr_D1 = cellfun(@numel,D1_el.PixelIdxList);
    [~, X_D1] = max(arsr_D1);

    D1_imm = false(dim1_f,dim2_f);
    D1_imm(D1_el.PixelIdxList{X_D1}) = 1;
    foot1_rect = D1_imm;
    results.foot1_imm{lu} = foot1_rect;

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
    results.foot2_imm{lu} = foot2_rect;

    %     figure,imshow(foot1_rect)

    %%contour
    [e_r_in_D,e_c_in_D] = find(edge(foot1_rect,'canny'));

    %inverting and scaling y-axis
    e_rr_in_D = abs(e_r_in_D-dim1_f)+1;
    e_cc_in_D = e_c_in_D;
    e_rr_D = e_rr_in_D;
    e_cc_D = e_cc_in_D;

    %%area
    [f_r_in_D,f_c_in_D] = find(D1_imm);
    %     figure,imshow(D1_imm)

    %inverting and scaling y-axis
    f_rr_in_D = abs(f_r_in_D-dim1_f)+1;
    f_cc_in_D = f_c_in_D;
    f_rr_D = f_rr_in_D;
    f_cc_D = f_cc_in_D;
    e_rr_in_D_temp = e_rr_in_D';
    e_cc_in_D_temp = e_cc_in_D';

    %finding the sole of the foot
    clearvars rr_in_low_D cc_in_low_D

    if lu == inizio
        for i=min(e_cc_in_D_temp):max(e_cc_in_D_temp)
            cand_D = e_rr_in_D_temp(e_cc_in_D_temp==i);
            cand_cut_D = find(diff(cand_D)<-ceil(0.75*length(cand_D)),1,'last');
            if isempty(cand_cut_D)
                ok_cand_D = cand_D;
            else
                ok_cand_D = cand_D(cand_cut_D+1:end);
            end

            if exist('rr_in_low_D','var')
                rr_in_low_D = [rr_in_low_D ok_cand_D];
            else
                rr_in_low_D = ok_cand_D;
            end

            if exist('cc_in_low_D','var')
                cc_in_low_D = [cc_in_low_D ones(1,length(ok_cand_D))*i];
            else
                cc_in_low_D = ones(1,length(ok_cand_D))*i;
            end
        end

        if fl_sc == 0
            cc_low_D = cc_in_low_D;
            rr_low_D = rr_in_low_D;
        else
            DD_f_in = [f_rr_in_D f_cc_in_D];
            DD_f_in_low = [rr_in_low_D' cc_in_low_D'];
            [~, ia_D, ~] = intersect(DD_f_in,DD_f_in_low,'rows');
            DD_f = [f_rr_D f_cc_D];
            rr_low_D = DD_f(ia_D,1)';
            cc_low_D = DD_f(ia_D,2)';
        end

        if side_type == 'L'
            cc_low_D = fliplr(cc_low_D);
            rr_low_D = fliplr(rr_low_D);
        end
        %% cutting the dynamic foot
        %computing the slope of the sole as the line best fitting the point between
        %heel and metatarsus

        ind_start_sole_D = round(thr_sole_start*length(cc_low_D));
        ind_end_sole_D = round(thr_sole_end*length(cc_low_D));
        rr_sole_D = rr_low_D(ind_start_sole_D:ind_end_sole_D);
        cc_sole_D = cc_low_D(ind_start_sole_D:ind_end_sole_D);
        fun_D = fit(cc_sole_D',rr_sole_D','poly1');
        coeff_temp_D = coeffvalues(fun_D);
        m_sole_D = coeff_temp_D(1);
        m_rect_D = -1/m_sole_D;
        q_sole_D = coeff_temp_D(2);

        alpha_D = atan(m_sole_D);
    elseif lu == inizio+1
        alpha_D = list_alpha_Reallignedrect(lu-1);
    elseif lu == inizio+2
        alpha_D = 2*list_alpha_Reallignedrect(lu-1)-list_alpha_Reallignedrect(lu-2);
    else
        previous = list_alpha_Reallignedrect(lu-3:lu-1);
        alpha_D = 2.5*previous(3)-2*previous(2)+0.5*previous(1);
    end

    %% cutting the dynamic foot
    %centering data for an easier identification of the sole
    R0_1_D = [cos(alpha_D) -sin(alpha_D);sin(alpha_D) cos(alpha_D)];
    data_pl_D = [e_rr_D';e_cc_D'];
    O1_1_D = -R0_1_D*[mean(data_pl_D(1,:));mean(data_pl_D(2,:))];
    T0_1_D = zeros(3);
    T0_1_D(1:2,1:2) = R0_1_D;
    T0_1_D(:,3) = [O1_1_D;1];
    data_1_D = T0_1_D*[data_pl_D;ones(1,length(e_rr_D))];
    data_1_D = data_1_D(1:2,:);

    check_heel_D = 0;
    heel_D = 0;
    cut_length = template.cut;
    switch side_type
        case 'R'
            while check_heel_D<=thr*length(e_cc_D)
                heel_D=heel_D-1;
                check_heel_D = sum(data_1_D(2,:)>heel_D);
            end
            cut_D = heel_D + cut_length;
            ind_to_cut_D = find(data_1_D(2,:)>cut_D);
        case 'L'
            while check_heel_D<=thr*length(e_cc_D)
                heel_D=heel_D+1;
                check_heel_D = sum(data_1_D(2,:)<heel_D);
            end
            cut_D = heel_D - cut_length;
            ind_to_cut_D = find(data_1_D(2,:)<cut_D);
    end

    data_pl_D(:,ind_to_cut_D) = [];
    e_c_pl_D = data_pl_D(2,:);
    e_r_pl_D = data_pl_D(1,:);
    cdm_pl_D = mean(data_pl_D,2);

    foot_perim_cut = data_pl_D;
    foot_perim_cut(1,:) = abs(foot_perim_cut(1,:)-dim1_f)+1;
    foot_cut = false(dim1_f,dim2_f);
    ind_foot_dist = mat2cell(foot_perim_cut', size(foot_perim_cut', 1), ones(1, size(foot_perim_cut', 2)));
    foot_cut(sub2ind(size(foot_cut), ind_foot_dist{:})) = true;

    %% ICP
    %centering data to an easier identification of the sole of the template
    m_sole_T = template.m_sole;
    alpha_T = atan(m_sole_T);
    q_sole_T = template.q_sole;
    alpha = alpha_D - alpha_T;
    rotaz = [cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
    data_pl_T_rot = rotaz*(data_pl_T-repmat(cdm_pl_T,1,size(data_pl_T,2)))+repmat(cdm_pl_T,1,size(data_pl_T,2));
    ML_in = template.mal';
    ML_pl = ML_in;
    ML_pl(1) = abs(ML_pl(1)-dim1_f)+1;
    ML_pl_rot = rotaz*(ML_pl-cdm_pl_T)+cdm_pl_T;
    x_rect_sole_T = min(e_c_pl_T)-5:max(e_c_pl_T)+5;
    rect_sole_T = [m_sole_T*(x_rect_sole_T)+q_sole_T;x_rect_sole_T];
    x_rect_sole_D = min(e_c_pl_D)-5:max(e_c_pl_D)+5;
    rect_sole_D = [m_sole_D*(x_rect_sole_D)+q_sole_D;x_rect_sole_D];

    % pre-rototranslation of the template according alpha and trasl
    rect_sole_T_rot = rotaz*(rect_sole_T-repmat(cdm_pl_T,1,size(rect_sole_T,2)))+repmat(cdm_pl_T,1,size(rect_sole_T,2));
    trasl = cdm_pl_T-cdm_pl_D;

    clearvars data_pl_tr_T_rot rect_sole_tr_T_rot
    data_pl_tr_T_rot(1,:) = data_pl_T_rot(1,:) - trasl(1);
    data_pl_tr_T_rot(2,:) = data_pl_T_rot(2,:) - trasl(2);
    ML_pl_tr_rot = ML_pl_rot - trasl;

    rect_sole_tr_T_rot(1,:) = rect_sole_T_rot(1,:) - trasl(1);
    rect_sole_tr_T_rot(2,:) = rect_sole_T_rot(2,:) - trasl(2);

    rigid_target = [data_pl_D(2,:)' data_pl_D(1,:)'];
    rigid_source = [data_pl_tr_T_rot(2,:)' data_pl_tr_T_rot(1,:)'];

    %ICP
    [~,Reallignedsource,transform] = rigidICP(rigid_target,rigid_source);
    list_foot_pl{lu} = [Reallignedsource(:,2),Reallignedsource(:,1)]';

    Reallignedrect = transform.b*[rect_sole_tr_T_rot(2,:)' rect_sole_tr_T_rot(1,:)']*transform.T + transform.c(1:size(rect_sole_tr_T_rot,2),:);
    fun_Reallignedrect = fit(Reallignedrect(:,1),Reallignedrect(:,2),'poly1');
    coeff_temp_Reallignedrect = coeffvalues(fun_Reallignedrect);

    ReallignedML = transform.b * [ML_pl_tr_rot(2) ML_pl_tr_rot(1)] * transform.T + transform.c(1,:);
    list_ML_pl(:,lu) = fliplr(ReallignedML)';
    list_ML_imm(:,lu) = [abs(ReallignedML(2)-dim1_f)+1 ReallignedML(1)]';

    list_alpha_Reallignedrect(lu) = atan(coeff_temp_Reallignedrect(1));
    if lu~=inizio && abs(list_alpha_Reallignedrect(lu)-list_alpha_Reallignedrect(lu-1))>pi/2
        switch side_type
            case 'R'
                list_alpha_Reallignedrect(lu) = list_alpha_Reallignedrect(lu)-pi;
            case 'L'
                list_alpha_Reallignedrect(lu) = list_alpha_Reallignedrect(lu)+pi;
        end
    end

    list_angles_D(lu) = rad2deg(alpha_D);
    scaling_icp(lu) = transform.b;
    scaling_edge_cut(lu) = size(data_pl_D,2)/size(data_pl_T,2);

    if show_save_all
        h2 = figure;
        plot(Reallignedsource(:,1),Reallignedsource(:,2),'og');
        hold on;
        plot(rigid_target(:,1),rigid_target(:,2),'or');
        hold on
        plot(ReallignedML(1),ReallignedML(2),'ob','Linewidth',2)
        axis equal
        %         saveas(h2,[d_fol 'matching_foot\' template.type '_match_' num2str(fr)],'bmp');
    end
end

results.ML_pl = list_ML_pl;
results.ML_imm = list_ML_imm;
results.foot_pl = list_foot_pl;

if side_type=='L'
    list_alpha_Reallignedrect = list_alpha_Reallignedrect + deg2rad(180);
    list_angles_D = list_angles_D + 180;
    alpha_T = alpha_T + deg2rad(180);
end

list_angles_Reallignedrect = rad2deg(list_alpha_Reallignedrect);
results.foot_alpha = list_alpha_Reallignedrect;
results.foot_angles = list_angles_Reallignedrect;
results.scaling = scaling_icp;
close all
end

