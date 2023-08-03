function subject_segmentation(d_fol,frames,B,side,~)
%Author: Diletta Balta
%Department of Electronics and Telecommunication
%Politecnico di Torino 
%diletta.balta@polito.it

%This function provides a segmentation mask of the subject for each frame
%following the methods described in the methodological paper (paragraph
%SUBJECT SEGMENTATION)

%inputs
%d_fol = folder containing the dynamic acquisition
%frames = list containing the number of frames of the gait cycle 
%B = background image
%side = R (right) or L (left)

set(0, 'DefaultFigureVisible', 'off')

for i = 1:frames
    display(['Segmenting the subject in frame ', num2str(i), '/',num2str(frames)]);
    
    if i-1<10
        I = imread([d_fol 'Video\000' num2str(i-1) '.bmp']);
        D = readmatrix([d_fol '\Misc\MATCH_d_raw\000' num2str(i-1) '.txt']);
        Depth = readmatrix([d_fol '\Misc\MATCH_d_raw\000'  num2str(i-1) '.txt']);
    elseif i-1>=10 && i-1<100
        I = imread([d_fol 'Video\00' num2str(i-1) '.bmp']);
        D = readmatrix([d_fol '\Misc\MATCH_d_raw\00' num2str(i-1) '.txt']);
        Depth = readmatrix([d_fol '\Misc\MATCH_d_raw\00'  num2str(i-1) '.txt']);
    elseif i-1>=100
        I = imread([d_fol 'Video\000' num2str(i-1) '.bmp']);
        D = readmatrix([d_fol 'Misc\MATCH_d_raw\0' num2str(i-1) '.txt']);
        Depth = readmatrix([d_fol 'Misc\MATCH_d_raw\0'  num2str(i-1) '.txt']);
    end

    Difference = abs(im2double(I)- im2double(B));
    Difference = sqrt(Difference(:,:,1).^2 + Difference(:,:,2).^2 + Difference(:,:,3).^2);
    %     figure,imshow(Difference)

    D_norm = uint8(255*(Difference-min(Difference(:)))/(max(Difference(:))-min(Difference(:))));
    [count,bins] = imhist(D_norm);
    th_start = round(sum(count.*bins)/sum(count));
    figure,imshow(D_norm)
    idx = D_norm>=th_start; %threshoding
    bw = zeros(size(D_norm));
    bw(idx) = 1;
    %     figure,imshow(bw)
    bw1 = bw;
    bw1 = imerode(bw1,strel('disk',1));
    bw1 = imdilate(bw1,strel('disk',2));

    %     figure,imshow(bw1)

    bw1(Depth>2800) = 0; %removal of residuals areas belonging to the background
    bw1(Depth==0) = 0;
    [D_imm_l_t, D_imm_r_t] = Feet_segmentation_dyn1(I); %feet segmentation
    figure,imshow(D_imm_r_t)
    D_imm_feet = D_imm_l_t | D_imm_r_t;
    cd (d_fol)
    mkdir Feet_segmentation
    h = figure;imshowpair(I,D_imm_feet),saveas(gcf,['Feet_segmentation\Feet', num2str(i-1) ,'.jpg']);
    cd ..
    % Feet segmentation refinement
    if side == 'R'
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

    %     figure,imshow(foot1_rect)
    %     figure,imshow(foot2_rect)
    [r1,~] = find(foot1_rect);
    [r2,~] = find(foot2_rect);
    minr1 = max(r1);
    minr2 = max(r2);
    bw1(max([minr1 minr2]):end,:) = 0; %removal of shadows under the feet
    %     figure,imshow(bw1)
    foot1_dil = imdilate(foot1_rect,strel('disk',2));
    foot2_dil = imdilate(foot2_rect,strel('disk',2));

    bw2 = bw1 & (not((logical(foot1_dil))|(logical(foot2_dil))));
    bw3 = bwareaopen(bw2,100);

    bw_final = bw3 | foot1_rect | foot2_rect;

    %% Green carpet segmentation (optional) 
    %Based on the environment/lights conditions, the thresholds for the color
    %filter have to be adjusted

    Carpet = zeros(size(bw_final));
    for r = 510:670
        for c = 1:size(I,2)
            pix = [I(r,c,1) I(r,c,2) I(r,c,3)];
            if  pix(1)<50 && pix(2)>40 && pix(3)<70 && foot1_dil(r,c) == 0 && foot2_dil(r,c) == 0
                Carpet(r,c) = true;
            end

        end
    end
    figure('Visible','off'),
    imshow(Carpet)
    bw_final = bw_final - Carpet;
    bw_final(bw_final == -1) = 0;
    figure('Visible','off'),imshow(bw_final)
    if i-1<10
        imwrite(bw_final,[d_fol 'Misc\Segm\000' num2str(i-1) '.bmp'])
    elseif i-1>=10 && i-1<100
        imwrite(bw_final,[d_fol 'Misc\Segm\00' num2str(i-1) '.bmp'])
    elseif i-1>=100
        imwrite(bw_final,[d_fol 'Misc\Segm\0' num2str(i-1) '.bmp'])
    end
    close all
end
end
