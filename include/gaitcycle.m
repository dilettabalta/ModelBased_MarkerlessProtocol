function [start1, stop1, start2, pos_IC1, pos_IC2, pos_IC1_c, stride_length, step_length, gait_speed] = gaitcycle(d_fol, start, stop, side, conv_factor)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino 
%diletta.balta@polito.it

%This function identifies the first and the last frame of the most central
%gait cycle along with the frame representing the contact of the background foot.
%This function provides the stride length, step length and gait speed
%following the methods described in the methodological paper (paragraph
%GAIT CYCLE IDENTIFICATION)

%inputs
%d_fol = folder containing the dynamic acquisition
%start = first frame of the dynamic acquisition where both feet must be visible and
%on the green carpet
%stop = last frame of the dynamic acquisition where both feet must be visible and
%on the green carpet
%side = R (right) or L (left)
%conv_factor = conversion factor (m/pixel) (see Appendix A of the paper for
%more details)

%outputs
%start1 = first frame of the gait cycle
%stop1 = last frame of the gait cycle
%start2 = frame representing the contact of the background foot 
%pos_IC1 = position of first initial contact of the foreground foot
%pos_IC2 = position of the following initial contact of the foreground foot
%pos_IC1_c = position of the first initial contact of the backgroung foot
%stride_length = stride length
%step_length = step length
%gait_speed = gait speed


v = start:stop;

for fr = 1:length(v)

    I = imread([d_fol '\color_stream\color' num2str(v(fr)) '.png']);
    %     figure,imshow(I)

    %% Feet segmentation through two color filters
    [D_imm_l_t, D_imm_r_t,~] = Feet_segmentation_dyn1 (I);

    %     h = figure; imshow(I), hold on,visboundaries(D_imm_r_t|D_imm_l_t)

    %Feet segmentation refinement through morphological operators

    if side == 'L'
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
    %     h = figure; imshow(I), hold on,visboundaries(foot1_rect|foot2_rect)

    %     if ~exist([d_fol '\FeetSegmentation'])
    %         mkdir([d_fol '\FeetSegmentation'])
    %     end

    %     saveas(h,[d_fol '\FeetSegmentation' filesep 'Piedi' num2str(v(fr)) '.bmp']);
    colonne = sum(foot1_rect);
    righe = sum(foot1_rect,2);
    col = find(colonne);
    rig = find(righe);

    if col(1)-30<1
        x = 30-(abs(col(1)-30))-1;
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-x:col(end)+30);
    elseif  col(end)+30>size(I,2)
        x = abs(size(I,2)- col(end));
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1):col(end)+x);
    elseif col(1)-30>=1|| col(end)+30<=size(I,2)
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+30);
    end

    %% foregound foot: Foot inclination wrt the ground
    alpha = regionprops(logical(foot1_rect),'Orientation');
    alpha = cell2mat(struct2cell(alpha));
    alpha = alpha(1);
    ang(fr) = alpha;

    %% background foot: Foot inclination wrt the ground
    foot2_rect = bwareaopen(foot2_rect,50);
    areas = regionprops(foot2_rect);
    if size(areas,1)==1
        foot1_rect = [];
        foot1_rect = foot2_rect;
        colonne = sum(foot1_rect);
        righe = sum(foot1_rect,2);
        col = find(colonne);
        rig = find(righe);

        if col(1)-30<1
            x = 30-(abs(col(1)-30))-1;
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-x:col(end)+30);
        elseif  col(end)+30>size(I,2)
            x = abs(size(I,2)- col(end));
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+x);
        elseif col(1)-30>=1 || col(end)+30<=size(I,2)
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+30);
        end

        alpha = regionprops(logical(foot1_rect),'Orientation');
        alpha = cell2mat(struct2cell(alpha));
        alpha = alpha(1);
        ang1(fr)= alpha;

    end
end

ang1 = fillmissing(ang1,'previous');
ang = fillmissing(ang,'previous');
% figure,plot(ang)
% figure,plot(ang1)
diffang1 = diff(ang1);
diffang = diff(ang);

ind = find(abs(diffang)>120);
ind1 = find(abs(diffang1)>120);
vetang = zeros(1,stop-start+1);
vetang1 = zeros(1,stop-start+1);

if ~isempty (ind1)
    for i=1:2:length(ind1)-1
        vetang1(ind1(i)+1:ind1(i+1)) = 1;
    end
end

if ~isempty (ind)
    for i=1:2:length(ind)-1
        vetang(ind(i)+1:ind(i+1)) = 1;
    end
end


for fr = 1:length(v)

    I = imread([d_fol '\color_stream\color' num2str(v(fr)) '.png']);
    %     figure,imshow(I)
    [D_imm_l_t, D_imm_r_t,~] = Feet_segmentation_dyn1 (I);
    %     h = figure; imshow(I), hold on,visboundaries(D_imm_r_t|D_imm_l_t)

    %Feet segmentation refinement through morphological operators

    if side == 'L'
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

    %% Foreground foot
    colonne = sum(foot1_rect);
    righe = sum(foot1_rect,2);
    col = find(colonne);
    rig = find(righe);

    if col(1)-30<1
        x = 30-(abs(col(1)-30))-1;
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-x:col(end)+30);
    elseif  col(end)+30>size(I,2)
        x = abs(size(I,2)- col(end));
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1):col(end)+x);
    elseif col(1)-30>=1|| col(end)+30<=size(I,2)
        foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+30);
    end

    alpha = regionprops(logical(foot1_rect),'Orientation');
    alpha = cell2mat(struct2cell(alpha));
    alpha = alpha(1);
    ang(fr) = alpha;

    if vetang(fr) ==1
        I2 = imrotate(double(foot1_rect),180-(alpha),'bilinear','crop'); %realignment with respect to the horizontal axis
        I2_toe = imrotate(double(foot1_rect),180-(alpha),'bilinear','crop'); %realignment with respect to the horizontal axis
    else
        I2 = imrotate(double(foot1_rect),-(alpha),'bilinear','crop'); %realignment with respect to the horizontal axis
        I2_toe = imrotate(double(foot1_rect),-(alpha),'bilinear','crop'); %realignment with respect to the horizontal axiso per metterlo orizzontale
    end

    Extrem = regionprops(I2_toe,'Extrema');

    h=1;
    for i=1:size(Extrem.Extrema,1)-1
        for j=i+1:size(Extrem.Extrema,1)
            distance(h) = sqrt((Extrem.Extrema(i,1)-Extrem.Extrema(j,1))^2 + (Extrem.Extrema(i,2)-Extrem.Extrema(j,2))^2); %#ok<AGROW>
            h=h+1;
        end
    end

    h = 0;
    [~,idx] = max(distance);
    if idx<=7
        r1 = 1;
        r2 = idx+r1;
    elseif idx>7 && idx<=13
        r1 = 2;
        r2 = (idx-7)+r1;
    elseif idx>13 && idx<=18
        r1 = 3;
        r2 = (idx-13)+r1;
    elseif idx>18 && idx<=22
        r1 = 4;
        r2 = (idx-18)+r1;
    elseif idx>22 && idx<=25
        r1 = 5;
        r2 = (idx-22)+r1;
    elseif idx>25&& idx<=27
        r1 = 6;
        r2 = (idx-25)+r1;
    elseif idx>27
        r1 = 7;
        r2 = (idx-27)+r1;
    end

    if side =='L'
        if Extrem.Extrema(r1,1)>Extrem.Extrema(r2,1)
            heel = Extrem.Extrema(r1,1);
            toe = Extrem.Extrema(r2,1);
            %             f = figure;
            %             imshow(I2_toe)
            %             hold on
            %             plot(toe,Extrem.Extrema(r2,2),'.r','Linewidth',3)
            %             hold on
            %             plot(heel,Extrem.Extrema(r1,2),'.r','Linewidth',3)
            I2_toe(round(Extrem.Extrema(r2,2)),round(toe)) = 100;
        else
            heel = Extrem.Extrema(r2,1);
            toe = Extrem.Extrema(r1,1);
            %             f = figure;
            %             imshow(I2_toe)
            %             hold on
            %             plot(toe,Extrem.Extrema(r1,2),'.r','Linewidth',3)
            %             hold on
            %             plot(heel,Extrem.Extrema(r2,2),'.r','Linewidth',3)
            I2_toe(round(Extrem.Extrema(r1,2)),round(toe)) = 100;
        end
    else
        if Extrem.Extrema(r1,1)>Extrem.Extrema(r2,1)
            heel = Extrem.Extrema(r2,1);
            toe = Extrem.Extrema(r1,1);
            %             f = figure;
            %             imshow(I2_toe)
            %             hold on
            %             plot(toe,Extrem.Extrema(r2,2),'.r','Linewidth',3)
            %             hold on
            %             plot(heel,Extrem.Extrema(r1,2),'.r','Linewidth',3)
            I2_toe(round(Extrem.Extrema(r1,2)),round(toe)) = 100;
        else
            heel = Extrem.Extrema(r1,1);
            toe = Extrem.Extrema(r2,1);
            %             f = figure;
            %             imshow(I2_toe)
            %             hold on
            %             plot(toe,Extrem.Extrema(r1,2),'.r','Linewidth',3)
            %             hold on
            %             plot(heel,Extrem.Extrema(r2,2),'.r','Linewidth',3)
            I2_toe(round(Extrem.Extrema(r2,2)),round(toe)) = 100;
        end
    end

    d = edge(I2_toe,'canny');
    [r,c] = find(d);
    c1 = c;
    xp = zeros(1,length(c1));
    yp = zeros(1,length(c1));

    for i = 1:length(c1)
        supp = d(:,c1(i));
        ind = find(supp);
        yp(i) = ind(end);
        xp(i) = c1(i);
        ind = [];
    end

    %Identification of Q point (see the methodological paper for more
    %details)
    yp1 = yp(yp>=max(yp)-25);
    xp1 = xp(yp>=max(yp)-25);
    coef = polyfit(xp1,yp1,1);
    y = polyval(coef,xp1(1)-60:xp1(end)+60);

    %     f = figure;
    %     imshow(d)
    %     hold on
    %     plot(xp1(1)-60:xp1(end)+60,y)

    r1 = r;
    xp2 = zeros(1,length(r1));
    yp2 = zeros(1,length(r1));

    for i = 1:length(r1)
        supp = d(r1(i),:);
        ind = find(supp);
        if side =='L'
            xp2(i) = ind(end);
            yp2(i) = r1(i);
            ind = [];
        else
            xp2(i) = ind(1);
            yp2(i) = r1(i);
            ind = [];
        end
    end

    %Identification of R point (see the methodological paper for more
    %details)

    if side=='L'
        yp3 = yp2(xp2>=max(xp2)-25);
        xp3 = xp2(xp2>=max(xp2)-25);
        coef1 = polyfit(yp3,xp3,1);
        y1 = polyval(coef1,yp3(1)-60:yp3(end)+60);
        [hx, hy] = polyxpoly(xp1(1)-60:xp1(end)+60,y,y1,yp3(1)-60:yp3(end)+60);
    else
        yp3 = yp2(xp2<=min(xp2)+1);
        xp3 = xp2(xp2<=min(xp2)+1);
        coef1 = polyfit(yp3,xp3,1);
        y1 = polyval(coef1,yp3(1)-60:yp3(end)+60);
        [hx, hy] = polyxpoly(xp1(1)-60:xp1(end)+60,y,y1,yp3(1)-60:yp3(end)+60); %MRF = intersaction between the sole of the foot and the posterior line of the foot
    end

    %     hold on
    %     plot(y1,yp3(1)-60:yp3(end)+60)
    %     hold on
    %     plot(hx,hy,'.','Linewidth',3)
    I2(fix(hy),fix(hx)) = 20;

    if vetang(fr) ==1
        I3 = imrotate(double(I2),-(180-(alpha)),'bilinear','crop'); %foot to the starting position
        I3_toe = imrotate(double(I2_toe),-(180-(alpha)),'bilinear','crop'); %foot to the starting position
    else
        I3 = imrotate(double(I2),(alpha),'bilinear','crop'); %foot to the starting position
        I3_toe = imrotate(double(I2_toe),(alpha),'bilinear','crop'); %foot to the starting position
    end

    [heelrow ,heelcol] = find(I3==max(max(I3)));
    heelrow = heelrow+rig(1)-30;

    [toerow ,toecol] = find(I3_toe==max(max(I3_toe)));
    toerow = toerow+rig(1)-30;

    if col(1)-30<1
        heelcol = heelcol+col(1)-x;
        toecol = toecol+col(1)-x;
    elseif  col(end)+30>size(I,2)
        heelcol = heelcol+col(1)-30;
        toecol = toecol+col(1)-30;
    elseif col(1)-30>=1 || col(end)+30<=size(I,2)
        heelcol = heelcol+col(1)-30;
        toecol = toecol+col(1)-30;
    end

    %Heel (MRF) and Toe position (FF) of the foreground foot for each frame
    heelr(fr) = heelrow(1);
    heelc(fr) = heelcol(1);
    toer(fr) = toerow(1);
    toec(fr) = toecol(1);

    %% Background foot
    foot2_rect = bwareaopen(foot2_rect,50);
    areas = regionprops(foot2_rect);
    if size(areas,1)==1
        foot1_rect = [];
        foot1_rect = foot2_rect;
        colonne = sum(foot1_rect);
        righe = sum(foot1_rect,2);
        col = find(colonne);
        rig = find(righe);
        if col(1)-30<1
            x = 30-(abs(col(1)-30))-1;
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-x:col(end)+30);
        elseif  col(end)+30>size(I,2)
            x = abs(size(I,2)- col(end));
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+x);
        elseif col(1)-30>=1 || col(end)+30<=size(I,2)
            foot1_rect = foot1_rect(rig(1)-30:rig(end)+30,col(1)-30:col(end)+30);
        end

        alpha = regionprops(logical(foot1_rect),'Orientation');
        alpha = cell2mat(struct2cell(alpha));
        alpha = alpha(1);
        ang1(fr)= alpha;

        if vetang1(fr) ==1
            I2 = imrotate(double(foot1_rect),180-(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
            I2_toe = imrotate(double(foot1_rect),180-(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
        else
            I2 = imrotate(double(foot1_rect),-(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
            I2_toe = imrotate(double(foot1_rect),-(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
        end
        Extrem = regionprops(I2_toe,'Extrema');
        h=1;
        for i=1:size(Extrem.Extrema,1)-1
            for j=i+1:size(Extrem.Extrema,1)
                distance(h) = sqrt((Extrem.Extrema(i,1)-Extrem.Extrema(j,1))^2 + (Extrem.Extrema(i,2)-Extrem.Extrema(j,2))^2);
                h=h+1;
            end
        end

        h = 0;
        [~,idx] = max(distance);
        if idx<=7
            r1 = 1;
            r2 = idx+r1;
        elseif idx>7 && idx<=13
            r1 = 2;
            r2 = (idx-7)+r1;
        elseif idx>13 && idx<=18
            r1 = 3;
            r2 = (idx-13)+r1;
        elseif idx>18 && idx<=22
            r1 = 4;
            r2 = (idx-18)+r1;
        elseif idx>22 && idx<=25
            r1 = 5;
            r2 = (idx-22)+r1;
        elseif idx>25&& idx<=27
            r1 = 6;
            r2 = (idx-25)+r1;
        elseif idx>27
            r1 = 7;
            r2 = (idx-27)+r1;
        end

        if side =='L'
            if Extrem.Extrema(r1,1)>Extrem.Extrema(r2,1)
                heel = Extrem.Extrema(r1,1);
                toe = Extrem.Extrema(r2,1);
                %                 f = figure;
                %                 imshow(I2_toe)
                %                 hold on
                %                 plot(toe,Extrem.Extrema(r2,2),'.r','Linewidth',3)
                %                 hold on
                %                 plot(heel,Extrem.Extrema(r1,2),'.r','Linewidth',3)
                I2_toe(round(Extrem.Extrema(r2,2)),round(toe)) = 100;
            else
                heel = Extrem.Extrema(r2,1);
                toe = Extrem.Extrema(r1,1);
                %                 f = figure;
                %                 imshow(I2_toe)
                %                 hold on
                %                 plot(toe,Extrem.Extrema(r1,2),'.r','Linewidth',3)
                %                 hold on
                %                 plot(heel,Extrem.Extrema(r2,2),'.r','Linewidth',3)
                I2_toe(round(Extrem.Extrema(r1,2)),round(toe)) = 100;
            end
        else
            if Extrem.Extrema(r1,1)>Extrem.Extrema(r2,1)
                heel = Extrem.Extrema(r2,1);
                toe = Extrem.Extrema(r1,1);
                %                 f=figure;
                %                 imshow(I2_toe)
                %                 hold on
                %                 plot(toe,Extrem.Extrema(r1,2),'.r','Linewidth',3)
                %                 hold on
                %                 plot(heel,Extrem.Extrema(r2,2),'.r','Linewidth',3)
                I2_toe(round(Extrem.Extrema(r1,2)),round(toe)) = 100;
            else
                heel = Extrem.Extrema(r1,1);
                toe = Extrem.Extrema(r2,1);
                %                 f = figure;
                %                 imshow(I2_toe)
                %                 hold on
                %                 plot(toe,Extrem.Extrema(r2,2),'.r','Linewidth',3)
                %                 hold on
                %                 plot(heel,Extrem.Extrema(r1,2),'.r','Linewidth',3)
                I2_toe(round(Extrem.Extrema(r2,2)),round(toe)) = 100;
            end
        end

        d = edge(I2_toe,'canny');
        [r,c] = find(d);
        c1 = c;
        xp = zeros(1,length(c1));
        yp = zeros(1,length(c1));
        for i = 1:length(c1)
            supp = d(:,c1(i));
            ind = find(supp);
            yp(i) = ind(end);
            xp(i) = c1(i);
            ind = [];
        end

        % Identification of Q point
        yp1 = yp(yp>=max(yp)-25);
        xp1 = xp(yp>=max(yp)-25);
        coef = polyfit(xp1,yp1,1);
        y = polyval(coef,xp1(1)-60:xp1(end)+60);

        %         f = figure;
        %         imshow(d)
        %         hold on
        %         plot(xp1(1)-60:xp1(end)+60,y)
        r1 = r;
        xp2 = zeros(1,length(r1));
        yp2 = zeros(1,length(r1));
        for i = 1:length(r1)
            supp = d(r1(i),:);
            ind = find(supp);
            if side =='L'
                xp2(i) = ind(end);
                yp2(i) = r1(i);
                ind = [];
            else
                xp2(i) = ind(1);
                yp2(i) = r1(i);
                ind = [];
            end
        end

        %Identification of R point
        if side=='L'
            yp3 = yp2(xp2>=max(xp2)-25);
            xp3 = xp2(xp2>=max(xp2)-25);
            coef1 = polyfit(yp3,xp3,1);
            y1 = polyval(coef1,yp3(1)-60:yp3(end)+60);
            [hx, hy] = polyxpoly(xp1(1)-60:xp1(end)+60,y,y1,yp3(1)-60:yp3(end)+60);
        else
            yp3 = yp2(xp2<=min(xp2)+1);
            xp3 = xp2(xp2<=min(xp2)+1);
            coef1 = polyfit(yp3,xp3,1);
            y1 = polyval(coef1,yp3(1)-60:yp3(end)+60);
            [hx, hy] = polyxpoly(xp1(1)-60:xp1(end)+60,y,y1,yp3(1)-60:yp3(end)+60);
        end

        %         hold on
        %         plot(y1,yp3(1)-60:yp3(end)+60)
        %         hold on
        %         plot(hx,hy,'.','Linewidth',3)
        I2(fix(hy),fix(hx)) = 20;

        if vetang1(fr) ==1
            I3 = imrotate(double(I2),-(180-(alpha)),'bilinear','crop'); %ruoto per metterlo orizzontale
            I3_toe = imrotate(double(I2_toe),-(180-(alpha)),'bilinear','crop'); %ruoto per metterlo orizzontale
        else
            I3 = imrotate(double(I2),(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
            I3_toe = imrotate(double(I2_toe),(alpha),'bilinear','crop'); %ruoto per metterlo orizzontale
        end

        [heelrow ,heelcol] = find(I3==max(max(I3)));
        heelrow = heelrow+rig(1)-30;
        [toerow ,toecol] = find(I3_toe==max(max(I3_toe)));
        toerow = toerow+rig(1)-30;

        if col(1)-30<1
            heelcol = heelcol+col(1)-x;
            toecol = toecol+col(1)-x;
        elseif  col(end)+30>size(I,2)
            heelcol = heelcol+col(1)-30;
            toecol = toecol+col(1)-30;
        elseif col(1)-30>=1 || col(end)+30<=size(I,2)
            heelcol = heelcol+col(1)-30;
            toecol = toecol+col(1)-30;
        end

        %Heel (MRF) and Toe position (FF) of the background foot for each frame
        heelr_c(fr) = heelrow(1);
        heelc_c(fr) = heelcol(1);
        toer_c(fr) = toerow(1);
        toec_c(fr) = toecol(1);

        %         f = figure;
        %         imshow(I)
        %         hold on
        %         plot(heelc(fr),heelr(fr),'.r','Linewidth',5)
        %         plot(toec(fr),toer(fr),'.r','Linewidth',5)
        %         hold on
        %         plot(heelc_c(fr),heelr_c(fr),'.r','Linewidth',5)
        %         plot(toec_c(fr),toer_c(fr),'.r','Linewidth',5)
        %         saveas(f,[d_fol '\FeetSegmentation' filesep 'Piedi' num2str(v(fr)) '.bmp']);
    else
        %         f = figure;
        %         imshow(I)
        %         hold on
        %         plot(heelcol,heelrow,'.r','Linewidth',5)
        %         plot(toecol,toerow,'.r','Linewidth',5)
        %         saveas(f,[d_fol '\FeetSegmentation' filesep 'Piedi' num2str(v(fr)) '.bmp']);
    end
    close all
end

%% Foregroung foot: Identification of IC#1 and IC#2

set(0, 'DefaultFigureVisible', 'on')
% figure,plot(heelc,'.-')
% hold on
% plot(heelr,'.-'),legend('y','x'),title('Heel coordinates')
% figure,plot(toec,'.-')
% hold on
% plot(toer,'.-'),legend('y','x'),title('Toe coordinates')

dheelc = diff(heelc); %first derivative of y-coordinates (Heel)
dheelr = diff(heelr); %first derivative of x-coordinates (Heel)
dtoec = diff(toec); %first derivative of y-coordinates (Toe)
dtoer = diff(toer); %first derivative of x-coordinates (Toe)
Thcolonne = 2; %this parameter could be tuned based on the gait of the subject
Thrighe = 2; %this parameter could be tuned based on the gait of the subject

% figure,plot(dheelc),title('first derivative of y-coordinates - Heel')
% yline(Thcolonne)
% yline(-Thcolonne)
% figure,plot(dheelr),title('first derivative of x-coordinates - Heel')
% yline(Thrighe)
% yline(-Thrighe)
% figure,plot(dtoec),title('first derivative of y-coordinates  - Toe')
% yline(Thcolonne)
% yline(-Thcolonne)
% figure,plot(dtoer),title('first derivative of x-coordinates  - Toe')
% yline(Thrighe)
% yline(-Thrighe)


%% Heel - columns
contact_heel_c = (abs(dheelc)<=Thcolonne);
contact_heel_c = bwareaopen(contact_heel_c,3);
ind_contact_heel_c = find(contact_heel_c);
d1 = diff(ind_contact_heel_c);
fine_heel_c = find((d1>1));
int_heel_c = [ind_contact_heel_c(1) ind_contact_heel_c(fine_heel_c) ind_contact_heel_c(fine_heel_c+1) ind_contact_heel_c(end)];
int_heel_c = unique(sort(int_heel_c)); %instants of stationarity of the y-coordinates (heel)

% figure,plot(heelc)
% hold on
% xline(int_heel_c),title('Heel coordinates - x axis')

%% Heel - rows
contact_heel_r = (abs(dheelr)<=Thrighe);
contact_heel_r = bwareaopen(contact_heel_r,4);
ind_contact_heel_r = find(contact_heel_r);
d2 = diff(ind_contact_heel_r);
fine_heel_r = find((d2>1));
int_heel_r = [ind_contact_heel_r(1) ind_contact_heel_r(fine_heel_r)+1 ind_contact_heel_r(fine_heel_r+1) ind_contact_heel_r(end)];
int_heel_r = unique(sort(int_heel_r));  %instants of stationarity of the x-coordinates (heel)

% figure,plot(heelr)
% hold on
% xline(int_heel_r),title('Heel coordinates - y axis')

%% Toe - columns
contact_toe_c = (abs(dtoec)<=Thcolonne);
contact_toe_c = bwareaopen(contact_toe_c,3);
ind_contact_toe_c = find(contact_toe_c);
d3 = diff(ind_contact_toe_c);
fine_toe_c = find((d3>3));
int_toe_c = [ind_contact_toe_c(1) ind_contact_toe_c(fine_toe_c)+1 ind_contact_toe_c(fine_toe_c+1) ind_contact_toe_c(end)+1];
int_toe_c = unique(sort(int_toe_c));

% figure,plot(toec)
% hold on
% xline(int_toe_c),title('Toe coordinates - x axis')

int_heel_c = int_heel_c(1:2:end);
int_heel_r = int_heel_r(1:2:end);
inizio_nuovo_heel = [];
int_heel_c(int_heel_c==1) = [];
int_heel_r(int_heel_r==1) = [];
[~, ind_min] = min([length(int_heel_c) length(int_heel_r)]);
if ind_min == 1
    for i=1:length(int_heel_c)
        if int_heel_c(i)>=int_heel_r(i)
            inizio_nuovo_heel(i) = int_heel_c(i);
        else
            inizio_nuovo_heel(i) = int_heel_r(i);
        end
    end
else
    for i=1:length(int_heel_r)
        if int_heel_c(i)>=int_heel_r(i)
            inizio_nuovo_heel(i) = int_heel_c(i);
        else
            inizio_nuovo_heel(i) = int_heel_r(i);
        end
    end
end

for i =1:length(int_toe_c)-1
    vet_toe_r(i,:) = [toer(int_toe_c(i)) toer(int_toe_c(i)+1)];
    d_r(i,:) = diff(vet_toe_r(i,:));
end

inizio_nuovo_toe = int_toe_c;
vet_toe_nuovo_r = vet_toe_r;


for i=1:size(d_r,1)
    check = nnz(abs(d_r(i,:))>2);
    h=1;
    while check>=1 && int_toe_c(i)+h+1<=length(toer)
        vet_toe_nuovo_r(i,:) = [toer(int_toe_c(i)+h) toer(int_toe_c(i)+h+1)];
        d_r_nuovo(i,:) = diff(vet_toe_nuovo_r(i,:));
        check = nnz(abs(d_r_nuovo(i,:))>2);
        inizio_nuovo_toe(i) = int_toe_c(i)+h;
        h=h+1;
    end
end

inizio_nuovo_toe = inizio_nuovo_toe(1:2:end);
inizio_nuovo_heel(inizio_nuovo_heel==1) = [];
inizio_nuovo_toe(inizio_nuovo_toe==1) = [];
contatti_heel = inizio_nuovo_heel + start -1;
contatti_toe = inizio_nuovo_toe + start - 1;

if length(contatti_heel) <=2
    heel = contatti_heel(1);
else
    heel = contatti_heel(2);
end
if length(contatti_toe) <=2
    toe = contatti_toe(1);
else
    toe = contatti_toe(2);
end

if heel<=toe
    start1 = heel; %first frame of the gait cycle (heel contact)
    pos_IC1 = heelc(start1-start+1); %heel position in the first frame of the gait cycle
    if length(contatti_heel) == 2
        stop1 = contatti_heel(2); %last frame of the gait cycle (heel contact)
    else
        stop1 = contatti_heel(3); %last frame of the gait cycle (heel contact)
    end
    pos_IC2 = heelc(stop1-start+1); %heel position in the last frame of the gait cycle
else
    start1 = toe; %first frame of the gait cycle (toe contact)
    pos_IC1 = toec(start1-start+1); %toe position in the first frame of the gait cycle
    if length(contatti_toe) == 2
        stop1 = contatti_toe(2); %last frame of the gait cycle (toe contact)
    else
        stop1 = contatti_toe(3); %last frame of the gait cycle (toe contact)
    end
    pos_IC2 = toec(stop1-start+1); %toe position in the last frame of the gait cycle
end
%
%% Background foot (identification of IC#1)
heelc = [];
toec = [];
heelr = [];
toer = [];
ind_start = start1-start+1;
ind_stop = stop1-start+1;
heelc = heelc_c(ind_start:ind_stop);
heelr = heelr_c(ind_start:ind_stop);
toec = toec_c(ind_start:ind_stop);
toer = toer_c(ind_start:ind_stop);
set(0, 'DefaultFigureVisible', 'on')
% figure,plot(heelc,'.-')
% hold on
% plot(heelr,'.-'),legend('y','x'),title('MRF coordinates')
% figure,plot(toec,'.-')
% hold on
% plot(toer,'.-'),legend('y','x'),title('RF coordinates')

dheelc = diff(heelc);
dheelr = diff(heelr);
dtoec = diff(toec);
dtoer = diff(toer);
Thcolonne = 2;
Thrighe = 2;
% figure,plot(dheelc),title('differenza nella posizione delle colonne - Heel')
% yline(Thcolonne)
% yline(-Thcolonne)
% figure,plot(dheelr),title('differenza nella posizione delle righe - Heel')
% yline(Thrighe)
% yline(-Thrighe)
% figure,plot(dtoec),title('differenza nella posizione delle colonne - Toe')
% yline(Thcolonne)
% yline(-Thcolonne)
% figure,plot(dtoer),title('differenza nella posizione delle righe - Toe')
% yline(Thrighe)
% yline(-Thrighe)

%% MRF - columns
contact_heel_c = (abs(dheelc)<=Thcolonne);
contact_heel_c = bwareaopen(contact_heel_c,3);
ind_contact_heel_c = find(contact_heel_c);
d1 = diff(ind_contact_heel_c);
fine_heel_c = find((d1>3));
int_heel_c = [ind_contact_heel_c(1) ind_contact_heel_c(fine_heel_c) ind_contact_heel_c(fine_heel_c+1) ind_contact_heel_c(end)];
int_heel_c = unique(sort(int_heel_c));

% figure,plot(heelc)
% hold on
% xline(int_heel_c),title('Heel coordinates - x axis')

%% MRF - rows
contact_heel_r = (abs(dheelr)<=Thrighe);
contact_heel_r = bwareaopen(contact_heel_r,2);
ind_contact_heel_r = find(contact_heel_r);
d2 = diff(ind_contact_heel_r);
fine_heel_r = find((d2>1));
int_heel_r = [ind_contact_heel_r(1) ind_contact_heel_r(fine_heel_r)+1 ind_contact_heel_r(fine_heel_r+1) ind_contact_heel_r(end)];
int_heel_r = unique(sort(int_heel_r));

% figure,plot(heelr)
% hold on
% xline(int_heel_r),title('Heel coordinates - y axis')

%% Toe - columns
contact_toe_c = (abs(dtoec)<=Thcolonne);
contact_toe_c = bwareaopen(contact_toe_c,4);
ind_contact_toe_c = find(contact_toe_c);
d3 = diff(ind_contact_toe_c);
fine_toe_c = find((d3>1));
int_toe_c = [ind_contact_toe_c(1) ind_contact_toe_c(fine_toe_c)+1 ind_contact_toe_c(fine_toe_c+1) ind_contact_toe_c(end)+1];
int_toe_c = unique(sort(int_toe_c));

% figure,plot(toec)
% hold on
% xline(int_toe_c),title('Toe coordinates - x axis')

int_heel_c = int_heel_c(1:2:end);
int_heel_c(int_heel_c == 1) = [];
int_heel_r = int_heel_r(1:2:end);
int_heel_r(int_heel_r == 1) = [];
inizio_nuovo_heel = [];
[~, ind_min] = min([length(int_heel_c) length(int_heel_r )]);

if ind_min == 1
    for i=1:length(int_heel_c)
        if int_heel_c(i)>=int_heel_r(i)
            inizio_nuovo_heel(i) = int_heel_c(i);
        else
            inizio_nuovo_heel(i) = int_heel_r(i);
        end
    end
else
    for i=1:length(int_heel_r)
        if int_heel_c(i)>=int_heel_r(i)
            inizio_nuovo_heel(i) = int_heel_c(i);
        else
            inizio_nuovo_heel(i) = int_heel_r(i);
        end
    end
end

for i =1:length(int_toe_c)-1
    vet_toe_r(i,:) = [toer(int_toe_c(i)) toer(int_toe_c(i)+1)];
    d_r(i,:) = diff(vet_toe_r(i,:));
end

inizio_nuovo_toe = int_toe_c;
vet_toe_nuovo_r = vet_toe_r;


for i=1:size(d_r,1)
    check = nnz(abs(d_r(i,:))>2);
    h=1;
    while check>=1 && int_toe_c(i)+h+1<=length(toer)
        vet_toe_nuovo_r(i,:) = [toer(int_toe_c(i)+h) toer(int_toe_c(i)+h+1)];
        d_r_nuovo(i,:) = diff(vet_toe_nuovo_r(i,:));
        check = nnz(abs(d_r_nuovo(i,:))>2);
        inizio_nuovo_toe(i) = int_toe_c(i)+h;
        h=h+1;
    end
end

inizio_nuovo_toe = int_toe_c;
inizio_nuovo_toe = inizio_nuovo_toe(1:2:end);
contatti_heel = inizio_nuovo_heel + start -1 + ind_start-1;
contatti_toe = inizio_nuovo_toe + start - 1+ ind_start-1;
contatti_toe(contatti_toe <= start1) = [];
contatti_heel(contatti_heel <= start1) = [];
heel = contatti_heel(1);
toe = contatti_toe(1);

if heel<=toe
    start2 = heel;
    pos_IC1_c = heelc(start2-start+1 - ind_start+1);
else
    start2 = toe;
    pos_IC1_c = toec(start2-start+1 - ind_start+1);
end

stride_length_pixel = abs(pos_IC1 - pos_IC2); %pixel
step_length_pixel = abs(pos_IC1-pos_IC1_c);
stride_length = abs(pos_IC1 - pos_IC2)*conv_factor; %m/s
step_length = abs(pos_IC1-pos_IC1_c)*conv_factor;
gait_speed = stride_length/(stop1-start1);
end
