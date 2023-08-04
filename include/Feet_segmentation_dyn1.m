function [D_imm_l_t, D_imm_r_t, I1] = Feet_segmentation_dyn1(I)
%Author: Diletta Balta
%Department of Electronics and Telecommunications
%Politecnico di Torino
%diletta.balta@polito.it

%This function provides two feet segmentation masks for each frame by
%implementing two color filters (red for right foot and blue for left foot)

%inputs
% I = RGB image

%outputs
%D_imm_l_t = left foot segmentation mask
%D_imm_r_t = right foot segmentation mask
%I1 = segmentation mask containing both feet

D_imm_r_t = false(size(I,1),size(I,2));
D_imm_l_t = false(size(I,1),size(I,2));


% figure
% imshow(I)

%Based on the environment/lights conditions or the colors of socks, the thresholds for the color filters have to be adjusted

for r = 480:600 %rows
    for c = 165:1100 % columns
        pix = [I(r,c,1) I(r,c,2) I(r,c,3)];
        if  pix(1)<60 && pix(2)>=40 && pix(2)<=100 && pix(3)>90 && pix(3)<=250 %left foot blue
            D_imm_l_t(r,c) = true;
        end
        if  pix(1)>30 && pix(2)<30 && pix(3)<40 %right foot red
            D_imm_r_t(r,c) = true;
        end
    end
end
% hold on
% imshow(D_imm_l_t)
I1 = zeros(size(I));
I1(D_imm_r_t)=1;
I1(D_imm_l_t)=1;

end
