%% GLaucoma Detection 
 
clc; 
clear all; 
close all; 
 
addpath('C:\Users\Kaushal\Desktop\glaucoma and healthy image dataset'); 
 
 
[f p] = uigetfile('*.jpg'); 
I=imread([p f]); 
figure 
imshow(I); 
title('INPUT FUNDUS IMAGE'); 
 
I1=im2bw(I,.7); 
figure 
imshow(I1); 
se = strel('disk',10); 
I1=imclose(I1,se); 
I1 = bwlabel(I1, 8); 
I1=bwareaopen(I1,400); 
figure 
imshow(I1); 
stats=regionprops(I1,'BoundingBox'); 
 
 hold on 
    for object = 1:length(stats) 
        bb = stats(object).BoundingBox; 
         
        rectangle('Position',bb,'EdgeColor','r','LineWidth',2) 
        bb(1:2)=bb(1:2)-70; 
        bb(3:4)=bb(3:4)+140; 
        rectangle('Position',bb,'EdgeColor','y','LineWidth',2) 
    end 
    hold off 
     
I2 = imcrop(I,bb); 
I2=imresize(I2,[220,220]); 
figure 
imshow(I2); 
title('Region of intrest ROI Image'); 
 
% Extracting the 3 channels red,green and blue 
Rd = I2(:,:,1); 
Gn = I2(:,:,2); 
Bl = I2(:,:,3); 
 
figure 
imshow(Rd); 
title('RED CHANNEL OF ROI'); 
figure 
imshow(Gn); 
title('GREEN CHANNEL OF ROI'); 
figure 
imshow(Bl); 
title('BLUE CHANNEL OF ROI'); 
 
%% Segmentation of ROI Images and measuring NRR 
 
% Segmentation of RED channel 
n = 2 
IDX = otsu(Rd,n); 
    
figure 
imagesc(IDX), axis image off 
title(['Segmented Red Channel with n=' int2str(n)],'FontWeight','bold') 
colormap(gray) 
  
% Segmentation of GREEN channel 
n = 2 
idx = otsu(Gn,n); 
    
figure 
imagesc(idx), axis image off 
title(['Segmented Green Channel with n = ' int2str(n)],'FontWeight','bold') 
colormap(gray) 
 
 
 
[R C] = size(Rd); 
  for i=1:R 
      for j=1:C 
          if(idx(i,j)==1) 
              idx(i,j)=0; 
          end 
          if(idx(i,j)==2) 
              idx(i,j)=1; 
          end 
      end 
  end 
figure 
imshow(idx); 
title('Binary image of idx'); 
 
 for i=1:R 
      for j=1:C 
          if(IDX(i,j)==1) 
              IDX(i,j)=0; 
          end 
          if(IDX(i,j)==2) 
              IDX(i,j)=1; 
          end 
      end 
  end 
figure 
imshow(IDX); 
title('Binary image of IDX'); 
 
 
 
% determining the centroids for ISNT generation 
 
 stats = regionprops(IDX,'centroid'); 
 centroids = cat(1,stats.Centroid); 
 hold on; 
 plot(centroids(1,1), centroids(1,2), 'r+', 'LineWidth', 1); 
 hold off; 
 
 cC = int32(centroids(1));  
 cR = int32(centroids(1)); 
 
cC = R/2;  
cR = C/2; 
 
disp('Centroid Row='); 
disp(cR); 
disp('Centroid Coloumn='); 
disp(cC); 
 

 
% Measuring the NeuroRetinalRim 
NRR = IDX - idx; 
 
figure 
imagesc(NRR); 
title('Neuro Retinal Rim Image'); 
colormap(gray) 
[R C] = size(NRR); 
 
 for i=1:R 
      for j=1:C 
          if(NRR(i,j)==1) 
              NRR(i,j)=1; 
          else 
              NRR(i,j)=0; 
          end 
      end 
  end 
figure 
imshow(NRR); 
title('Binary image of NRR'); 
 
 
 
%% generating the ISNT Mask And Measuring ISNT NRR 
 
% Inferior region mask 
[R C] = size(NRR); 
disp('ISNT Row='); 
disp(R); 
disp('ISNT coloumn='); 
disp(C); 
 
c1 = cC; 
c2 = cR; 
c1=round(c1); 
c2=round(c2); 
% c1 = 93; 
% c2 = 95; 
z=zeros(R,C); 
for i = c2:R 
    for j = 1:C 
         
        if((i+j)<(c1+c2)||(i<j)) 
     
                z(i,j) = 0; 
   
        elseif ((i+j)<(c1+c2)&&((j-1)>(c1-c2))) 
                 
                z(i,j) = 0; 
             
        elseif (((i+j)<(c1+c2))||(i<j)||(i==j)||((i-j)<(c2-c1))) 
                     
                 z(i,j) = 0; 
        else 
                z(i,j) = 1; 
        
        end 
    end 
end 
 
figure 
imshow(z); 
title('Inferior Region MASK'); 
 
% inferior region of NRR 
X = NRR.*z; 
figure 
imagesc(X), axis image off 
title('Inferior Region NRR'); 
colormap(gray); 
 
% Superior region of mask 
z1 = zeros(R,C); 
 
for i = 1:c2 
    for j = 1:C 
         
        if ((i+j)>(c1+c2)||(i>j)) 
             
           z1(i,j) = 0; 
        elseif (((i+j)>(c1+c2))||(i>j)||(i==j)||((j-1)<(c1-c2))) 
             
            z1(i,j) = 0; 
             
        elseif(((i+j)>(c1+c2))||(i>j)&((i-j)>(c2-c1))) 
             
            z1(i,j) = 0; 
             
        else 
             
            z1(i,j) = 1; 
             
        end 
    end 
end 
 
 
figure 
imshow(z1); 
title('Superior Region MASK'); 
 
% superior region of NRR 
X1 = NRR.*z1; 
figure 
imagesc(X1), axis image off 
title('Superior Region NRR'); 
colormap(gray); 
 
% Temporal region of mask 
z2 = zeros(R,C); 
 
for i = 1:R 
    for j = 1:c1 
         
        if (((i+j)>(c1+c2))||(i<j)) 
             
            z2(i,j) = 0; 
             
        elseif(((i+j)>(c1+c2))||(i<j)&((j-1)>(c1-c2))) 
             
             
            z2(i,j) = 0; 
             
        elseif(((i+j)>(c1+c2))||(i<j)||(i==j)||((i-j)<(c2-c1))) 
             
            z2(i,j) = 0; 
             
        else 
             
            z2(i,j) = 1; 
             
        end 
    end 
end 
 
figure 
imshow(z2); 
title('Temporal Region MASK'); 
 
% Temporal region of NRR 
X2 = NRR.*z2; 
figure 
imagesc(X2), axis image off 
title('Temporal Region NRR'); 
colormap(gray);           
             
% Nasal Region of mask 
z3 = zeros(R,C); 
 
for i= 1:R 
    for j = c1:C 
         
        if(((i+j)<(c1+c2))||(i>j)) 
             
            z3(i,j) = 0; 
             
        elseif(((i+j)<(c1+c2))||(i>j)||(i==j)||((j-1)<(c1-c2))) 
             
            z3(i,j) = 0; 
             
        elseif(((i+j)<(c1+c2))||(i>j)&((i-j)>(c2-c1))) 
             
            z3(i,j) = 0; 
             
        else 
             
            z3(i,j) = 1; 
             
        end 
    end 
end 
 
figure 
imshow(z3); 
title('Nasal Region MASK'); 
 
% Nasal region of NRR 
X3 = NRR.*z3; 
figure 
imagesc(X3), axis image off 
title('Nasal Region NRR'); 
colormap(gray);  
 

%% Determining the CDR and RDR 
 
% Calculating the area of optic Disc 
count = 0; 
for i = 1:R 
    for j =1:C 
         
        if(IDX(i,j)==1) 
             
            count = count+1; 
             
        end 
    end 
end 
 
disp('Optic disc area='); 
disp(count); 
 
% Calculating the area of optic cup 
count1 = 0; 
 
for i = 1:R 
    for j =1:C 
         
        if(idx(i,j)==1) 
             
            count1 = count1+1; 
             
        end 
    end 
end 
 
disp('Area of optic Cup='); 
disp(count1); 
 
% Measuring the cup to disc ratio 
CDR = (count1/count); 
disp('CUP TO DISC RATIO='); 
disp(CDR); 
 
% Calculating the inferior region area 
inf_ar =0; 
 
for i = 1:R 
    for j =1:C 
         
        if(X(i,j)==1) 
            
            inf_ar = (inf_ar+1); 
        end 
         
    end 
end 
 
disp('inferior region Area='); 
disp(inf_ar); 
 
% calculating the temporal region area 
temp_ar = 0; 
 
for i =1:R 
    for j=1:C 
         
        if(X2(i,j)==1) 
             
            temp_ar = (temp_ar+1); 
         
    end 
end 
end 
 
disp('Temporal region area='); 
disp(temp_ar); 
 
% inferior-temporal region area 
inte_ar = (inf_ar+temp_ar); 
disp('infer_temporal area='); 
disp(inte_ar); 
 
nasal_ar = 0; 
 
for i =1:R 
    for j=1:C 
         
        if(X3(i,j)==1) 
             
            nasal_ar = (nasal_ar+1); 
         
    end 
end 
end 
disp('Nasal area='); 
disp(nasal_ar); 
 
inna_ar = (inf_ar+nasal_ar); 
 
% Measuring the RIM TO DISC RATIO 
RDR = (inte_ar/count); 
disp('CUP TO DISC RATIO='); 
disp(CDR); 
disp('RIM TO DISC RATIO='); 
disp(RDR); 
 
RDR1 = inna_ar/count; 
disp(RDR1); 
 
% Displaying output as Glaucomatic or Not 
if ((CDR>0.3)&(RDR<=0.35)) 
     
    disp('Patient is Glaucomatic'); 
else 
    disp('Patient is Non_Glaucomatic'); 
     
end 
if ((CDR>0.3)&(RDR1<0.35)) 
     
    disp('Patient is Glaucomatic'); 
else 
    disp('Patient is Non_Glaucomatic'); 
     
end 
 
 