function [ PSD, residue, SPD ]= index_model ( Filtered_IFG)
% Quality evaluation index
%  [ PSD, RMS, EPI, residue, residuez, SPD ]= index ( Filtered_IFG )
%
%   Guanya Wang, 2017/11/09
%
% ================= Input ====================
% Filtered_IFG: filtered interferogram
%
%
% ================ Output ===================
% u       : The mean value of the Filtered_IFG.
% sigma   : The variance of the Filtered_IFG.
% ENL     : The equivalent number of looks.
% PSD     : Phase standard deviation.
% RMS     : Root mean square.
% EPI     : Edge preserve index.
% residue : the number of residue points.
% residuez: the number of residue > 0.
% SPD     : the sum of phase gradient.

I = Filtered_IFG;
I_angle = angle(I);

[rows,cols] = size(I);


% % =============== u / sigma ================
% u = mean2(I_abs);
% sigma = std2(I_abs);
%
% % ============== ENL =======================
% ENL = (u/sigma)^2;


% ============== PSD =======================
B=zeros(rows-2,cols-2);
phai = I_angle;
for a=2:rows-1
    for b=2:cols-1
        temp1=(1/9)*(phai(a-1,b-1)+phai(a-1,b)+phai(a-1,b+1)...
            +phai(a,b-1)+phai(a,b)+phai(a,b+1)...
            +phai(a+1,b-1)+phai(a+1,b)+phai(a+1,b+1));
        temp2=phai(a,b)-temp1;
        B(a,b)=temp2.^2;
    end
end

A=zeros(rows-2,cols-2);
for c=3:rows-2
    for d=3:cols-2
        temp3=B(c-1,d-1)+B(c-1,d)+B(c-1,d+1)...
            +B(c,d-1)+B(c,d)+B(c,d+1)...
            +B(c+1,d-1)+B(c+1,d)+B(c+1,d+1);
        A(c,d)=sqrt(temp3/8);
    end 
end
PSD = A;

% 
% % ============= RMS =======================
% tmp_D = (I_angle - I0_angle).^2;
% tmp_DD = sum(tmp_D(:));
% RMS = sqrt(tmp_DD/(rows*cols-1));


% % ============= EPI =======================
% tmp_4 = 0;
% tmp_5 = 0;
% 
% for i=1:rows-1
%     for j=1:cols-1
%         
%         tmp_4 = tmp_4 + abs(I_abs(i,j)-I_abs(i+1,j)) + abs(I_abs(i,j)-I_abs(i,j+1));
%         tmp_5 = tmp_5 + abs(I0_abs(i,j)-I0_abs(i+1,j)) + abs(I0_abs(i,j)-I0_abs(i,j+1));
%         
%     end
% end
% 
% EPI = tmp_4/tmp_5;


% =========== residue =====================
tmp_6 = 0;
tmp_8 = 0;

for m = 1:rows-1
    for n = 1:cols-1
         tmp_7 = wrap(I_angle(m+1,n)-I_angle(m,n))...
               + wrap(I_angle(m+1,n+1)-I_angle(m+1,n))...
               + wrap(I_angle(m,n+1)-I_angle(m+1,n+1))...
               + wrap(I_angle(m,n)-I_angle(m,n+1));
%          tmp_7 = (I_angle(m+1,n)-I_angle(m,n))...
%                + (I_angle(m+1,n+1)-I_angle(m+1,n))...
%                + (I_angle(m,n+1)-I_angle(m+1,n+1))...
%                + (I_angle(m,n)-I_angle(m,n+1));

         if (abs(tmp_7) > 0.0000001)
             tmp_6 = tmp_6 +1;
         end
         
         if (tmp_7 > 0.0000001)
             tmp_8 = tmp_8 +1;
         end
         
    end
end

residue = tmp_6;
residuez = tmp_8;

% =================== SPD ====================
SPD = 0;
for a = 2:rows-1
    for b = 2:cols-1
        APD=(1/8)*(abs(I_angle(a,b)-I_angle(a-1,b-1))...
            +abs(I_angle(a,b)-I_angle(a-1,b))...
            +abs(I_angle(a,b)-I_angle(a-1,b+1))...
            +abs(I_angle(a,b)-I_angle(a,b-1))...
            +abs(I_angle(a,b)-I_angle(a,b))...
            +abs(I_angle(a,b)-I_angle(a,b+1))...
            +abs(I_angle(a,b)-I_angle(a+1,b-1))...
            +abs(I_angle(a,b)-I_angle(a+1,b))...
            +abs(I_angle(a,b)-I_angle(a+1,b+1)));
        SPD = SPD + APD;
    end
end



end



      
        
        
        
        
        
        
        
        
        
        
        