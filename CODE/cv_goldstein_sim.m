function [ Filtered_IFG, B, CV ]= cv_goldstein_sim( IFG, left, CC_map, block, step )
% Modified Goldstein filter using CV 
%
% Guanya Wang, 2017/11/09  2017/12/15
%
% ================= Input ====================
% IFG     :  Input interferogram
% left    :  Input the interferogram without fringes
% CC_map  :  Input coherence map
% block   :  Filtering window, should be even a number, such as 32 or 64,s default 32
% step    :  the step of sliding window
%
% ================ Output ====================
% Filtered_IFG   : Filtered interferogram (in complex format)
% B              : the coefficients of poly_fit
% CV             : the coefficient variation map of IFG


% ================ Handle Inputs ==============
if nargin < 4;      block= 32;        end    % set the filtering window if the 3rd parameter not set.
if nargin < 5;      step= block/8;    end    % set the window sliding step if the 4th parameter not set.
if isreal(IFG);     IFG= exp(1j*IFG); end    % real image to complex image.
if isreal(left);    left= exp(1j*left); end    % real image to complex image.

flag= fix(log2(block))-log2(block);          % checking whether the filtering window is 2*N.
if flag; error('The filtering window is not 2*N! Please Check it!'); end
if step > block; error('ERROR: the step larger than the filtering window!'); end


% ================ CV ==============================
[ CV, cv ] = cv_main( left, block, step );

% ================ CV_alpha_fit ====================
[Size_IFG, Size_CC, B, R, J, Cov, MSE ] = cv_alpha_fit (left, CC_map, block, step );

B(2) = B(2) ;
B(1) = roundn(B(1),-2);
B(2) = roundn(B(2),-2);
B(3) = roundn(B(3),-2);


% ================ Set parameter ===============
lines = Size_IFG(1);
width = Size_IFG(2);

Filtered_IFG = IFG;
           

% ================ Calculate filtering image size ======
overlap= (block-step)/2;          % if N_win is 32 and step is 4, the overlap is 14.

Start_Col= overlap+1;
End_Col= width-(block-overlap)+1;

Start_Row= overlap+1;
End_Row= lines-(block-overlap)+1;


% ================ Kernel ===============================
N=5;
K= ones(N,N);
K=K/sum(K(:));
K=fftshift(fft2(K));

% =======================================================
% ================ Filter ===============================
% =======================================================

% ====== up-left corner ======
window = IFG(1:block,1:block);          % cuting IFG
H = fft2(window);                       % 2D FFT
H = fftshift(H);
S = conv2(abs(H),K,'same');             % 2D Convolution: smoothing response spectum with Kernel K.
% S = S ./max(S(:));                       % scale of the S [0 1].
S = (S - min(S(:)))./(max(S(:))-min(S(:)));                       % scale of the S [0 1].

% cv = CV(1,1) 
% cv = cv_calculation (S);
tmp= reshape(S, block*block, 1);
tmp = abs(tmp);
u= mean(tmp);                       % Calculate the mean value in window
sigma= std(tmp);                    % Calculate the standard deviation of window
cv= sigma/u;
alpha = B(1)*((cv+B(3))^(-1)) + B(2);      

S = S.^alpha;                           % filtering
H = H.*S;
H = ifftshift(H);                       % inverse FFT
H1 = ifft2(H);                            
Filtered_IFG(1:block,1:block) = H1;       % recode
% ===============================


% ====== up-right corner ======
window = IFG(1:block,End_Col-overlap:width);           % cuting IFG
H = fft2(window);                                       % 2D FFT
H = fftshift(H);
S = conv2(abs(H),K,'same');                        % 2D Convolution: smoothing response spectum with Kernel K.
% S = S./max(S(:));                                  % scale of the S [0 1].
S = (S - min(S(:)))./(max(S(:))-min(S(:))); 

% cv = CV(1,n+2);
% cv = cv_calculation (S);
tmp= reshape(S, block*block, 1);
tmp = abs(tmp);
u= mean(tmp);                       % Calculate the mean value in window
sigma= std(tmp);                    % Calculate the standard deviation of window
cv= sigma/u;
alpha = B(1)*((cv+B(3))^(-1)) + B(2);   

S = S.^alpha; 
H = H.*S;
H = ifftshift(H);
H1 = ifft2(H);
Filtered_IFG(1:block,End_Col-overlap:width) = H1;
% ===============================


% ====== down-left corner ======
window = IFG(End_Row-overlap:lines,1:block);            % cuting IFG
H = fft2(window);                                       % 2D FFT
H = fftshift(H);
S = conv2(abs(H),K,'same');                        % 2D Convolution: smoothing response spectum with Kernel K.
% S = S./max(S(:));                                  % scale of the S [0 1].
S = (S - min(S(:)))./(max(S(:))-min(S(:))); 

% cv = CV(m+2,1);
% cv = cv_calculation (S);
tmp= reshape(S, block*block, 1);
tmp = abs(tmp);
u= mean(tmp);                       % Calculate the mean value in window
sigma= std(tmp);                    % Calculate the standard deviation of window
cv= sigma/u;
alpha = B(1)*((cv+B(3))^(-1)) + B(2);   

S = S.^alpha; 
H = H.*S;
H = ifftshift(H);
H1 = ifft2(H);
Filtered_IFG(End_Row-overlap:lines,1:block) = H1;
% ===============================


% ====== down-right corner ======
window = IFG(End_Row-overlap:lines,End_Col-overlap:width);            % cuting IFG
H = fft2(window);                                       % 2D FFT
H = fftshift(H);
S = conv2(abs(H),K,'same');                        % 2D Convolution: smoothing response spectum with Kernel K.
% S = S./max(S(:));                                  % scale of the S [0 1].
S = (S - min(S(:)))./(max(S(:))-min(S(:))); 

% cv = CV(m+2,n+2);
% cv = cv_calculation (S);
tmp= reshape(S, block*block, 1);
tmp = abs(tmp);
u= mean(tmp);                       % Calculate the mean value in window
sigma= std(tmp);                    % Calculate the standard deviation of window
cv= sigma/u;
alpha = B(1)*((cv+B(3))^(-1)) + B(2);   

S = S.^alpha;
H = H.*S;
H = ifftshift(H);
H1 = ifft2(H);
Filtered_IFG(End_Row-overlap:lines,End_Col-overlap:width) = H1;
% ===============================

% ===========================================
% ====== Filter the first and last row ======
% ===========================================
for jj = Start_Col:step:End_Col

% ====== first row ======
    window = IFG(1:block,jj-overlap:jj+(block-overlap)-1);  
    H = fft2(window);
    H = fftshift(H);
    S = conv2(abs(H),K,'same');
%     S = S./max(S(:));
    S = (S - min(S(:)))./(max(S(:))-min(S(:))); 
    
%     cv = CV(1,nn);
%     cv = cv_calculation (S);
    tmp= reshape(S, block*block, 1);
    tmp = abs(tmp);
    u= mean(tmp);                       % Calculate the mean value in window
    sigma= std(tmp);                    % Calculate the standard deviation of window
    cv= sigma/u;
    alpha = B(1)*((cv+B(3))^(-1)) + B(2);   
    
    S = S.^alpha; 
    H = H.*S;
    H = ifftshift(H);
    H1 = ifft2(H);
    Filtered_IFG(1:block,jj:jj+step-1) = H1(:,overlap+1:block-overlap);
    
% ====== last row ======
    window = IFG(lines-block+1:lines,jj-overlap:jj+(block-overlap)-1);
    H = fft2(window);
    H = fftshift(H);
    S = conv2(abs(H),K,'same');
%     S = S./max(S(:));
    S = (S - min(S(:)))./(max(S(:))-min(S(:))); 
    
%     cv = CV(m+2,nn);
%     cv = cv_calculation (S);
    tmp= reshape(S, block*block, 1);
    tmp = abs(tmp);
    u= mean(tmp);                       % Calculate the mean value in window
    sigma= std(tmp);                    % Calculate the standard deviation of window
    cv= sigma/u;
    alpha = B(1)*((cv+B(3))^(-1)) + B(2);   
    
    S = S.^alpha;
    H = H.*S;
    H = ifftshift(H);
    H1 = ifft2(H);
    Filtered_IFG(lines-block+1:lines,jj:jj+step-1) = H1(:,overlap+1:block-overlap);

end


% ===========================================
% ====== Filter the first and last column ===
% ===========================================
for ii = Start_Row:step:End_Row
    
% ====== first column ======
    window = IFG(ii-overlap:ii+(block-overlap)-1,1:block);  
    H = fft2(window);
    H = fftshift(H);
    S = conv2(abs(H),K,'same');
%     S = S./max(S(:));
    S = (S - min(S(:)))./(max(S(:))-min(S(:)));  
    
%     cv = CV(mm,1);
%     cv = cv_calculation (S);
    tmp= reshape(S, block*block, 1);
    tmp = abs(tmp);
    u= mean(tmp);                       % Calculate the mean value in window
    sigma= std(tmp);                    % Calculate the standard deviation of window
    cv= sigma/u;
    alpha = B(1)*((cv+B(3))^(-1)) + B(2);   
    
    S = S.^alpha; 
    H = H.*S;
    H = ifftshift(H);
    H1 = ifft2(H);
    Filtered_IFG(ii:ii+step-1,1:block) = H1(overlap+1:block-overlap,:);
    
% ====== last column ======
    window = IFG(ii-overlap:ii+(block-overlap)-1,width-block+1:width);
    H = fft2(window);
    H = fftshift(H);
    S = conv2(abs(H),K,'same');
%     S = S./max(S(:));
    S = (S - min(S(:)))./(max(S(:))-min(S(:))); 
    
%     cv = CV(mm,n+2);
%     cv = cv_calculation (S);
    tmp= reshape(S, block*block, 1);
    tmp = abs(tmp);
    u= mean(tmp);                       % Calculate the mean value in window
    sigma= std(tmp);                    % Calculate the standard deviation of window
    cv= sigma/u;
    alpha = B(1)*((cv+B(3))^(-1)) + B(2);   
    
    S = S.^alpha;
    H = H.*S;
    H = ifftshift(H);
    H1 = ifft2(H);
    Filtered_IFG(ii:ii+step-1,width-block+1:width)=H1(overlap+1:block-overlap,:);

end

% ===========================================
% ====== Filter the center part =============
% ===========================================
for ii = Start_Row:step:End_Row
    for jj = Start_Col:step:End_Col

    window = IFG(ii-overlap:ii+(block-overlap)-1,jj-overlap:jj+(block-overlap)-1);  
    H = fft2(window);
    H = fftshift(H);
    S = conv2(abs(H),K,'same');
%     S = S./max(S(:));
    S = (S - min(S(:)))./(max(S(:))-min(S(:))); 
    
%     cv = CV(mm,nn);
%     cv = cv_calculation (S);
    tmp= reshape(S, block*block, 1);
    tmp = abs(tmp);
    u= mean(tmp);                       % Calculate the mean value in window
    sigma= std(tmp);                    % Calculate the standard deviation of window
    cv= sigma/u;
    alpha = B(1)*((cv+B(3))^(-1)) + B(2);   
    
    S = S.^alpha; 
    H = H.*S;
    H = ifftshift(H);
    H1 = ifft2(H);
    Filtered_IFG(ii:ii+step-1,jj:jj+step-1)=H1(overlap+1:block-overlap,overlap+1:block-overlap);    
    end  
end


end
