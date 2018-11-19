function [ cv ]= cv_calculation ( window, block )
% Calculate coefficient of desperision in one window by a same-size block
%   [ cv ]= cv_calculation ( window, block )
%
%   Guanya Wang, 2017/11/07
%
% ================= Input ====================
%
% window  :  window_image in interferogram for cv_calculation
% block   :  Filtering window, should be even a number, such as 32, 64, …�?, default 32
%
%
% ================ Output ===================
% cv   : coefficient of desperision



% Handle Inputs
if nargin < 2;      block = 32;        end    % set the filtering window if the 2nd parameter not set.
if nargin < 3;      L = 5 ;            end    % set the look number if the 3rd parameter not set.
if isreal(window);  window= exp(1j*window); end    % real image to complex image.

flag= fix(log2(block))-log2(block);          % checking whether the filtering window is 2*N.
if flag; error('The filtering window is not 2*N! Please Check it!'); end

[lines, width]=size(window);          % get the size of window


% ======= Kernel ======
H= window(1:block,1:block);         % cuting window
H= fft2(H);                         % 2D FFT
H= fftshift(H);

N=3;
K= ones(N,N);
K=K/sum(K(:));
K=fftshift(fft2(K));

S= conv2(abs(H),K,'same');          % 2D Convolution: smoothing response spectum with Kernel K.
% S = S ./max(S(:));   
S = (S - min(S(:)))./(max(S(:))-min(S(:)));                    % scale of the S [0 1].

tmp= reshape(S, block*block, 1);
tmp = abs(tmp);

u= mean(tmp);                       % Calculate the mean value in window
sigma= std(tmp);                    % Calculate the standard deviation of window

cv= sigma/u;


end

