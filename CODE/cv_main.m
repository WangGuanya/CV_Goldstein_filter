function [ CV, cv ] = cv_main( left, block, step )
% Calculate CV for the whole interferogram
%   [ CV ]= CV_main( IFG, block, step )
%
%   Guanya Wang, 2017/11/08
%
% ================= Input ====================
%
% IFG     :  Input interferogram (in complex format)
% block   :  Filtering window, should be even a number, such as 32, 64, â€¦â?, default 32
% step    :  the step of sliding window
%
%
% ================ Output ===================
% CV   : a matrix of the coefficient of despresion for the whole IFG



% Handle Inputs
% if nargin < 2;      block = 32;        end    % set the filtering window if the 2nd parameter not set.
% if nargin < 3;      step = block/8;    end    % set the window sliding step if the 4th parameter not set.
% if nargin < 4;      L = 5 ;            end    % set the look number if the 3rd parameter not set.
% if isreal(IFG);     IFG = exp(1j*IFG); end    % real image to complex image.

flag = fix(log2(block))-log2(block);          % checking whether the filtering window is 2*N.
if flag; error('The filtering window is not 2*N! Please Check it!'); end
if step > block; error('ERROR: the step larger than the filtering window!'); end

[lines, width] = size(left);          % get the size of input inmage

overlap = (block - step)/2;            % if N_win is 32 and step is 4, the overlap is 14.

Start_Col = overlap+1;
End_Col = width-(block-overlap)+1;

Start_Row = overlap+1;
End_Row = lines-(block-overlap)+1;

m = (End_Col - Start_Col)/4 + 1;
n = (End_Row - Start_Col)/4 + 1;

CV = zeros(m+2, n+2);

% ====== up-left corner ===================================================
window = left ( 1:block,1:block );
[ cv ] = cv_calculation( window, block );
CV(1,1) = cv;
% =========================================================================

% ====== up-right corner ==================================================
window = left (1:block,End_Col-overlap:width);
[ cv ] = cv_calculation( window, block );
CV(1,n+2) = cv;
% =========================================================================

% ====== down-left corner =================================================
window = left (End_Row-overlap:lines,1:block);
[ cv ] = cv_calculation( window, block );
CV(m+2,1) = cv;
% =========================================================================

% ====== down-right corner ================================================
window = left (End_Row-overlap:lines,End_Col-overlap:width);
[ cv ] = cv_calculation( window, block );
CV(m+2,n+2) = cv;
% =========================================================================

% ====== Filter the first and last row ====================================
nn = 2;

for jj=Start_Col:step:End_Col
      
    % ====== first row ======
      window = left(1:block,jj-overlap:jj+(block-overlap)-1);
      [ cv ] = cv_calculation( window, block );
      CV(1,nn) = cv;
         
    % ====== last row ======
      window = left(1:block,jj-overlap:jj+(block-overlap)-1);
      [ cv ] = cv_calculation( window, block );
      CV(m+2,nn) = cv;
      
      nn = nn+1;
 end
% =========================================================================

% ====== Filter the first and last column =================================
mm = 2;
 
for ii=Start_Row:step:End_Row
    
   % ====== first column ======
    window = left(ii-overlap:ii+(block-overlap)-1,1:block);  
    [ cv ] = cv_calculation( window, block );
    CV(mm,1) = cv;
    
   % ====== last column ======
    window = left(ii-overlap:ii+(block-overlap)-1,width-block+1:width);
    [ cv ] = cv_calculation( window, block );
    CV(mm,n+2) = cv;
    
    mm = mm+1;
end
% =========================================================================

% ====== Calculate the center part ========================================
mm = 2;
for ii=Start_Row:step:End_Row
    nn = 2;
    for jj=Start_Col:step:End_Col

    window = left(ii-overlap:ii+(block-overlap)-1,jj-overlap:jj+(block-overlap)-1);  
    cv = cv_calculation ( window, block );
    CV(mm,nn) = cv;
    
    nn = nn+1;    
    end  
    mm = mm+1;
end
% =========================================================================

cv = reshape (CV, (m+2)*(n+2), 1);

end
