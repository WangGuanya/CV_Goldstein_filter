% 条纹重构函数
function [F_IFG] = fringe_reconstruction_corase(IFG, block, overlap)
%ibgoldstein -- Script to filter interferogram using modified
%goldstein algorythm.
%
%e.g. [F_IFG] = ibgoldstein(1, IFG, coh, 32, 15);
%
%mode 1/0 - 1 - modified Goldstein, 0 - unmodified Goldstein
%IFG - interferogram (complex)
%coh - coherence map (float)
%block - patch in the interferogram (IFG) to calculate the spectrum e.g. 32 ( 32x32 pixels)
%overlap - in pixels e.g. 10, but max: (block/2)-1
%F_IFG filtered interferogram (complex)
%
%This script does have the edge effect that needs to be improved
%**********************************************************
%After Goldstein and Werner,
%Radar interferogram Estimate for geophisical
%applications, pp. 4035,-4038, 1998.
%Irek's idea is to change the constant character of alpha
%and depend it from coherence
%//Irek Baran 10/08/02
%//Modified:  28/05/03
%**********************************************************
1;
%%% Handle input.
if isreal(IFG)
    IFG = exp(j*IFG);
end

if nargin < 3
    help ibgoldstein;
    return;
end
%---------------------------------------------------------

%Goldstein filter-----------------------------------------
p=1;
disp('Start Estimate...');
[L P] = size(IFG);
all = (L*P);
all_step = floor(all/10);

F_IFG = IFG;
first_F_IFG = overlap + 1;
step = block - (2 * overlap);
IFG_L_stop = floor(L - (block - overlap) + 1);
IFG_P_stop = floor(P - (block - overlap) + 1);


%processing of last corner
jj = P - (block-overlap) + 1;
ii = L - (block-overlap) + 1;

H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
[H1 , ~, ~] = fringe_estimate_coarse(H);
F_IFG(ii:L,jj:P) = H1(overlap+1:block, overlap+1:block);


%precessing of first row
disp('Estimate First row...');
ii = first_F_IFG;
for jj = first_F_IFG:step:IFG_P_stop
    H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
    [H1 , ~, ~] = fringe_estimate_coarse(H);
    if jj==first_F_IFG
        F_IFG(1:ii+step-1,1:jj+step-1) = H1(1:block-overlap, 1:block-overlap);
    else
        F_IFG(1:ii+step-1,jj:jj+step-1) = H1(1:block-overlap, overlap+1:block-overlap);
    end
end


%precessing of last row
disp('Estimate last row...');
ii = L-(block-overlap)+1;
for jj = first_F_IFG:step:IFG_P_stop
    H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
    [H1 , ~, ~] = fringe_estimate_coarse(H);
    if jj==first_F_IFG
        F_IFG(ii:L,1:jj+step-1) = H1(overlap+1:block, 1:block-overlap);
    else
        F_IFG(ii:L,jj:jj+step-1) = H1(overlap+1:block, overlap+1:block-overlap);
    end
end


%processing of first column
disp('Estimate first column...');
jj = first_F_IFG;
for ii = first_F_IFG+step:step:IFG_L_stop 				%filtered interferogram
    H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
    [H1 , ~, ~] = fringe_estimate_coarse(H);
    F_IFG(ii:ii+step-1,1:jj+step-1) = H1(overlap+1:block-overlap, 1:block-overlap);
end


%processing of last column
disp('Estimate last column...');
jj = P-(block-overlap)+1;
for ii = first_F_IFG:step:IFG_L_stop 				%filtered interferogram
    H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
    [H1 , ~, ~] = fringe_estimate_coarse(H);
    if ii==first_F_IFG
        F_IFG(1:ii+step-1,jj:P) = H1(1:block-overlap, overlap+1:block);
    else
        F_IFG(ii:ii+step-1,jj:P) = H1(overlap+1:block-overlap, overlap+1:block);
    end
end


disp('Estimate main part...');
first_F_IFG=first_F_IFG+step;
for ii = first_F_IFG:step:IFG_L_stop 				%filtered interferogram
    for jj = first_F_IFG:step:IFG_P_stop
        H = IFG(ii-overlap:ii + (block-overlap) - 1, jj-overlap:jj + (block - overlap) - 1);  % cuting IFG
        [H1 , ~, ~] = fringe_estimate_coarse(H);
        F_IFG(ii:ii+step-1,jj:jj+step-1) = H1(overlap+1:block-overlap, overlap+1:block-overlap);
        if (ii*jj) > all_step * p;
            disp(['progress: ', num2str(10*p),'%']);
            p = p+1;
        end;
    end
end

disp('progress: 100%');
disp('End Process! Have a nice day. Bye;-)');
%------------------------------------------------------------------------
%END=====================================================================



