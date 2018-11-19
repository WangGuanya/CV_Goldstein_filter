%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Guanya Wang %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% wangguanya@csu.edu.cn %%%%%%%%%%%%%%
%%%%%%%%%%%%% 2017/11/07 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Revise: 2018/10/01 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Reference hardware conditions %%%%%%
%%%%%%%%%%%%% 4 cores of CPU i7-4790 3.6GHZ %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% ===============================================================
% =============================== INPUT ===================================
% =========================================================================

% ======================= Input Interferogram =============================
% =========================================================================
% Load existing data 
sim_phase_with_noise_2 = load('sim_phase_with_noise_two_look.dat');
% Euler Trasformation 
ifg_2 =  exp(1j*sim_phase_with_noise_2);
% Change Fringe Characters
IFG_2 = sim_phase(ifg_2);

% ======================= Input Coherence_map =============================
% =========================================================================
coherence_map = load('Data_sim/Coherence_map.dat');
cc_map = coherence_map;

% ======================== Set Key parameters =============================
% =========================================================================
block=32;
step=4;
L_2 = 2;



%% ===============================================================
% ======================= Fringe Reconstruction ===========================
% =========================================================================
fringe_2 = fringe_reconstruction_corase(IFG_2, 12, 5);
% figure; imagesc(angle(fringe_2)); axis image; colorbar; colormap jet;
left_2 = IFG_2.*conj(fringe_2);

%****************************** Display ***********************************
figure; imagesc(angle(left_2)); axis image; set(gca,'FontName','Times New Roman','FontSize',10);
colorbar; colormap jet; title('Interferogram without fringe'); 
saveas(gcf,'Left_2.tif');



%% ===============================================================
% ============================= Filter ====================================
% =========================================================================
[cvFiltered_IFG_2, cvB2, cvCV2] = cv_goldstein_sim( IFG_2, left_2, cc_map, block, step );



%% ===============================================================
% ===================== Image quality assessment ==========================
% =========================================================================
[ cvPSD2, cvresidue2, cvSPD2 ]= index_model ( cvFiltered_IFG_2);
residue = cvresidue2;
psd = mean2(cvPSD2);



%% ===============================================================
% ====================== Display ==========================================
% =========================================================================
% ====== filtered interferogram =============
figure;
subplot(1,2,1); imagesc(angle(IFG_2)); axis image; colormap jet; title('Original interferogram');
subplot(1,2,2); imagesc(angle(cvFiltered_IFG_2)); axis image; colormap jet;title('Filtered interferogram');
set(gca,'XTick',[]); set(gca,'XTickLabel',[]); 
set(gca,'YTick',[]); set(gca,'YTickLabel',[]);
saveas(gcf,'Result.tif');

% ====== cv_alpha_fit plot =============
fit_2 = [];
x = 0.6:0.01:4;
fit_2 = (cvB2(1)./(x+cvB2(3))) + cvB2(2);  
figure;
plot(x,fit_2,'-.','Color',[0.9 0.17 0.17],'linewidth',2);
legend('L=2');
axis([0.5 4 0 1]);
xlabel('CV'); ylabel('¦Á');
set(gca,'FontName','Times New Roman','FontSize',10);
saveas(gcf,'cv_alpha.tif');






