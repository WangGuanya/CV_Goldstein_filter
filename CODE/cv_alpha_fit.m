function [ Size_IFG, Size_CC, B, R, J, Cov, MSE ] = cv_alpha_fit( left, CC_map, block, step )
% Polynomial fit for CV vs.Coherence
%   [ B, R, J, Cov, MSE ] = cv_c_fit( IFG, x, y, L )
%
%   Guanya Wang, 2017/11/09
%
% ================= Input ====================
% IFG     :  Input interferogram
% CC_map  :  Input coherence map
% block   :  Filtering window, should be even a number, such as 32, 64, â€¦ï¿½?, default 32
% step    :  the step of sliding window
% 
%
%
% ================ Output ===================
% Size_IFG:  The size of input interferogram
% Size_CC :  The size of CC_map
% B       : Polynomial coefficient
% R       : Residual of poly_fit
% J       : Jacobian matrix
% Cov     : Estimated variance-covariance matrix
% Mse     : Mean Square Error in scalar value



[Size_IFG(1), Size_IFG(2)] = size(left);
[Size_CC(1), Size_CC(2)] = size(CC_map);

% ====== Calculate CV ======
[ CV, cv ] = cv_main( left, block, step );

% ====== Calculate Alpha ======
[ Alpha ] = alpha_main( CC_map, block, step )



% ====== Alpha vs. CV ======
for i = 1:length(Alpha)
    x_cv(i) = cv(i);
    y_alpha(i) = Alpha(i);
end
% figure; scatter(x_cv, y_alpha, 'r'); title('CV vs.Alpha'); 


% % ====== Initial polynomial_fit ======
% p = polyfit (x_cv, y_alpha, 2);
% y_fit = polyval(p, x_cv);
% hold on; plot(x_cv,y_fit);legend('Scatters','regress' );
% 
% % ====== Regression Analysis ======
% model = inline('B(1)*x.^2+B(2)*x+B(3)','B','x');             % definite fit function.
% [B, R, J, Cov, MSE]= nlinfit(x_cv, y_alpha, model, p);


% ====== Initial polynomial_fit ======
p = [1 0 0];

% ====== Regress polynomial_fit ========
 model = inline('a(1).*((x+a(3)).^(-1))+a(2)','a','x');             % definite fit function.
 [B, R, J, Cov, MSE]= nlinfit(x_cv, y_alpha, model, p);
 
 y_fit = B(1).*((x_cv+B(3)).^(-1)) + B(2);
 
 b1 = num2str(B(1));
 b2 = num2str(B(2));
 b3 = num2str(B(3));
 
% hold on; plot(x_cv, y_fit); legend('Scatters', 'y = b1/(x+b3)+b2');
% xlabel('CV '); ylabel('¦Á');
% set(get(gca,'XLabel'),'FontSize',15,'FontName','Times New Roman');
% set(get(gca,'YLabel'),'FontSize',15,'FontName','Times New Roman');
% set(gca,'FontName','Times New Roman','FontSize',10);

end


