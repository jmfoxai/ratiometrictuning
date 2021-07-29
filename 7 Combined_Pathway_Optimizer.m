%FAS Optimizer Function:
%Optimizes model parameters toward matching model results (least squares) to data or
%specified outputs. Utilizes multi-parameter approaches (constrained
%or unconstrained)
%   Input:
%       p_vec0: initial guess for the values to be scaled in p_vec
%   Output:
%       p_vec_sol: optimized scaling values found by fminsearch (either
%       until objective tolerance was met or the number of function
%       evaluations or number of iterations exceeded the maximum value,
%       specified below). 
function [p_vec_sol] = Combined_Pathway_Optimizer(p_vec0)


fun_evals_max = 1000;%max number of function evaluations
max_iter = 180;%max number of optimiztion iterations

%set up optimization call

%Replace values to be scaled with components of vector x.
fitfunc = @(x) Combined_Pathway_Handler([142473.7238 7597.676912 4.276689943 40213.92919 88.88525384 0.005388274 4.645634978 0.006677519 0.284982219 -0.285700283 3.348915642 2.886607673 x(1) 2180.050007 x(2) x(3) x(4) x(5)]);

options = optimset('MaxFunEvals',fun_evals_max,'Display','iter','MaxIter',max_iter);
[p_vec_sol,~,~,~] = fminsearch(fitfunc,p_vec0,options);
