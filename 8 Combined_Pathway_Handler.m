%FAS Handler Function
%Calculates objectives to be returned to Combined_Pathway_Optimizer.m for
%optimization of scaling parameters to experimental data. 
%   Input:
%       p_vec: 18 value vector containing scaling parameters for the model.
%   Output:
%       total_obj: the total objective to be returned to the Optimizer,
%       which uses fminsearch to lower this objective by chaning the
%       variable values in p_vec. 

%All units are uM and seconds.

function [total_obj] = Combined_Pathway_Handler(p_vec)


% if p_vec(2) < 100
%     total_obj = 1E8;
%     return
% end


param_names = {{'k1_1f';'k1_1r';'k1_2f';'k1_2r';'kcat1_1';'k1_3f';'k1_3r';...
'kcat1_2';'k2_1f';'k2_1r';'k2_2f';'k2_2r';'k2_3f';'k2_3r';'k2_4f';...
'k2_4r';'k3_1f';'k3_1r';'k3_2f';'k3_2r';'k3_3f';'k3_3r';'k3_4f';'k3_4r';'k3_5f';'k3_5r';'kcat3';...
'k10_4f';'k10_4r';'k10_5f';'k10_5r';'k10_6f';'k10_6r';'kcat10_H';...
'k8_4f';'k8_4r';'k8_5f';'k8_5r';'k8_6f';'k8_6r';'kcat8_H';...
'k10_7f';'k10_7r';'k10_8f';'k10_8r';'k10_9f';'k10_9r';'kcat10_CO2';...
'k8_7f';'k8_7r';'k8_8f';'k8_8r';'k8_9f';'k8_9r';'kcat8_CO2';...
'k4_1f';'k4_1r';'k4_2f';'k4_2r';'kcat4';'k5_1f';'k5_1r';...
'kcat5';'k5_2r';'k6_1f';'k6_1r';'k6_2f';'k6_2r';'kcat6';'k7_1f';'k7_1r';...
'kcat7';'k8_1f';'k8_1r';'k8_2f';'k8_2r';'k8_3f';'k8_3r';'kcat8';...
'k9_1f';'k9_1r';'kcat9';'k9_2r';'k10_1f';'k10_1r';'k10_2f';'k10_2r';'k10_3f';'k10_3r';'kcat10';...
}};


%set of parameter values that are varied (or likely subject to change) in
%structure sent to parameterization function param_func.m
acp_bind = p_vec(12);%parameter "e"

field5 = 'opt_name'; val5 = {{'k2_4f','k3_4f','k3_4r','k3_5f','k3_5r','k7_1f','kcat7','k8_1f','k8_2f','k10_1f','k10_2f','kcat5','kcat8','kcat9','kcat10','k5_2r','k9_2r','k5_1f','k9_1f'}};%parameters with unique modifications in param_func
field6 = 'param_names'; val6 = param_names;
field10 = 'scaling_factor_init'; val10 = p_vec(1);%parameter "a1"
field11 = 'scaling_factor_elon'; val11 = p_vec(2);%parameter "a2"
field12 = 'scaling_factor_kcat'; val12 = p_vec(5);%parameter "c2"
field13 = 'scaling_factor_term'; val13 = p_vec(3);%parameter "a3"
field14 = 'scaling_factor_fabf'; val14 = p_vec(2);%parameter "a2" (option here to modify FabF scaling seperately)
field15 = 'scaling_factor_kcat_term'; val15 = p_vec(6);%parameter "c3"
%ACP_inh lists kon,koff (in this order) for inhibitory ACP binding
%[FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB]
field16 = 'ACP_inh'; val16 = [acp_bind*2.41E-04,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,acp_bind*2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02,2.41E-03,4.22*2.17E-02];
%field17 = 'enzyme_conc'; val17 = enz_conc;%Initial enzyme conentrations (ACC,FabD,FabH,FabG,FabZ,FabI,TesA,FabF,FabA,FabB)
field18 = 'inhibition_kds';val18 = (1/acp_bind).*[4335.05,23.67;824.7,30.1;967.5,7.55;251.52,8.484;128.78,2.509];%(Acyl-ACP binding FabH Kd values in pairs (binding to FabH and FabH*), first two values are Kd for 4-12, subsequent values are 14-20)
field19 = 'inhibition_on_rates';val19 = [0.3088,1.552];%On rates of acyl-ACP binding FabH (binding to FabH and FabH*)
field20 = 'kd_fits';val20 = [p_vec(8),p_vec(9),p_vec(4),p_vec(9)]; %[b2,b3,b1,b3] %Keq or Kd values that are fit
field21 = 'lin_param';val21=[p_vec(10) p_vec(11)];% [d1 d2]; TesA linear free energy slope and intercept (as a function of chain length)
field22 = 'scaling_factor_kcat_init';val22 = p_vec(7);%parameter "c1"
field23 = 'scaling_factor_FabA_unsat';val23 = p_vec(13);%parameter "f"
field24 = 'scaling_factor_FabAZ_kcat';val24 = p_vec(14);%parameter "c4"
field25 = 'TesA_fitting_source';val25 = 'Pf';%source of thioesterase fitting ('Pf'; Pfleger group measurements, 'Fox'; Fox group measurements', 'Non-native'; For alternative thioesterases)
field26 = 'scaling_factor_kcat8_CO2'; val26 = p_vec(15); %scaling for FabF kcatCO2
field27 = 'scaling_factor_kcat10_CO2'; val27 = p_vec(16); %scaling for FabB kcatCO2
field28 = 'scaling_factor_aCoA_8'; val28 = p_vec(17); %scaling for FabF binding with aCoA in decarboxylation
field29 = 'scaling_factor_aCoA_10'; val29 = p_vec(18); %scaling for FabB binding with aCoA in FabH-style


opt_struct = struct(field5,val5,field6,val6,...
    field10,val10,field11,val11,field12,val12,field13,val13,field14,val14,field15,val15,...
    field16,val16,field18,val18,field19,val19,field20,val20,field21,val21,field22,val22,field23,val23,field24,val24,field25,val25,...
    field26,val26,field27,val27,field28,val28,field29,val29);

%this loads the parameter data in test_data file in a tables stucture
tables = struct;
tables.('param_table') = readtable('est_param.csv','ReadRowNames',true);
tables.('km_table') = readtable('km_est.csv','ReadRowNames',true);

%generates the parameter values in paramter structure from the parameter
%function
param_struct = struct;
for i = 1:length(param_names{1})
    param_struct.(param_names{1}{i}) = param_func(param_names{1}{i},i,opt_struct,tables);
end

% Various checks to make sure that parameters aren't becoming unreasonable.
tp_vec = p_vec;
tp_vec(10) = abs(tp_vec(10));
tp_vec(11) = abs(tp_vec(11));
if any(tp_vec(tp_vec<=0))
    total_obj = 1E8;
    return
elseif tp_vec(4) > 6.29E4 || tp_vec(5) > 240 || tp_vec(6) > 1 || tp_vec(7) > 15
    total_obj = 1E8;
    return
elseif tp_vec(1) < .001 || tp_vec(2) < .001 || tp_vec(3) < .001 || tp_vec(5) < 0.01 || tp_vec(6) < 7.41E-05 || tp_vec(7) < 0.001 || tp_vec(8) < 2.46E-7 || tp_vec(9) < 2.46E-7
    total_obj = 1E8;
    return
end

for i = 1:length(param_names{1})
    if ismember(param_names{1}{i},{'k2_1f','k2_3f','k3_1f','k3_3f','k4_1f','k4_2f','k5_1f','k6_1f','k6_2f','k7_1f','k8_1f','k8_3f'})
        if ismember(param_names{1}{i},{'k2_1f','k2_3f','k3_1f'})
            if param_struct.(param_names{1}{i}) > 1650
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k3_3f','k4_2f','k5_1f','k6_2f','k8_1f','k8_3f'})
            if param_struct.(param_names{1}{i}) > 629
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k4_1f','k6_1f'})
            if param_struct.(param_names{1}{i}) > 1650
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            elseif param_struct.(param_names{1}{i+1})/param_struct.(param_names{1}{i}) < .01
                total_obj = 1E8;
                return
            end
        elseif ismember(param_names{1}{i},{'k7_1f'})
            if max(param_struct.(param_names{1}{i})) > 629
                total_obj = 1E8;%param_struct.(param_names{1}{i});
                return
            end
        end
    end
end

tic
%load('JpMat.mat','JpMatPrime')

%% Calculation of Obj1 (New data with various combinations of FabH, FabH, FabB, AcCoA removed)
%Add enz_conc vector, initial conditions, time, and make a function call to
%the solver. The resulting concentrations from the solver can be used for calculating the objective here.  
range_init = [0 150];

init_cond_1 = zeros(315,1);
init_cond_1(3) = 500;%s3 (Acetyl-CoA, 0.5mM)
init_cond_1(4) = 10;%s6 (holo ACP)
init_cond_1(5) = 1000;%s7 (NADPH)
init_cond_1(6) = 1000;%s8 (NADH)
init_cond_1(8) = 500;%p2 (malonyl-CoA)

enz_conc_ref = [0 1 1 1 1 1 10 1 1 1];
rate(1) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_ref,range_init,init_cond_1);
enz_conc_H = [0 1 0 1 1 1 10 1 1 1];
rate(2) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_H,range_init,init_cond_1);
enz_conc_HF = [0 1 0 1 1 1 10 0 1 1];
rate(3) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_HF,range_init,init_cond_1);
enz_conc_HB = [0 1 0 1 1 1 10 1 1 0];
rate(4) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_HB,range_init,init_cond_1);

init_cond_2 = zeros(315,1);
init_cond_2(3) = 0;%s3 (Acetyl-CoA, 0.5mM)
init_cond_2(4) = 10;%s6 (holo ACP)
init_cond_2(5) = 1000;%s7 (NADPH)
init_cond_2(6) = 1000;%s8 (NADH)
init_cond_2(8) = 500;%p2 (malonyl-CoA)

rate(5) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_H,range_init,init_cond_2);
rate(6) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_HF,range_init,init_cond_2);
rate(7) = Combined_Pathway_Solver_init(p_vec,param_struct,enz_conc_HB,range_init,init_cond_2);

rate_exp = [4.655129552 2.41772584 1.796165966 2.976304272 2.91667542 0.893110994 2.428676471];

obj1 = sum((rate_exp - rate).^2);
%% obj2
range_init_2 = [0 720];
obj23 = Combined_Pathway_Solver(p_vec,param_struct,enz_conc_ref,range_init_2,init_cond_1); 

%%
total_obj = obj1*obj23(1)*obj23(2);
end

