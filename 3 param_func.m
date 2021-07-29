function param = param_func(name,index,opt_struct,tables)
%param_func Outputs kinetic parameters given the pathway step and specified
%options (including fit parameters)
%Returns the kinetic parameters for each label in param_names
%The kinetic parameters are specified as a single
%value for the initiation steps, and as a vector for the elongation and
%termination steps.
%   Input:
%       name: the label (a string) being passed
%       index: the index in the list of the label
%       opt_struct: structure containing inputs that determine
%       parameterization
%       tables: tables containing data on initial parameter values
%   Output:
%       param: A vector or double (depending on the label being passed) of
%       kinetic parameters. Enzymes in elongation steps need a vector as a
%       different paramter value for each subsequent elongation step is
%       an important possibility


opt_name = opt_struct.('opt_name');
param_names = opt_struct.('param_names');
scaling_factor_init = opt_struct.('scaling_factor_init');%parameter "a1"
scaling_factor_elon = opt_struct.('scaling_factor_elon');%parameter "a2"
scaling_factor_term = opt_struct.('scaling_factor_term');%parameter "a3"
scaling_factor_fabf = opt_struct.('scaling_factor_fabf');%parameter "a2"
scaling_factor_kcat_term = opt_struct.('scaling_factor_kcat_term');%parameter "c3"
scaling_factor_kcat = opt_struct.('scaling_factor_kcat');%parameter "c2"
scaling_factor_kcat_init = opt_struct.('scaling_factor_kcat_init');%parameter "c1"
scaling_factor_FabAZ_kcat = opt_struct.('scaling_factor_FabAZ_kcat');%parameter "c4"
inhibition_kds = opt_struct.('inhibition_kds');
inhibition_on_rates = opt_struct.('inhibition_on_rates');
kd_fits = opt_struct.('kd_fits');%parameter [b2,b3,b1,b3] 
lin_param = opt_struct.('lin_param');
TesA_source = opt_struct.('TesA_fitting_source');
scaling_factor_kcat8_CO2 = opt_struct.('scaling_factor_kcat8_CO2');
scaling_factor_kcat10_CO2 = opt_struct.('scaling_factor_kcat10_CO2');
scaling_factor_aCoA_8 = opt_struct.('scaling_factor_aCoA_8');
scaling_factor_aCoA_10 = opt_struct.('scaling_factor_aCoA_10');

lin_slope = lin_param(1);%parameter "d1"
lin_int = lin_param(2);%parameter "d2"




%this loads the parameter data from the tables structure
param_table = tables.('param_table');
km_table = tables.('km_table');



%searches for the index in which the first elongation step occurs (this is
%why param_names needs to have all initiation steps listed first)
for i = 1:length(param_names)
    param_val = char(param_names(i));
    first_char = param_val(1:2);
    if first_char == char('k4')
        split_index = i;
        break
    end
end


%creates structure to hold the data from the param_table, with the table
%variable names assigned as fields
param_struct = struct;
for i = 1:length(param_names)
    %if the parameter is kcat, then scale each estimated kcat value from
    %the table using the appropriate kcat scaling terms
    if strcmp('kcat',param_names{i}(1:4))
        if strcmp('kcat7',param_names{i})%termination scaling for TesA
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat_term;%termination scaling for TesA
        elseif ismember(param_names{i},{'kcat4','kcat6','kcat8','kcat9','kcat10'})%elongation scaling for FabG,FabI,FabF,FabA,FabB
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat;
        elseif ismember(param_names{i},{'kcat5'})%FabZ (c4) scaling
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_FabAZ_kcat;
        elseif ismember(param_names{i},{'kcat10_H','kcat8_H'}) %new initiation by decarboxylation
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}; % If we don't want to scale kcat_H
%             if ismember(param_names{i},{'kcat10_H'}) %scale for FabB separately
%                 param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat10_H;
%             else
%                 param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat8_H;
%             end
        elseif ismember(param_names{i},{'kcat10_CO2','kcat8_CO2'}) %new initiation by decarboxylation
            %param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcatCO2;
            if ismember(param_names{i},{'kcat10_CO2'}) %scale for FabB separately
                 param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat10_CO2;
            else %scale for FabF separately
                param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat8_CO2;
            end 
        else %initiation scaling for FabD, FabH
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat_init;
        end
        
        %if the parameter is kon or koff, then scale both using the
        %appropriate scaling factors. Kon is further modified using the
        %estimated Kd values such that kd_est = koff_final/kon_final
        
        %For FabD,FabH use "a1" scaling
    elseif ismember(param_names{i},{'k2_1f','k2_1r','k2_3f','k2_3r','k3_1f','k3_1r','k3_3f','k3_3r',...
             'k10_4f','k10_4r','k10_6f','k10_6r','k8_4f','k8_4r','k8_6f','k8_6r',...
             'k10_7f','k10_7r','k10_8f','k10_8r','k8_7f','k8_7r','k8_8f','k8_8r'})
        if ismember(param_names{i},{'k2_1f','k2_1r','k2_3f','k2_3r','k3_1f','k3_1r','k3_3f','k3_3r'})
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_init;
%         elseif ismember(param_names{i},{'k10_8r'})
%             param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_aCoA_CO2;
        elseif ismember(param_names{i},{'k8_8r'})
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_aCoA_8;
        else
            param_struct.(param_names{i}) = param_table{{param_names{i}},:};
        end 
        %FabD and FabH are seperated here in case different scalings are
        %desired
        %{i+1} is used because  i+1 index is the reverse rate constant. kf = kr/kd
        if ismember(param_names{i},{'k2_1f','k2_3f'})
            kd_est = km_table{{param_names{i}},:};
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_init;
        elseif ismember(param_names{i},{'k3_1f','k3_3f'})
            kd_est = km_table{{param_names{i}},:};
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_init;
        %Add more elseif statements if you want to be able to adjust
        %various initiation parameters for FabB and FabF
        end
        
        
        %For TesA use parameter "a3"
    elseif strcmp('k7_1f',param_names{i}) || strcmp('k7_1r',param_names{i})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_term;
        
        %For FabF and FabB use paramter "a2" (seperated as scaling_factor_fabf in case
        %different scalings are desired, note scaling_factor_fabf=a2 here)
        %can also use scaling_vactor_elon = a2.
    elseif ismember(param_names{i},{'k8_1f','k8_1r','k10_1f','k10_1r'})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_fabf;
        if ismember(param_names{i},{'k8_1f','k10_1f'})
            kd_est = km_table{{param_names{i}},:};
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_fabf;
        end
        
        %This paramter modification step differs from other steps in that 
        %it is not a binding step, but a reverse catalytic step.
        %For FabZ and FabA reverse use parameter "c4" and "c2"
        %Note that FabZ and FabA are reversible reactions, the forward
        %reaction is denoted by kcat5 and kcat9, the reverse reaction by
        %k5_2r and k9_2r. As both forward and reverse are scaled by the
        %same constant the ratio (Keq) is maintained.
    elseif ismember(param_names{i},{'k5_2r','k9_2r'})
        if ismember(param_names{i},{'k5_2r'})
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_FabAZ_kcat;
        else
            param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_kcat;
        end
        
        %For FabD,FabH,FabF,FabB forward and reverse intermediate reaction
        %steps (the intermediate reaction of the ping-pong mechanism) use
        %scaling "b1"
    elseif ismember(param_names{i},{'k2_2f','k2_4f','k3_2f','k8_2f','k10_2f','k2_2r','k2_4r','k3_2r','k8_2r','k10_2r'})
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*kd_fits(3);   
    elseif ismember(param_names{i},{'k10_5f','k10_5r','k8_5f','k8_5r','k10_9f','k10_9r','k8_9f','k8_9r'})
         if ismember(param_names{i},{'k10_5r'})
             param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_aCoA_10;
         else
             param_struct.(param_names{i}) = param_table{{param_names{i}},:};
         end
  
        %For FabZ use parameter "a2"
    elseif ismember(param_names{i},{'k5_1f'})
        kd_est = km_table{{param_names{i}},:};
        param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_elon;
        
        %For FabG,FabI,FabF,FabA,FabB use parameter "a2"
    else
        param_struct.(param_names{i}) = param_table{{param_names{i}},:}*scaling_factor_elon;
        if ismember(param_names{i},{'k4_1f','k4_2f','k6_1f','k6_2f','k8_3f','k9_1f','k10_3f'})
            kd_est = km_table{{param_names{i}},:};
            param_struct.(param_names{i}) = (param_table{{param_names{i+1}},:}/kd_est)*scaling_factor_elon;
        end
    end
end


%if the parameter being looked for has index lower than split_index then it
%is an initiation step and a double will be returned from the structure
if index < split_index
    param = param_struct.(name);
end


num_elong_steps = 9;%number of elongation steps


%After the initial paramter assignment, additional parameters are
%assigned, and modified (for example incorporating substrate specificity)

%For FabH inhibition use appropriate on/off rates (note that inhibition has
%two types, noncompetitive (k3_4*) with respect to acetyl-CoA
%and competitive (k3_5*)with respect to malonyl-ACP

%noncompetitive inhibition with respect to acetyl-CoA (binding FabH)
%we would consider wanting to add FabB and FabF inhibition similar to that
%of FabH if we use FabB and FabF to inhibit. 
if strcmp(name,'k3_4f') || strcmp(name,'k3_4r')
    param = zeros(1,num_elong_steps);
    
    %off rate calculation from Kds (listed in inhibition_kds)
    if strcmp(name,'k3_4r')
        for i = 1:length(param)
            %inhibition for acyl-ACPs of 4-12 has the same value
            if i <= 5
                param(i) = inhibition_kds(1,1)*inhibition_on_rates(1);
            else
                param(i) = inhibition_kds(i-4,1)*inhibition_on_rates(1);
            end
        end
        
    %on rate calculation from Kds (same value for all chain lengths)
    else
        for i = 1:length(param)
            if i <= 5
                param(i) = inhibition_on_rates(1);
            else
                param(i) = inhibition_on_rates(1);
            end
        end
    end
    
%competitive inhibition with respect to malonyl-ACP (binding FabH*)
elseif strcmp(name,'k3_5f') || strcmp(name,'k3_5r')
    param = zeros(1,num_elong_steps);
    if strcmp(name,'k3_5r')
        for i = 1:length(param)
            %inhibition for acyl-ACPs of 4-12 has the same value
            if i <= 5
                param(i) = inhibition_kds(1,2)*inhibition_on_rates(2);
            else
                param(i) = inhibition_kds(i-4,2)*inhibition_on_rates(2);
            end
        end
    else
        for i = 1:length(param)
            %inhibition for acyl-ACPs of 4-12 has the same value
            if i <= 5
                param(i) = inhibition_on_rates(2);
            else
                param(i) = inhibition_on_rates(2);
            end
        end
    end
    
%Specify chain length dependence of kon and kcat for TesA
elseif strcmp(name,'k7_1f') || strcmp(name,'kcat7')
    param = zeros(1,num_elong_steps);
    
    %specify source of measurments, 'Pf'; Pfleger group measurements, 'Fox'; Fox group measurements
    if strcmp(name,'k7_1f')
        if strcmp(TesA_source,'Pf')
            kd_12 = exp(lin_slope*(12) + lin_int);%kd estimated from linear free energy relationship for chain lengths 12-20
            ratio_val = kd_12/(0.519*14.79);%ratio used to match kd at chain length 12
            kd_est = (ratio_val).*[473 293.9 52.986 14.79];%Kds for 4,6 and 8-12, estimated Kd is scaled to match linear free energy values
            for i = 1:num_elong_steps
                kd_long = exp(lin_slope*(i*2+2) + lin_int);
                if i<5
                    param(i) = param_struct.('k7_1r')/kd_est(i);%for 4-12 use linear free energy values
                else
                    param(i) = param_struct.('k7_1r')/kd_long;%for 12-20 use scaled estimates
                end
            end
        elseif strcmp(TesA_source,'Fox')
            kd_12 = exp(lin_slope*(12) + lin_int);
            ratio_val = kd_12/(1.0416*48);
            kd_est = (ratio_val).*[4815.2	4815.2	434	48];
            for i = 1:num_elong_steps
                kd_long = exp(lin_slope*(i*2+2) + lin_int);
                if i<5
                    param(i) = param_struct.('k7_1r')/kd_est(i);%for 4-12 use linear free energy values
                else
                    param(i) = param_struct.('k7_1r')/kd_long;%for 12-20 use scaled estimates
                end
            end
        elseif strcmp(TesA_source,'R3M4') %for TesA R3M4 from Pfleger
%             kd_12 = exp((-0.3978*(12)) + 4.8292);%kd estimated from linear free energy relationship for chain lengths 12-20
%             ratio_val = kd_12/(1.085541561*1.294673371);%ratio used to match kd at chain length 12
%             kd_est = (ratio_val).*[26.83111624 10.13381933 7.611809486 1.294673371];%Kds for 4,6 and 8-10, estimated Kd is scaled to match linear free energy values
            kd_est=[56.91208755 35.36250007 0.294555191 1.779555549 0.92358933 0.521582239 0.92358933 0.166345313 0.093940844];
            for i = 1:num_elong_steps
%                 kd_long = exp((-0.3978*(i*2+2)) + 4.8292);
%                 if i<5
                param(i) = param_struct.('k7_1r')/kd_est(i);%for 4-10 scaled estimates 
%                 else
%                     param(i) = param_struct.('k7_1r')/kd_long;%for 12-20 use scaled estimates use linear free energy values
%                 end
            end
        elseif strcmp(TesA_source,'R3M1') %for TesA R3M1 from Pfleger
            kd_est=[56.91208755 35.36250007 0.92358933 1.779555549 0.093940844 0.521582239 0.92358933 0.166345313 0.093940844];
            for i = 1:num_elong_steps
                param(i) = param_struct.('k7_1r')/kd_est(i);%for 4-10 scaled estimates 
            end
        elseif strcmp(TesA_source,'Non-native')
            kd_est = [1,1,1,1,1,1,1,1,1];
            for i = 1:num_elong_steps
                param(i) = param_struct.('k7_1r')/kd_est(i);%for all chain lengths
            end
        end

    %Implements TesA kcat chain length dependence as relative scaling (of
    %base value which is fit)
    elseif strcmp(name,'kcat7')
        if strcmp(TesA_source,'Pf')
            kcat_scaling = [0.0568,0.0509,0.1035,0.0158,0.25256,0.45819,1,1.221,1.5368];
        elseif strcmp(TesA_source,'Fox')
            kcat_scaling = [0.47895,0.47895,0.361111111,0.110185185,1.212962963,1.453703704,1.444444444,3.064814815,3.064814815];
        elseif strcmp(TesA_source,'R3M4') %for TesA R3M4
            kcat_scaling = [0.0568 0.0509 1 0.0158 0.25256 0.45819 0.25256 1.221 1.5368];
        elseif strcmp(TesA_source,'R3M1') %for TesA R3M1
            kcat_scaling = [0.0568 0.0509 0.25256 0.0158 1.5368 0.45819 0.25256 1.221 1.5368];
        elseif strcmp(TesA_source,'Non-native')
            kcat_scaling = [0,0,0,0,0,0,1,1,1];
        end
        for i = 1:num_elong_steps
                param(i) = param_struct.(name)*kcat_scaling(i);
        end
    end
    
%FabF and FabB on rates can be further modified by restricting elongation
%here by changing the value of i (i = 5 is chain length 12) and modifying
%the if/else statement
elseif ismember(name,{'k8_1f','k10_1f'})
    param = zeros(1,num_elong_steps);
    if strcmp(name,'k8_1f')
        for i = 1:num_elong_steps
            if i>5
                param(i) = param_struct.('k8_1f');
            else
                param(i) = param_struct.('k8_1f');%change if using restricted elongation mutant
            end
        end
    else
        for i = 1:num_elong_steps
            if i>5
                param(i) = param_struct.('k10_1f');
            else
                param(i) = param_struct.('k10_1f');%change if using restricted elongaiton mutant
            end
        end
    end
    
    
%Acyl transfer step k_fwd and k_rvs are determined by the fit parameters b2
%and b3, which are the Keq values for the transfer step
%b2 is Keq for FabD (first step) only
%b3 is Keq for FabD (second step), FabH, FabF, and FabB
elseif ismember(name,{'k2_2f','k3_2f','k8_2f','k10_2f'})
    if strcmp(name,'k2_2f') || strcmp(name,'k3_2f')
        
        %FabD first transfer step
        if strcmp(name,'k2_2f')
            param = param_struct.('k2_2r')/kd_fits(1);%b2
            
        %FabH transfer step
        elseif strcmp(name,'k3_2f')
            param = param_struct.('k3_2r')/kd_fits(2);%b3
        end
    else
        %FabF transfer step
        if strcmp(name,'k8_2f')
            for i = 1:num_elong_steps
                param(i) = param_struct.('k8_2r')/kd_fits(2);%b3
            end
        else
            %FabB transfer step
            for i = 1:num_elong_steps
                param(i) = param_struct.('k10_2r')/kd_fits(2);%b3
            end
        end
    end
    
    %FabD second transfer step
elseif strcmp(name,'k2_4f')
    param = param_struct.('k2_4r')/kd_fits(4);%b2
    
    
%chain length specificities for FabZ,FabF,FabA,FabB
elseif ismember(name,{'kcat5','kcat8','kcat9','kcat10'})
    
    %FabZ chain length kcat scaling (none)
    if strcmp(name,'kcat5')
        kcat_scaling_fabZ = [1,1,1,1,1,1,1,1,1];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabZ(i);
        end
        
    %FabF chain length kcat scaling
    elseif strcmp(name,'kcat8')
        kcat_scaling_fabF = [0.914,0.914,0.901,1,0.9,0.289,0.0222,0.0222,.0222];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabF(i);
        end
        
    %FabA chain length kcat scaling (none)
    elseif strcmp(name,'kcat9')
        kcat_scaling_fabA = [1,1,1,1,1,1,1,1,1];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabA(i);
        end
        
    %FabB chain length kcat scaling
    elseif strcmp(name,'kcat10')
        kcat_scaling_fabB = [0.855,0.855,0.975,0.967,1,0.125,0.0208,0.0208,.0208];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabB(i);
        end
    end
    
%For FabZ and FabA chain length specificities of k_rvs (reverse reaction
%rate 'k5_2r' and 'k9_2r') and kon
elseif ismember(name,{'k5_2r','k9_2r','k5_1f','k9_1f'})
    
    %FabZ chain length k_rvs scaling (none)
    if strcmp(name,'k5_2r')
        kcat_scaling_fabZ = [1,1,1,1,1,1,1,1,1];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabZ(i);
        end
    %FabA chain length k_rvs scaling (none)
    elseif strcmp(name,'k9_2r')
        kcat_scaling_fabA = [1,1,1,1,1,1,1,1,1];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kcat_scaling_fabA(i);
        end
    %FabZ chain length kon scaling
    elseif strcmp(name,'k5_1f')
        kon_scaling_fabZ = [0.469,1,0.296,0.372,0.2,0.0551,0.105,0.105,0.105];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kon_scaling_fabZ(i);
        end
    %FabA chain length kon scaling
    elseif strcmp(name,'k9_1f')
        kon_scaling_fabA = [0.0847,0.322,0.717,1,0.751,0.0847,0.0373,0.0373,0.0373];
        for i = 1:num_elong_steps
            param(i) = param_struct.(name)*kon_scaling_fabA(i);
        end
    end
end

%Can add more elseif statements if we want to be able to adjust more
%parameters separately.


%this returns the vector of values (all the same at this point) for kinetic
%parameters in elongation steps (as they should have indicies higher than
%split index)
if index >= split_index && ~any(strcmp(opt_name,name))
    param = zeros(1,num_elong_steps);
    for i = 1:num_elong_steps
        param(i) = param_struct.(name);
    end
end



end

