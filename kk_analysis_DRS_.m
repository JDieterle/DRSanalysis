% kk_analysis_DRS. Written by Johannes Dieterle
% description see end of the file.

function kk_analysis_DRS_(settings)
global pathname
global pathname_fit

try 
    run=evalin('base','run');
catch
    handle=figure;
assignin('base', 'handle', handle)
end
evalin('base', 'clear e2_fitted_e2c');
%##########################################################################
%modelparameters, fitparameters etc.:

%pathname for DRS raw data
% pathname=...
% 'O:\DRS\2014_08_07_pen1pic4_ch\';
[pathname_]=(pwd);
pathname=fullfile(pathname_, 'example');

%pathname for initial e2 parameters: if possible use *.xy.fitparms files
 pathname_initial_parameter=...
     fullfile(pathname);
filename_initial='initial';
filename_initial_end='.txt';

ema=struct(         'use_ema_yes_no',0,...
                    'thickness',0);
                
model_parm=struct(  'thickness',0,...
                    'substrate',1,...% 1: glas, 2: si native
                    'e1e2_substrate',0,...
                    'dOxide',0,...
                    'dstart', 25,...
                    'dstop',1,...
                    'ema', ema,...
                    'aoi', 0); % aoi in degree
                
                
initial_parms=struct('thickness', [],... % [] if no number, 0 for different initial parameters for each spectrum
                     'use', 0,...); % 0 means no everything else yes
                     'use_parms_of_previous_thickness', 'no');
                 
fit_parm=struct(    'extrapol', struct('E_k', []), ...
                    'lower_limit_no_absorption_eV_index',[1.5,0],...
                    'upper_limit_no_absorption_eV_index',[1.6,0],...
                    'upper_limit_small_absorption_eV_index',[1.65,0],...
                    'upper_limit_convex',[1.7,0],... % [] means: not used, not implemented yet
                    'cutspectrum_off_lower',1.5,...
                    'cutspectrum_off_upper',3, ...
                    'dEnergy', 0.005,...
                    'parameterization', 'ppp',...
                    'higher_energy_extrapolation', 'step',...
                    'initial_parms', initial_parms,...
                    'fit_ready', 0);

                    % weights for delta.... the larger the values the
                    % stronger the effect on:
                    fit_parm.weights.st=0.000005;     % short distance curvature in e2
                    fit_parm.weights.st_diff_c=0.0001; % difference between e2 
                                                  % at the end and e2 extraploated
                    fit_parm.weights.st_del=0.002;    % correlated distance and 
                                                  % difference between fit and experiment
                    fit_parm.weights.diff_exp_mod_transparent_range=1;
                    fit_parm.weights.e2_curv_transparent_range=10;
                    fit_parm.weights.dev_transparent_range='hhh';                
                                                  % deviation of e2 in

                    fit_parm.weight_e2_small_factor=00.7;       
                    
                    fit_parm.weights.bad_fit_in_transparent_range_select=2;
                    
%##########################################################################
% correct inputs:
fit_parm.cutspectrum_off_lower=max(fit_parm.cutspectrum_off_lower, ...
                                   fit_parm.lower_limit_no_absorption_eV_index);
                               
% take settings from GUI:
if nargin==1
    try
        pathname=settings.pathname_data;
    catch
        disp('no data file specified by GUI')
    end
    try
        initial_parms=struct('thickness', [],...
                     'use', ...
                     settings.initial_parms.use,...
                     'use_parms_of_previous_thickness',...
                     settings.initial_parms.use_parms_of_previous_thickness);
        fit_parm.initial_parms=initial_parms;
    catch
    end
    try
         [a1,a2,a3]=fileparts(settings.pathname_initial);
         pathname_initial_parameter=a1;
         filename_initial=a2;
         filename_initial_end=a3;
    catch
    end
    
    fit_parm.lower_limit_no_absorption_eV_index=str2num(settings.transparent_range.min);
    fit_parm.upper_limit_no_absorption_eV_index=str2num(settings.transparent_range.max);
    fit_parm.upper_limit_small_absorption_eV_index=str2num(settings.transparent_range.max_small_abs);
    
    fit_parm.cutspectrum_off_lower=str2num(settings.energy.min);
    fit_parm.cutspectrum_off_upper=str2num(settings.energy.max);
    fit_parm.dEnergy=str2num(settings.energy.d);
    model_parm.aoi=str2num(settings.aoi);

    model_parm.dstart=str2num(settings.thickness_range.min);
    model_parm.dstop=str2num(settings.thickness_range.max);

    model_parm.substrate=settings.substrate.substrate;
    model_parm.thickness=str2num(settings.substrate.d_oxide);    
    
    fit_parm.weights.st=str2num(settings.fit_parm.weights.st);
    fit_parm.weights.st_diff_c=str2num(settings.fit_parm.weights.st_diff_c);
    fit_parm.weights.st_del=str2num(settings.fit_parm.weights.st_del);   
    fit_parm.weights.diff_exp_mod_transparent_range=str2num(settings.fit_parm.weights.diff_exp_mod_transparent_range);
    fit_parm.weights.e2_curv_transparent_range=str2num(settings.fit_parm.weights.e2_curv_transparent_range);
    fit_parm.weight_e2_small_factor=str2num(settings.fit_parm.weight_e2_small_factor);
end
                               
%### create file for fit results containing time and date:
try 
    pathname_fit=evalin('base', 'pathname_fit');
    pathname_import=fullfile(pathname_fit, 'imported_initial_parms');
    pathname_extrapolation=fullfile(pathname_fit, 'extrapolation');
    pathname_weights=fullfile(pathname_fit, 'weights');
catch
    filenumber=clock;
    pathname_fit=fullfile(pathname, ['analysis', num2str(filenumber(1)-2000),'_',...
    num2str(filenumber(2)),'_',num2str(filenumber(3)),'_',...
    num2str(filenumber(4)),'_',num2str(filenumber(5))], '\');
    mkdir(pathname_fit)
    pathname_import=fullfile(pathname_fit, 'imported_initial_parms');
    mkdir(pathname_import)
    pathname_extrapolation=fullfile(pathname_fit, 'extrapolation');
    mkdir(pathname_extrapolation)
    assignin('base', 'pathname_fit', pathname_fit);
    pathname_weights=fullfile(pathname_fit, 'weights');
    mkdir(pathname_weights)
end

if nargin==1
    save(fullfile(pathname_fit, 'settings'), 'settings');
end
%### ...file created

global laufindex
laufindex=0;

% load DRS raw data:
data=prepare_data(model_parm, fit_parm);

% load, if necessary, dielectric function from substrate:
%model_parm.e1e2_substrate=load_substrate_optical_constants_xx(model_parm, data);
model_parm.ex_sub=load_ex_sub(model_parm, data);

%### prepare set e2=0 in transparent region:
% calculate index of energies with e2=0:
fit_parm.lower_limit_no_absorption_eV_index(2)=...
    energy_to_index(data.Energy, fit_parm.lower_limit_no_absorption_eV_index(1));
fit_parm.upper_limit_no_absorption_eV_index(2)=...
    energy_to_index(data.Energy, fit_parm.upper_limit_no_absorption_eV_index(1));
fit_parm.upper_limit_small_absorption_eV_index(2)=...
    energy_to_index(data.Energy, fit_parm.upper_limit_small_absorption_eV_index(1));
if isempty(fit_parm.upper_limit_convex)==0
fit_parm.upper_limit_convex(2)=...
    energy_to_index(data.Energy, fit_parm.upper_limit_convex(1));
end
% indexes of Energys with e2=0:
fit_parm.e2_eqal_zero(:,1)=...
    fit_parm.lower_limit_no_absorption_eV_index(2):1:fit_parm.upper_limit_no_absorption_eV_index(2);
%### ...set e2=0 ready.

% total number of e2-values:
fit_parm.num_e2_total=numel(data.Energy);

%### prepare extrapolation of e2:
if exist('pathname_extrapolation','var')==0
    pathname_extrapolation=pathname_fit;
end
  extrapolation_data=...
  prepare_extrapolation_old(data.Energy,pathname_extrapolation, fit_parm);

%put extrapolation data in fit_parm:
E_e2_extra=extrapolation_data.E_e2_extra;
e1_extra=extrapolation_data.e1_extra;
num_en_experiment=numel(data.Energy);
fit_parm.extrapol.num_summands_e2=size(E_e2_extra,2)-1;
fit_parm.extrapol.num_coeff_e2=fit_parm.extrapol.num_summands_e2;
fit_parm.extrapol.e1=e1_extra(1:num_en_experiment,:);
%### ...extrapolation preparations finished.

%### prepare weights for smoothening:
%prepare weigth data for smoothening:
Energy_range(:,1)=max(abs(extrapolation_data.lmax-...
                          extrapolation_data.lmin),1);
weight(:,1)=ones(numel(data.Energy),1);
weight(numel(data.Energy)+1:numel(data.Energy)+numel(extrapolation_data.lmin),1)=1./Energy_range;
weight(max(1,numel(data.Energy)-5*numel(extrapolation_data.lmin)):numel(data.Energy),1)=2;
weight(numel(data.Energy),1)=2;
weight(numel(data.Energy)+numel(extrapolation_data.lmin),1)=0;
weight(1:round(numel(data.Energy)/5))=1;
weight(1:fit_parm.upper_limit_small_absorption_eV_index(2))=1;
weight(1:fit_parm.upper_limit_no_absorption_eV_index(2))=1;

E_weight(:,1)=extrapolation_data.E_lower_limit:extrapolation_data.dEnergy:...
              2*extrapolation_data.E_upper_limit-...
              extrapolation_data.E_lower_limit;
E_weight(2:numel(weight)+1,2)=weight;
write_file(fullfile(pathname_weights,'E_weight_curv.dat'), E_weight, 'Energy\t weight_curv');
fit_parm.weight=weight;
clear extrapolation_data;

fit_parm.weight_e2_small(:,1)=1:-fit_parm.upper_limit_no_absorption_eV_index(2)+fit_parm.upper_limit_small_absorption_eV_index(2);
fit_parm.weight_e2_small(:,1)=sqrt(exp(-fit_parm.weight_e2_small(:,1)/5)*1000);
fit_parm.weight_e2_small(:,1)=max(fit_parm.weight_e2_small(:,1))-fit_parm.weight_e2_small(:,1);
fit_parm.weight_e2_small(:,1)=fit_parm.weight_e2_small(numel(fit_parm.weight_e2_small(:,1)):-1:1,1);
if max(fit_parm.weight_e2_small(:,1))>0
fit_parm.weight_e2_small(:,1)=fit_parm.weight_e2_small(:,1)/max(fit_parm.weight_e2_small(:,1))*fit_parm.weight_e2_small_factor;
end

%weights for datapoints:
fit_parm.weight_drs=ones(size(data.Energy));
%### ...weights finished.


% declare variables for KKT as visible in main function. i.e. persistent
E_ij_1=[];
E_ij_2=[];
energies_odd_number=[];
energies_even_number=[];
k_odd_number=[];
E_k_j_odd=[];
k_even_number=[];
E_k_j_even=[];
Delta_E=[];
e2x=[];
e10=[];
n=[];
k=[];
Energy_KKT=[];
Energy=data.Energy;



%##########################################################################
%###fit e2:
% analyse spectra:
if data.start>data.stop
    h=-1;
else
    h=1;
end

for k_=data.start:h:data.stop %running through the ML�s
    laufindex=0;
    fit_parm.fit_ready=0;

    % create folder for fit results of the monolayer:
    str=num2str((data.DRS_thickness(k_,1)+100), '%3.1f');
    fit_parm.pathname_fit_ml=strcat(pathname_fit,str(2:3), '_', str(5),'nm\');
    mkdir(fit_parm.pathname_fit_ml)
    mkdir(strcat(fit_parm.pathname_fit_ml,'fitparms\'))
        
    % put thickness in model_parm.thickness:
    model_parm.thickness=data.DRS_thickness(k_,1)

    num_points_without_absorption=numel(fit_parm.e2_eqal_zero(:,1));
    fit_parm.num_e2_fitted=numel(data.DRS_ML(:,k_))-num_points_without_absorption;
    num_summands_e2=fit_parm.extrapol.num_summands_e2;
    num_coeff_e2=fit_parm.extrapol.num_coeff_e2;
    DRS_exp(:,1)=data.DRS_ML(:,k_);
    assignin('base','DRS_ML_recent', DRS_exp)
    fit_parm.drs=DRS_exp;
    assignin('base','fit_parm', fit_parm);    
    
    % startparameter: e1_offset_parms=1, e2=0 �berall
    % Anzahl der Fitparameter: l ohne offset
    l=numel(DRS_exp)-num_points_without_absorption;
    fit_parm.num_e2_fitted=numel(DRS_exp)-num_points_without_absorption;
    num_parms=l+num_coeff_e2;
   
    %### load initial parameters:
    
    %- without file containing data
    if fit_parm.initial_parms.use==0
        e2_and_offset_initial(l+1-10:l+num_coeff_e2,1)=0.0;
             e2_and_offset_initial(l:l+num_coeff_e2,1)=0.0;    
             e2_and_offset_initial(l+8:l+num_coeff_e2,1)=0.5;    
             e2_and_offset_initial(round(3/4*l):l,1)=0.0;
             e2_and_offset_initial(round(7/8*l):l,1)=0.0;
        coeff_e2_initial=[];
        e2_initial=[];
        
    %- with file containing data   
    % preparation for taking initial parameters of previous spectra:
    else
        try
            e2_and_offset_initial=[];
            e2_fitted_e2c_all=[];
            e2_fitted_e2c_all2=[];
            fom_list1=[];
            fom_list2=[];
            e2_and_offset_initial(:,1)=evalin('base', 'e2_fitted_e2c');
            e2_fitted_e2c_all=evalin('base', 'e2_fitted_e2c_all');
            e2_fitted_e2c_all2=evalin('base', 'e2_fitted_e2c_all2');            
            fom_list1=evalin('base', 'fom_list1');            
            fom_list2=evalin('base', 'fom_list2');            
            disp('take starting parms from previous layer')
                        
            % select best initial parameters: 
            fom1=0.5*fom_list1(max(1,size(fom_list1,1)-2))+0.5*fom_list1(max(1,size(fom_list1,1)-2));
            fom2=0.5*fom_list2(max(1,size(fom_list2,1)-2))+0.5*fom_list2(max(1,size(fom_list2,1)-2));
            fom3=0.5*fom_list1(max(1,size(fom_list1,1)-1))+0.5*fom_list1(max(1,size(fom_list1,1)-1));
            fom4=0.5*fom_list2(max(1,size(fom_list2,1)-1))+0.5*fom_list2(max(1,size(fom_list2,1)-1));
            fom5=0.5*fom_list1(size(fom_list1,1))+0.5*fom_list1(size(fom_list1,1));
            fom6=0.5*fom_list2(size(fom_list2,1))+0.5*fom_list2(size(fom_list2,1));
            
            fom=[fom1*2.5 , fom2*2.5, fom3*2, fom4*2, fom5, fom6*1.3];
            [~,fmin]=min(fom);
            switch fmin
                case 1
                    e2_and_offset_initial=e2_fitted_e2c_all(:,max(1,size(fom_list1,1)-2));
                    disp('select 1. try of pre pre previous sprectrum')
                    fmin2_=0;
                case 2
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,max(1,size(fom_list2,1)-2));
                    disp('select 2. try of pre pre previous sprectrum')
                    fmin2_=2;
                case 3
                    e2_and_offset_initial=e2_fitted_e2c_all(:,max(1,size(fom_list1,1)-1));
                    disp('select 1. try of pre previous sprectrum')
                    fmin2_=1;
               case 4
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,max(1,size(fom_list2,1)-1));
                    disp('select 2. try of pre previous sprectrum')
                    fmin2_=4;
                case 5
                    e2_and_offset_initial=e2_fitted_e2c_all(:,size(fom_list1,1));
                    disp('select 1. try of previous sprectrum')
                    fmin2_=3;
                case 6
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,size(fom_list2,1));
                    disp('select 2. try of previous sprectrum')
                    fmin2_=6;
            end
            different_initial_parameters_present=1;
            
        catch
            e2_and_offset_initial=[];
            e2_fitted_e2c_all=[];
            e2_fitted_e2c_all2=[];
            fom_list1=[];
            fom_list2=[];
            assignin('base', 'e2_fitted_e2c_all', e2_fitted_e2c_all);
            assignin('base', 'e2_fitted_e2c_all2', e2_fitted_e2c_all2);
            assignin('base', 'fom_list1', fom_list1);
            assignin('base', 'fom_list2', fom_list2);
            disp('load initial parms from file...')
        end
        
        if fit_parm.initial_parms.use_parms_of_previous_thickness==0 ...
                || isempty(e2_and_offset_initial)
        
            % same initial parameters for all thicknesses:
            if isempty(fit_parm.initial_parms.thickness) || (fit_parm.initial_parms.thickness>0)
             if isempty(fit_parm.initial_parms.thickness)
                initial_parameter_file_name=...
                    fullfile(pathname_initial_parameter, [filename_initial,...
                    filename_initial_end]);
                initial_parameter_file_name2=...
                    fullfile(pathname_initial_parameter, [filename_initial,...
                    filename_initial_end]);
             else
                str=num2str((fit_parm.initial_parms.thickness+100), '%3.1f');
                initial_parameter_file_name=...
                    strcat(pathname_initial_parameter, filename_initial,...
                    str(2:3), '_', str(5), filename_initial_end);
                initial_parameter_file_name2=...
                    strcat(pathname_initial_parameter, filename_initial,...
                    num2str(fit_parm.initial_parms.thickness), filename_initial_end);
             end
             
             try
                imported_data=import_E_e1_e2_new(initial_parameter_file_name,fit_parm,data,pathname_import);
             catch           
             imported_data=import_E_e1_e2_new(initial_parameter_file_name2,fit_parm,data,pathname_import);
             end
             e2_initial_=imported_data.E_e1_e2(:,3);
             coeff_e2_initial=imported_data.e2extrapolation_coefficients;
            else
            % different initial parameters:
                str=num2str((data.DRS_thickness(k_,1)+100), '%3.1f');
                initial_parameter_file_name=...
                strcat(pathname_initial_parameter,'E_DR_e1_e2_e2c_DRex',...
                str(2:3), '_', str(5),'.dat');
                initial_parameter_file_name2=...
                strcat(pathname_initial_parameter,'E_DR_e1_e2_e2c_DRex',...
                num2str(round(data.DRS_thickness(k_,1))),'.dat'); 
                try
                    imported_data=import_E_e1_e2_new(initial_parameter_file_name,fit_parm,data,pathname_import);
                catch           
                    imported_data=import_E_e1_e2_new(initial_parameter_file_name2,fit_parm,data,pathname_import);
                end
                e2_initial_=imported_data.E_e1_e2(:,3);
                coeff_e2_initial=imported_data.e2extrapolation_coefficients;
            end

        assignin('base', 'filename_imp', initial_parameter_file_name);
        assignin('base', 'fit_parm_imp', fit_parm);
        assignin('base', 'data_imp', data);
        assignin('base', 'pathname_import', pathname_import);
    
        e2_initial(:,1)=cut_out_points_without_absorption(e2_initial_, fit_parm);
        e2_and_e2coefficients_initial(:,1)=[e2_initial; coeff_e2_initial];
        e2_and_offset_initial(:,1)=e2_and_e2coefficients_initial;
        end
    end
    %### initial parameters loaded.
    
    % limits for fitparameters:
    lb(:,1)=zeros(1,num_parms);
    lb(l+1:num_parms-3,1)=0.001;
    lb(l+1:num_parms-1,1)=0.0;
    ub(:,1)=ones(1,num_parms)*5;
    ub(l+1:num_parms,1)=10;
    
    k;
    
    % check bounds etc.
        %assignin('base', 'coeff_e2_initial', coeff_e2_initial)
        %assignin('base', 'e2_initial', e2_initial)
        assignin('base', 'ub', ub)
        assignin('base', 'lb', lb)
        assignin('base', 'DRS_exp_zeros', [DRS_exp;zeros(13,1)])
        assignin('base', 'Energy', data.Energy)
        assignin('base', 'e2_and_offset_initial', e2_and_offset_initial)
               
    if size(lb,1)~=size(ub,1)
        disp('lb ub not same length')
        s_lb=size(lb)
        s_ub=size(ub)
        error=xxx_error
    end
    if size(lb,1)~=size(e2_and_offset_initial,1)
        disp('bound and parms not same length')
        s_lb=size(lb)
        s_pa=size(e2_and_offset_initial)
        error=xxx_error
    end
    if size(data.Energy,1)~=size(DRS_exp,1)
        disp('size(Energy)~=size(DRS_exp)')
        s_E=size(data.Energy)
        s_DRS=size(DRS_exp)
        error=xxx_error
    end
    
    if sum(ub<e2_and_offset_initial)>0
        disp('upper bound lower than initial parameter!')
    end
    
e2_and_offset_initial=max(e2_and_offset_initial,lb);

    %### call fit algorithm:
        options = optimset('Algorithm','trust-region-reflective','DiffMinChange', 0.0001,...
        'ScaleProblem','jacobian','FinDiffType','central',...
        'Display', 'iter','MaxFunEvals', 50000);
    
        model_parm.R_ss_0=[];
        thickness=model_parm.thickness;
        model_parm.thickness=0;
        laufindex=laufindex-2;
        [~]=modelDRS_KKT_e1_e2_nested(e2_and_offset_initial,data.Energy, model_parm, fit_parm);
        model_parm.R_ss_0=evalin('base', 'R_ss_0');
        model_parm.thickness=thickness;
        
        E_DRS_0=modelDRS_KKT_e1_e2_nested(e2_and_offset_initial,data.Energy, model_parm, fit_parm);
        size_restraints=size(E_DRS_0,1)-size(data.Energy,1);

        e2_fitted_e2c=...
            lsqcurvefit(@(e2_and_offset,Energy)...
            modelDRS_KKT_e1_e2_nested...
            (e2_and_offset, Energy, model_parm, fit_parm),...
            e2_and_offset_initial, data.Energy, [DRS_exp*0;zeros(size_restraints,1)], ...
            lb,ub, options);

    assignin('base', 'e2_fitted_e2c', e2_fitted_e2c);   
    %### ...fit done.
    
    %### save fit result:
    fit_parm.fit_ready=1;        
    % e1 und e2 berechnen und mit E,DRSmodel in Matrix schreiben:
    E_DRS_e1_e2_e2c=modelDRS_KKT_e1_e2_nested(e2_fitted_e2c,data.Energy, model_parm, fit_parm);
    %save data:
    str=num2str((data.DRS_thickness(k_,1)+100), '%3.1f');
    write_file(fullfile(pathname_fit,['E_DR_e1_e2_e2c_DRex_', ...
        str(2:3), '_', str(5),'.dat']), E_DRS_e1_e2_e2c, 'Energy\t DRSmod\t e1\t e2\t extcoeff\t DRSexp\t dexpmod\t evale2')
    
    % save e2,  for higher energies:   
    save_e2_for_high_energies(fit_parm, data, pathname_fit, E_DRS_e1_e2_e2c,E_e2_extra, e1_extra, k_,''); 

    %### ...fit results (1st fit) saved. ##################################
    %######################################################################
    
    fit_parm.fit_ready=0;

    %######################################################################
    %### 2nd fit with other limits: #######################################
            options = optimset('Algorithm','trust-region-reflective','DiffMinChange', 0.001,...
        'ScaleProblem','jacobian','FinDiffType','central',...
        'Display', 'iter','MaxFunEvals', 5000);
    
    % limits for fitting parameters:
    lb(:,1)=zeros(1,num_parms);
    lb(l+1:num_parms-3,1)=0.0;
    lb(l+1:num_parms-1,1)=0.0;
    ub(:,1)=5*ones(1,num_parms);
    ub(l+1:num_parms,1)=10;
    
    d_exp_mod=sum(E_DRS_e1_e2_e2c(:,7).^2);
     fit_parm.weights.st=sum(E_DRS_e1_e2_e2c(1,8).^2)/d_exp_mod;
     fit_parm.weights.st_diff_=min(sum(E_DRS_e1_e2_e2c(2,9).^2)/d_exp_mod, d_exp_mod);
     fit_parm.weights.st_del=sum(E_DRS_e1_e2_e2c(3,10).^2)/d_exp_mod;
     fit_parm.weights.bad_fit_in_transparent_range=0.3;
    
    % e2=0 am unteren rand
    
    for m=1:numel(e2_fitted_e2c)
    factorr(m,1)=0.5+1.5/numel(e2_fitted_e2c)*m;
    end
    
    factorr(l+1:l+num_coeff_e2,1)=ones(num_coeff_e2,1);
    
    if exist('different_initial_parameters_present', 'var')
        % select best initial parameters: 
            fom1=0.5*fom_list1(max(1,size(fom_list1,1)-2))+0.5*fom_list1(max(1,size(fom_list1,1)-2));
            fom2=0.5*fom_list2(max(1,size(fom_list2,1)-2))+0.5*fom_list2(max(1,size(fom_list2,1)-2));
            fom3=0.5*fom_list1(max(1,size(fom_list1,1)-1))+0.5*fom_list1(max(1,size(fom_list1,1)-1));
            fom4=0.5*fom_list2(max(1,size(fom_list2,1)-1))+0.5*fom_list2(max(1,size(fom_list2,1)-1));
            fom5=0.5*fom_list1(size(fom_list1,1))+0.5*fom_list1(size(fom_list1,1));
            fom6=0.5*fom_list2(size(fom_list2,1))+0.5*fom_list2(size(fom_list2,1));
            
            fom=[fom1*3.5, fom2*3, fom3*2, fom4*3, fom5, fom6*2];
            if fmin2_~=0
            fom(fmin2_)=100; % exclude starting parameters of previous try
            end
            [~,fmin2]=min(fom);          
            switch fmin2
                case 1
                    e2_and_offset_initial=e2_fitted_e2c_all(:,max(1,size(fom_list1,1)-2));
                    disp('select 1. try of pre previous sprectrum')
                case 2
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,max(1,size(fom_list2,1)-2));
                    disp('select 2. try of pre pre previous sprectrum')         
                case 3
                    e2_and_offset_initial=e2_fitted_e2c_all(:,max(1,size(fom_list1,1)-1));
                    disp('select 1. try of previous sprectrum')
                case 4
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,max(1,size(fom_list2,1)-1));
                    disp('select 2. try of pre previous sprectrum')
                case 5
                    e2_and_offset_initial=e2_fitted_e2c_all(:,size(fom_list1,1));
                    disp('select 1. try of this sprectrum')
                case 6
                    e2_and_offset_initial=e2_fitted_e2c_all2(:,size(fom_list2,1));
                    disp('select 2. try of previous sprectrum')
            end
    end
    
    if data.DRS_thickness(k_,1)<5
    %e2= Calc_e2_thin([data.Energy, DRS_exp],model_parm.thickness);
    %e2_c=cut_out_points_without_absorption(e2, fit_parm);
    %e2_fitted_e2c(1:numel(e2_c))=e2_c;
    if data.DRS_thickness(k_,1)<5
    e2_fitted_e2c=e2_fitted_e2c*1.1;%.*factorr;
    disp('initial parameters =0 for d<3nm')
    lb=min(e2_fitted_e2c, lb);
    ub=max(e2_fitted_e2c, ub);
    end
    end
    
    index1=max(1,round(numel(e2_fitted_e2c)/4));
    e2_fitted_e2c_changed(1:numel(e2_fitted_e2c),1)=...
        [zeros(index1,1); e2_fitted_e2c(index1+1:numel(e2_fitted_e2c),1)];
    
    
    %e2_fitted_e2c_changed(l+1-3:l+num_coeff_e2-5,1)=0.15;
    %e2_fitted_e2c_changed(l+num_coeff_e2-5:num_parms,1)=1.6;
    %e2_fitted_e2c_changed(l+8:l+num_coeff_e2,1)=0.0; 
    
    E_DRS_0=modelDRS_KKT_e1_e2_nested(e2_fitted_e2c_changed,data.Energy, model_parm, fit_parm);
    size_restraints=size(E_DRS_0,1)-size(data.Energy,1);
    
    e2_fitted_e2c=...
    lsqcurvefit(@(e2_and_offset,Energy)...
    modelDRS_KKT_e1_e2_nested...
    (e2_and_offset, Energy, model_parm, fit_parm),...
    e2_fitted_e2c_changed, data.Energy, [DRS_exp*0;zeros(size_restraints,1)], ...
    lb,ub, options);

%################################

    %profile viewer
    %p = profile('info');
    %profsave(p,'profile_results')
    fit_parm.fit_ready=1; %
    
    % calculate e1 and e2 and out in matrix with E,DRStheo:
    E_DRS_e1_e2_e2c=modelDRS_KKT_e1_e2_nested(e2_fitted_e2c,data.Energy, model_parm, fit_parm);
    
    %save data. E,DRStheo,n,k,DRSexp...
    E_DRS_n_k_fit_DRSexp=[E_DRS_e1_e2_e2c];%, data.DRS_ML(:,k)];
    str=num2str((data.DRS_thickness(k_,1)+100), '%3.1f');
    %save (fullfile(pathname_fit,['E_DR_e1_e2_e2c_DRex', str(2:3), '_', str(5),'.dat']), 'E_DRS_n_k_fit_DRSexp', '-ascii');
  
    write_file(fullfile(pathname_fit,['E_DR_e1_e2_e2c_DRex_2_', ...
        str(2:3), '_', str(5),'.dat']), E_DRS_n_k_fit_DRSexp, 'Energy\t DRSexp\t e1\t e2\t extcoeff\t evale2\t DRSfit')
    
    
    
   % save e2,  for higher energies:   
    save_e2_for_high_energies(fit_parm, data, pathname_fit, E_DRS_e1_e2_e2c,E_e2_extra, e1_extra, k_,'_end');  

    
end
%### ...all fits done. ###################################################
%##########################################################################

%### save all results of last fit in one file:
    E_DRexp_all(:,1)=data.Energy;
    E_DRexp_mod(:,1)=data.Energy;
    E_e1_all(:,1)=data.Energy;
    E_e2_all(:,1)=data.Energy;
    E_e2_all2(:,1)=data.Energy;
for k_=data.start:h:data.stop %running through the ML�s
    str=num2str((data.DRS_thickness(k_,1)+100), '%3.1f');
    [~,data_k]=hdrload(strcat(pathname_fit,'E_DR_e1_e2_e2c_DRex_',...
        str(2:3), '_', str(5),'.dat'));
    E_DRexp_all(:, k_+1)=data_k(1:numel(data.Energy),6);
    E_DRmod_all(:, k_+1)=data_k(1:numel(data.Energy),2);
    E_e1_all(:,k_+1)=data_k(1:numel(data.Energy),3);
    E_e2_all(:,k_+1)=data_k(1:numel(data.Energy),4);
    
    [~,data_k]=hdrload(strcat(pathname_fit,'E_DR_e1_e2_e2c_DRex_2_',...
        str(2:3), '_', str(5),'.dat'));
%     E_DRexp_all(:, k+1)=data_k(:,7);
%     E_DRmod_all(:, k+1)=data_k(:,2);
%     E_e1_all(:,k+1)=data_k(:,3);
    E_e2_all2(:,k_+1)=data_k(1:numel(data.Energy),4);
    
end
    %save (strcat(pathname_fit,'E_DRexp_all', '.dat'), 'E_DRexp_all', '-ascii');
    write_file(strcat(pathname_fit,'E_DRexp_all', '.dat'), E_DRexp_all, 'Energy', data.DRS_thickness(1:data.start))
    %save (strcat(pathname_fit,'E_DRmod_all', '.dat'), 'E_DRmod_all', '-ascii');
    write_file(strcat(pathname_fit,'E_DRmod_all', '.dat'), E_DRmod_all, 'Energy', data.DRS_thickness(1:data.start))
    %save (strcat(pathname_fit,'E_e1_all', '.dat'), 'E_e1_all', '-ascii');
    write_file(strcat(pathname_fit,'E_e1_all', '.dat'), E_e1_all, 'Energy', data.DRS_thickness(1:data.start));    
    %save (strcat(pathname_fit,'E_e2_all', '.dat'), 'E_e2_all', '-ascii');
    write_file(strcat(pathname_fit,'E_e2_all', '.dat'), E_e2_all, 'Energy', data.DRS_thickness(1:data.start));
    %save (strcat(pathname_fit,'E_e2_all2', '.dat'), 'E_e2_all2', '-ascii');
    write_file(strcat(pathname_fit,'E_e2_all2', '.dat'), E_e2_all2, 'Energy', data.DRS_thickness(1:data.start));
%### ...data saved.
%##########################################################################
%### all fits finished and saved #############################################
%##########################################################################

%##########################################################################
% nested functions:
% modelDRS_KKT_e1_e2:
% transforms vector with fitparameters in e2 
% makes KKT
% calls ModelDRS_e1_e2
% appends additional goodness of the fit values to simulated DRS data
% saves fit results.

% DRS_ev=modelDRS_KKT_e1_e2(e2_and_offset, E, model_parm, fit_parm)
% e2_and_offset: vector containing e2 and e2 at higher energies
% E: vector Energy
% model_parm, fit_parm: structure containing parameters
% DRS_ev: if fit in progress: vector containing simulated DRS
%       : if fit finished: matrix containing DRS, e1, e2...  

function DRS_ev=modelDRS_KKT_e1_e2_nested(e2_and_offset, E, model_parm, fit_parm)

num_fit_parms=size(e2_and_offset,1);
%fit_parm.num_k_fitted
%fit_parm.num_k_total
%fit_parm.extrapol.num_summands_k
%k(:,1)=zeros(fit_parm.num_k_total,1);


% ### transforms vector with fitparameters in e2 and extract coefficients...
% for transitions at higher energies: 
% insert 0 for energies with e2=0 in e2-Vector:
e2(:,1)=[e2_and_offset(1:fit_parm.num_e2_fitted,1);...
    zeros(num_points_without_absorption,1)];

for a=1:num_points_without_absorption
    %translate all e2 with index > index e2=0 by +1
    e2(fit_parm.e2_eqal_zero(a,1)+1:fit_parm.num_e2_fitted+a)=...
        e2(fit_parm.e2_eqal_zero(a,1):fit_parm.num_e2_fitted+a-1);
    %set e2()=0
    e2(fit_parm.e2_eqal_zero(a,1),1)=0;
end

coeff_e2_fitted=abs(e2_and_offset(fit_parm.num_e2_fitted+1:...
    fit_parm.num_e2_fitted+fit_parm.extrapol.num_coeff_e2));
coeff_e2(:,1)=coeff_e2_fitted;
% ### ...e2 und coefficients for e2 at higher energies finished.


% ### KKT:
%##############################
%e2alt=e2;
e2_end=e2(numel(e2));%xxxx

e2x=[e2;e2_end*0.25];
e2x(numel(e2))=e2_end-coeff_e2(1)*0.25;

k=e2x;
Energy_KKT=[Energy;fit_parm.dEnergy+max(Energy)];
KKTxfast_;
% % e10=KKTxfast_([E;fit_parm.dEnergy+max(E)],e2x,1); %e1_offset=1
e10x=e10;
e10=n(1:numel(Energy)); % remove last element
% transitions at higher energies:
e1(:,1)=e10+sum(fit_parm.extrapol.e1*coeff_e2,2);
E_e1_e2=[E, e1, e2];
% ### ...KKT finished.

% calculate DRS
%DRS=modelDRS_e1_e2(E_e1_e2, model_parm);
if model_parm.substrate==1 % film, glas
num_layers=2;
d_=zeros(numel(E),num_layers);
d_(:,1)=model_parm.thickness*10^-9;
d_(:,2)=1000*10^-9;   
elseif model_parm.substrate==2 % film, sio, si
num_layers=3;
d_=zeros(numel(E),num_layers);
d_(:,1)=model_parm.thickness*10^-9;
d_(:,2)=model_parm.dOxide*10^-9;
d_(:,3)=1000*10^-9;
else
    disp('modelDRS_KKT_e1_e2: error no substrate specified')
end
ex_film=E_e1_e2(:,2)-1i*E_e1_e2(:,3);
ex_all=[ex_film,model_parm.ex_sub]; 

DRS=model_drs_data_all(E, ex_all,d_, model_parm.R_ss_0, model_parm.aoi);


% change weight of data points:
diff_DRSexp_mod=DRS-fit_parm.drs;

diff_DRSexp_mod=diff_DRSexp_mod.*fit_parm.weight_drs;

% add some additional evaluations
ev=evaluate_e2(E_e1_e2(:,3), coeff_e2,fit_parm.weight,fit_parm, DRS);

DRS_ev(:,1)=[diff_DRSexp_mod;ev.st; ev.diff;ev.drs_st; ev.e2small; ev.e2_transparent_range; ev.e2_curv_transparent_range];

% output, plotting etc...
laufindex=laufindex+1;


if mod(laufindex,num_fit_parms*6+1)==0 || laufindex==1
    try 
    run=evalin('base','run');
    catch
    run=1;   
    end
    if run==0
        xxx=stopped
    end
    
    iteration_=laufindex;
    sst=[' cur ', num2str(sum(ev.st), '%6.2s'), ' extr. c. diff.', num2str(sum(ev.diff), '%6.2s'),...
        ' curv. & dev. ', num2str(sum(ev.drs_st),'%6.2s'),...
        ' low abs.' num2str(sum(ev.e2small),'%6.2s'),...
        ' diff trans', num2str(sum(ev.e2_transparent_range), '%6.2s'),...
        ' diff trans c', num2str(sum(ev.e2_curv_transparent_range), '%6.2s')];
    
    disp([...' # eval. ', num2str(iteration_, '%10.2e')
        '       ', sst]);
    
    E_e1_e2_DRS=[E_e1_e2, DRS];
    %save (strcat(fit_parm.pathname_fit_ml,'E_e1_e2_DRS_iteration', ...
    %    num2str(laufindex),'.dat'), 'E_e1_e2_DRS', '-ascii');
    
    fitparms_xy=[E_e1_e2, [coeff_e2;zeros(numel(E)-numel(coeff_e2),1)], DRS(1:size(E_e1_e2,1))];
    %save (fullfile(fit_parm.pathname_fit_ml, 'fitparms', [num2str(laufindex),'.xy.fitparms']), ...
    %    'fitparms_xy', '-ascii');
    write_file(fullfile(fit_parm.pathname_fit_ml, 'fitparms', [num2str(laufindex),'.xy.fitparms']), ...
        fitparms_xy,...
        'Energy\t e1\t e2\t extcoeff\t DRSfit');

    try
    subplot(evalin('base', 'subplot_gII'));
    catch
    end
    
    positionVector1 = [0.1, 0.1, 0.3, 0.8];
    g=subplot('Position', positionVector1);
    assignin('base', 'subplot_gII', g);
    p = get(g,'position');
    p(4) = p(4)*0.8;
    p(2) = p(2)+p(4)*0.2;
    set(g, 'position', p);
    
    min_e1=max(0,(max(E_e1_e2(:,2))-1.5));
    plot(E_e1_e2(:,1), [E_e1_e2(:,2)-min_e1,E_e1_e2(:,3)])
    numstring=[num2str(min_e1),'0000'];
    legend(['e1xy-',numstring(1:3)],'e2xy')    
    axis([min(E_e1_e2(:,1)) max(E_e1_e2(:,1)) 0 1.5])
    xlabel('Energy [eV]','FontSize',11)
    ylabel(['\epsilon _{2},\epsilon _{1}-', numstring(1:3)],'FontSize',11)
    
    coeff=coeff_e2;
    coeff_index(:,1)=1:size(coeff,1);
    positionVector2 = [0.45, 0.1, 0.1, 0.8];
    %subplot(1,7,[3])
    g=subplot('Position', positionVector2);
        p = get(g,'position');
    p(4) = p(4)*0.8;
    p(2) = p(2)+p(4)*0.2;
    set(g, 'position', p);
    
    plot(coeff_index, coeff);
    axis([1 size(coeff,1) 0 1.5])
    xlabel('# coefficient','FontSize',11)
    
    positionVector3 = [0.65, 0.1, 0.3, 0.8];
    %subplot(1,7,[4,5])
    g=subplot('Position', positionVector3);
        p = get(g,'position');
    p(4) = p(4)*0.8;
    p(2) = p(2)+p(4)*0.2;
    set(g, 'position', p);
    drs_exp(:,1)=fit_parm.drs;
    plot(E_e1_e2(:,1), [DRS,drs_exp(:,1)])
    legend('DRS_{mod}','DRS_{exp}')
    xlabel('Energy [eV]','FontSize',11)
    ylabel('DRS','FontSize',11)
    title(['thickness ', num2str(model_parm.thickness), 'nm'])
    drawnow
    
    try 
    run=evalin('base','run');
    catch
    run=1;   
    end
    if run==0
        xxx=stopped
    end
end
    
if fit_parm.fit_ready==1 % fit finished, save all data, return results
    E_DRS_e1_e2_e2c(:,[1,3,4])=E_e1_e2(:,:);
    E_DRS_e1_e2_e2c(:,2)=DRS;    
    E_DRS_e1_e2_e2c(1:fit_parm.extrapol.num_summands_e2,5)=coeff_e2;
    drs_exp(:,1)=fit_parm.drs;
    E_DRS_e1_e2_e2c(:,6)=drs_exp;
    d_drs=(drs_exp-DRS);    
    E_DRS_e1_e2_e2c(:,7)=d_drs;
    diff_length=numel(ev.st)-numel(E_e1_e2(:,1))
    E_DRS_e1_e2_e2c(:,8)=ev.st(diff_length+1:numel(ev.st));
    E_DRS_e1_e2_e2c(1:numel(ev.diff),9)=ev.diff;
    E_DRS_e1_e2_e2c(1:numel(ev.drs_st),10)=ev.drs_st;
    
    DRS_ev=[];
    DRS_ev=E_DRS_e1_e2_e2c;
    
    % put in base workspace for initial parameter in following spectra:
    %fit_parm.upper_limit_no_absorption_eV_index(2)
    
    
    % compare difference between fit and exp in transparent range for
    % selection of initial parameters of following spectra
    index=fit_parm.upper_limit_no_absorption_eV_index(2);
    index2=max(index+10, round(numel(drs_exp)-(numel(drs_exp)-index)./3));
    difference_in_transparent_range=abs(mean((drs_exp(1:index)-DRS(1:index))));
    difference_in_rest_of_spectrum=abs(mean(drs_exp(index+1:index2)-DRS(index+1:index2)));
    bad_fit_in_transparent_range.q=difference_in_transparent_range/sqrt(min(difference_in_rest_of_spectrum,0.5))*...
                                   fit_parm.weights.bad_fit_in_transparent_range_select;
    bad_fit_in_transparent_range.a=sum(abs(drs_exp(index+1:index2)-DRS(index+1:index2)));
        
    b_q=bad_fit_in_transparent_range.q;
    b_a=bad_fit_in_transparent_range.a;
    if fit_parm.initial_parms.use==1
    e2_fitted_e2c_all=evalin('base', 'e2_fitted_e2c_all');
    e2_fitted_e2c_all2=evalin('base', 'e2_fitted_e2c_all2');
    fom_list1=evalin('base', 'fom_list1');
    fom_list2=evalin('base', 'fom_list2');
    if size(e2_fitted_e2c_all,2)>size(e2_fitted_e2c_all2,2)
    e2_fitted_e2c_all2(:,size(e2_fitted_e2c_all2,2)+1)=e2_and_offset;
    fom_list2=[fom_list2; [b_a, b_q]];
    assignin('base', 'e2_fitted_e2c_all2', e2_fitted_e2c_all2);
    assignin('base', 'fom_list2', fom_list2);
    else
    e2_fitted_e2c_all(:,size(e2_fitted_e2c_all,2)+1)=e2_and_offset;
    fom_list1=[fom_list1; [b_a, b_q]];
    assignin('base', 'e2_fitted_e2c_all', e2_fitted_e2c_all);
    assignin('base', 'fom_list1', fom_list1);
    end
    end
    
    save (fullfile(fit_parm.pathname_fit_ml,['E_DR_e1_e2_e2c', 'finished', ...
        '.dat']), 'E_DRS_e1_e2_e2c', '-ascii');
    disp('...analysis finished, data saved.')
    
    E_e2x=[[E;max(E)+fit_parm.dEnergy], e2x, [e10x;0]];
    save (fullfile(fit_parm.pathname_fit_ml,['E_e2x', ...
        '.dat']), 'E_e2x', '-ascii');
    
    fitparms_xy=[E_e1_e2, [coeff_e2;zeros(numel(E)-numel(coeff_e2),1)]];
    save (fullfile(fit_parm.pathname_fit_ml, 'fitparms', ['finished','.xy.fitparms']), 'fitparms_xy', '-ascii');
end



function KKTxfast_
n=zeros(numel(Energy_KKT),1);
if isempty(E_ij_1);
energies_odd_number(:,1)=Energy_KKT(1:2:numel(Energy_KKT),1);
energies_odd_number_quadrat(:,1)=energies_odd_number.^2;
energies_even_number(:,1)=Energy_KKT(2:2:numel(Energy_KKT),1);
energies_even_number_quadrat(:,1)=energies_even_number.^2;
E_ij_1=1./(repmat(energies_even_number_quadrat',numel(energies_odd_number),1)...
            -repmat(energies_odd_number_quadrat,1,numel(energies_even_number)));
E_ij_2=1./(repmat(energies_odd_number_quadrat',numel(energies_even_number),1)...
            -repmat(energies_even_number_quadrat,1,numel(energies_odd_number)));
end
Delta_E=(Energy_KKT(2)-Energy_KKT(1));

k_odd_number(:,1)=k(1:2:numel(k),1);
E_k_j_odd(:,1)=energies_odd_number.*k_odd_number;
k_even_number(:,1)=k(2:2:numel(k),1);
E_k_j_even(:,1)=energies_even_number.*k_even_number;

n(1:2:numel(Energy_KKT),1)=E_ij_1*E_k_j_even;
n(2:2:numel(Energy_KKT),1)=E_ij_2*E_k_j_odd;

n=1+4 / pi *Delta_E*n;

end

end


end

%load data: L�d Daten und bereitet DRS_ML (�quidistante Energie...), DRS_thickness
%vor
%parameter: model_parm
%Ergebnis:
%      struct: data.DRS_ML
%              data.DRS_thickness
%              ..............
function data=prepare_data(model_parm, fit_parm)
global pathname;
global pathname_fit;
% load experimental data (third row: thickness,the first 4 rows do not
% contain DRS data!):
disp(fullfile(pathname,'DRS_ML.txt'))
DRS_ML=load(fullfile(pathname,'DRS_ML.txt'));   
[row, col]=size(DRS_ML);

index=SelectThickness(DRS_ML, model_parm.dstart, model_parm.dstop);
stop=index.min
start=index.max

%extract thickness and cut first 4 rows:
DRS_thickness(:,1)=DRS_ML(3,:); %row vector
assignin('base', 'DRS_thickness', DRS_thickness)
DRS_MLorg=DRS_ML(5:row,:);

%check, whether first column is energy or wavelength:
if DRS_MLorg(10,1)<10
    Energy(:,1)=DRS_MLorg(:,1); % column vector
else
    Energy(:,1)=DRS_MLorg(:,2); % column vector
end

%%Make Energy have equal spacing and initialize DRS_ML matrix for interpolation later%%
Energyorg(:,1)=Energy(:,1); % column vector
E_max=max([Energyorg(1),Energyorg(numel(Energyorg))]);
E_min=min([Energyorg(1),Energyorg(numel(Energyorg))]);
Energy=[];
Energy(:,1)=E_max:-fit_parm.dEnergy(1):E_min;    
DRS_ML=zeros(numel(Energy),col);

%### Interpolate spectra of all monolayers for equal energy spacing:
if  start>stop
    h=-1;
else
    h=1;
end
for k=start:h:stop %running through the ML�s

    DRS_ML(:,k)=interp1(Energyorg,DRS_MLorg(:,k),Energy);
end
%### ...spectra interpolated.

%### make sure Energy increases with index and cut data outside limits:
if  Energy(1)>Energy(2)
    DRS_ML_reverse(size(DRS_ML,1):-1:1,:)=DRS_ML;
    Energy_reverse(size(Energy,1):-1:1,:)=Energy;
else
    DRS_ML_reverse=DRS_ML;
    Energy_reverse=Energy;    
end

Energy_reverse_cut_l=...
    cutoff_lower(Energy_reverse, ...
        Energy_reverse, fit_parm.cutspectrum_off_lower);
DRS_ML_reverse_cut_l=...
    cutoff_lower(Energy_reverse, ...
        DRS_ML_reverse, fit_parm.cutspectrum_off_lower);
Energy_reverse_cut=...
    cutoff_upper(Energy_reverse_cut_l, ...
        Energy_reverse_cut_l, fit_parm.cutspectrum_off_upper);
DRS_ML_reverse_cut=...
    cutoff_upper(Energy_reverse_cut_l, ...
        DRS_ML_reverse_cut_l, fit_parm.cutspectrum_off_upper); 
assignin('base','DRS_ML', DRS_ML_reverse_cut)
%### ...data ordered and.
    

data=struct('DRS_ML',DRS_ML_reverse_cut,...
            'Energy',Energy_reverse_cut,...
            'DRS_thickness',DRS_thickness,...
            'start',start,...
            'stop',stop);
        
% save DRS_ML and DRS_thickness:
save (fullfile(pathname_fit,'data_DRS_ML.txt'), 'DRS_ML_reverse', '-ascii');
save (fullfile(pathname_fit,'data_DRS_thickness.txt'), 'DRS_thickness', '-ascii');

clear DRS_ML_reverse Energy_reverse DRS_ML Energy DRS_ML_reverse_cut...
        DRS_ML_reverse_cut
end


% wandelt eine Energieposition in den Index der n�chsth�heren Energie im Energievektor um:
% Parameter:    Energy: vektor mit Energie
%               Energyposition: Position des Punktes im Spektrum (eV)
% Ergebnis:     index: Index des Energiewertes im Vektor Energy, der knapp
%               unterhalb von Energyposition liegt
function index=energy_to_index(Energy, Energyposition)
    a=1;
    while Energy(a)<Energyposition 
        % doesn�t work if energy is to big/small
        a=a+1;
    end
    index=a;
end

% remove all energies up to a given index
% parameter: energievektor, fit_parm.energieposition, daten
function cut=cutoff_lower(Energy, data, Energyposition)
    index=energy_to_index(Energy, Energyposition);
    cut=data(index:size(data,1),:);
end

% remove all energies above a given index
% parameter: energievektor, fit_parm.energieposition, daten
function cut=cutoff_upper(Energy, data, Energyposition)
    index=energy_to_index(Energy, Energyposition);
    cut=data(1:index,:);
end

% removes all elements of e2 vector for e2=0
function cut_out=cut_out_points_without_absorption(e2, fit_parm)
    e2_eqal_zero=fit_parm.e2_eqal_zero(:,1);
    num_points_without_absorption=numel(e2_eqal_zero);
    cut_out(:,1)=e2;
    for a=1:num_points_without_absorption
        cut_out(e2_eqal_zero(a)-a+1)=[];
    end
end

function ex_sub=load_ex_sub(model_parm, data)
Energy_=data.Energy;
if model_parm.substrate == 1 % film, glas
    h=6.6260755*10^(-34); %[Js]
    c=299792458; %[m/s]
    Joule_conversion=1.6*10^(-19);
    lambda_mu(:,1)=h*c./(Energy_(:,1)*Joule_conversion)*10^(6); %[�m]
    A_c=1.5147;
    B_c=0.0045386;
    C_c=-0.000022659;
    n_sub(:,1)=A_c+B_c./lambda_mu(:,1).^2+C_c./lambda_mu(:,1).^4; %real part
    ex_sub=n_sub.^2-1i*0;
end
if model_parm.substrate == 2 % film, sio2, si
    [~,si]=hdrload('Mat\si.dat');
    e1_e2_si(:,1)=interp1(si(:,1),si(:,2),Energy_);
    e1_e2_si(:,2)=interp1(si(:,1),si(:,3),Energy_);
    epsilon_si=e1_e2_si(:,1)-1i*e1_e2_si(:,2);

    [~,sio2]=hdrload('Mat\sio2.dat');
    e1_e2_sio2(:,1)=interp1(sio2(:,1),sio2(:,2),Energy_);
    e1_e2_sio2(:,2)=interp1(sio2(:,1),sio2(:,3),Energy_);
    epsilon_sio2=e1_e2_sio2(:,1)-1i*e1_e2_sio2(:,2);
    ex_sub=[epsilon_sio2,epsilon_si];
end
end

function e1e2_substrate=load_substrate_optical_constants_xx(model_parm, data)

% load, if necessary, dielectric function from substrate:
% in principle
% Si_jaw and SiO2_jaw are not necessary for fit on Au, but if one wants to
% extend the fitting to some additional semitransparent film, this would be
% the case
% pathname for substrate optical constants:
pathname2=...
    'R:\Wvase\Mat\Semicond\'; 
e1e2_substrate(:,1)=data.Energy;
if model_parm.substrate == 2 || model_parm.substrate == 2 || model_parm.substrate == 3
    [~,Si_jaw]=hdrload(fullfile(pathname2,'Si_jaw.mat'));
    [~,SiO2_jaw]=hdrload(fullfile(pathname2,'SiO2_jaw.mat'));
    e1e2_substrate(:,2)=interp1(Si_jaw(:,1),Si_jaw(:,2),data.Energy);
    e1e2_substrate(:,3)=interp1(Si_jaw(:,1),Si_jaw(:,3),data.Energy);
    e1e2_substrate(:,4)=interp1(SiO2_jaw(:,1),SiO2_jaw(:,2),data.Energy);
    e1e2_substrate(:,5)=interp1(SiO2_jaw(:,1),SiO2_jaw(:,3),data.Energy);
    if model_parm.substrate ==3
        load(fullfile(pathname2,'Au.txt'));
        e1e2_substrate(:,6)=interp1(Au(:,1),Au(:,2),data.Energy);
        e1e2_substrate(:,7)=interp1(Au(:,1),Au(:,3),data.Energy);
    end
elseif model_parm.substrate == 4
    load(fullfile(pathname2,'ito2.txt'));
    ito=ito2;
    e1e2_substrate(:,2)=interp1(ito(:,1),ito(:,2),data.Energy);
    e1e2_substrate(:,3)=interp1(ito(:,1),ito(:,3),data.Energy);
end
e1e2_substrate=e1e2_substrate(:,2:size(e1e2_substrate,2));
end

% save e2 for high energies
function save_e2_for_high_energies(fit_parm, data, pathname_fit, E_DRS_e1_e2_e2c, E_e2_extra, e1_extra, k_, filename)
    num_en_experiment=numel(data.Energy);
    % coefficients
    coeff_e2=E_DRS_e1_e2_e2c(1:fit_parm.extrapol.num_summands_e2,5);
    %-- koeffizienten*e2
    for i=1 : fit_parm.extrapol.num_summands_e2
        e2_extra_(:,i)=  [zeros(num_en_experiment,1);E_e2_extra(num_en_experiment+1:numel(E_e2_extra(:,1)),1+i)].* coeff_e2(i);
        e1_extra_(:,i)= e1_extra(:,i).* coeff_e2(i);
    end
    
    e2_extra_sum=e2_extra_;
    e2_extra_sum(:,size(e2_extra_,2)+1)=sum(e2_extra_,2);
    e1_extra_sum=e1_extra_;
    e1_extra_sum(:,size(e1_extra_,2)+1)=sum(e1_extra_,2);

    clear E_e2;
    clear E_e1;
    E_e2=[E_e2_extra(:,1), e2_extra_sum];
    E_e1=[E_e2_extra(:,1), e1_extra_sum];
    
  %-- save
  save (fullfile(pathname_fit,['E_e2_extra', filename,num2str(round(data.DRS_thickness(k_,1))),'.dat']), 'E_e2', '-ascii');
  save (fullfile(pathname_fit,['E_e1_extra', filename,num2str(round(data.DRS_thickness(k_,1))),'.dat']), 'E_e1', '-ascii');    
end


% kramers kronig constrained variational analysis of differential reflectance spectra (DRS) 
% see A. B. Kuzmenko, Rev. Sci. Instrum. 76, 083108 (2005); 
% http://dx.doi.org/10.1063/1.1979470
% and R. Forker, Annu. Rep. Prog. Chem., Sect. C: Phys. Chem., 2012,108, 34-68 
% DOI:  10.1039/C2PC90002E 
% In contrast to the two appraoches mentioned above this program is able to
% iteratively account for absorptions at energies above the measured range.

% kk_analysis_DRS is able to iteratively account for arbitrary absorptions 
% at higher energies. This is achieved by an extrapolation of e2
% parameterized in an optimized way to provide flexiblitily of the spectrum
% and a minmal number of parameters. In each function evaluation a full 
% kramers kronig transformation over a large energy range is performed
% which is possible very efficently in this concept. This procedure was developed 
% and implemented for the work reported in: Johannes Dieterle, Ph.D. thesis, Eberhard Karls
% Universit�t T�bingen, T�bingen (2016): The Influence of Dilution on
% Opitical and Structural Properties of Organic Semiconductors.

%##########################################################################
% needed files:
% DRS ML.txt:
% first column: Energy, Second: wavelength,
% other columns: DRS
% first 4 rows: header:
% third row: thickness
% first, second, 4. row: not important

% initial parameter file (if not *.fitparms):
% 1. col: Energy
% 2. col: e1
% 3. col: e2

%##########################################################################
% output files:
% pathname stored in pathname_fit

%E_DR_e1_e2_e2c_DRex_xx: file containing result for thickness xx: 
% 1. col Energy
% 2. col DRS_exp
% 3. col e1
% 4. col e2
% 5. col extrapolation coeff
% 6. col result of evaluate e2
% 7. col DRS_fit

% E_e1_extra_endxx for thickness xx:
% 1. col Energy
% other cols: extrapolated e1

% E_e2_extra_endxx for thickness xx:
% 1. col Energy
% other cols: extrapolated e2

% finished.xy.fitparms (in subfolder)
% initial parms for next fits
% 1. col Energy
% 2. col e1
% 3. col e2
% 4. col extrapolation coeff
%--------------------------------------------------------------------------
% E_e2_all.dat: e2 for all thicknesses
% first col: Energy
% other cols: e2
% first row: thickness

% E_DRexp_all, E_DRmod_all, E_e1_all : same as E_e2_all.dat for other
% quantities
%--------------------------------------------------------------------------

%##########################################################################