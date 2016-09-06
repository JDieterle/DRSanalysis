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

function DRS_ev=modelDRS_KKT_e1_e2(e2_and_offset, E, model_parm, fit_parm)
num_points_without_absorption=numel(fit_parm.e2_eqal_zero(:,1));

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

e10=KKTx([E;fit_parm.dEnergy+max(E)],e2x,1); %e1_offset=1
e10x=e10;
e10(numel(E)+1)=[]; % remove last element
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

DRS=model_drs_data_all(E, ex_all,d_, model_parm.R_ss_0);


% data=[E_e1_e2];
% write_E_e1_e2_in_VASE_readable_file(data, 'R:\film.mat');
% data=[E_e1_e2(:,1),real(model_parm.ex_sub(:,1)), -imag(model_parm.ex_sub(:,1))];
% write_E_e1_e2_in_VASE_readable_file(data, 'R:\glas.mat');
% data=[E_e1_e2(:,1), DRS2, DRS];
% save('R:\drs.dat', 'data', '-ascii');
% diff_xx=sum(abs(DRS2-DRS))

% change weight of data points:
diff_DRSexp_mod=DRS-fit_parm.drs;

diff_DRSexp_mod=diff_DRSexp_mod.*fit_parm.weight_drs;

% add some additional evaluations
ev=evaluate_e2(E_e1_e2(:,3), coeff_e2,fit_parm.weight,fit_parm, DRS);

DRS_ev(:,1)=[diff_DRSexp_mod;ev.st; ev.diff;ev.drs_st; ev.e2small; ev.e2_transparent_range; ev.e2_curv_transparent_range];

% output, plotting etc...
global laufindex 
laufindex=laufindex+1;


if mod(laufindex,num_fit_parms*6+1)==0 || laufindex==1
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
    
    
    if laufindex<5
        sub=[E,model_parm.e1e2_substrate];
        save (fullfile(fit_parm.pathname_fit_ml,['substrate', ...
        num2str(laufindex),'.dat']), 'sub', '-ascii');
    end

    fitparms_xy=[E_e1_e2, [coeff_e2;zeros(numel(E)-numel(coeff_e2),1)], DRS(1:size(E_e1_e2,1))];
    %save (fullfile(fit_parm.pathname_fit_ml, 'fitparms', [num2str(laufindex),'.xy.fitparms']), ...
    %    'fitparms_xy', '-ascii');
    write_file(fullfile(fit_parm.pathname_fit_ml, 'fitparms', [num2str(laufindex),'.xy.fitparms']), ...
        fitparms_xy,...
        'Energy\t e1\t e2\t extcoeff\t DRSfit');
    
    positionVector1 = [0.1, 0.1, 0.3, 0.8];
    %subplot(1,7,[1,2])
    subplot('Position', positionVector1)
    min_e1=max(0,(max(E_e1_e2(:,2))-1.5));
    plot(E_e1_e2(:,1), [E_e1_e2(:,2)-min_e1,E_e1_e2(:,3)])
    legend(['e1xy-',num2str(min_e1)],'e2xy')    
    axis([min(E_e1_e2(:,1)) max(E_e1_e2(:,1)) 0 1.5])
    xlabel('Energy [eV]','FontSize',11)
    ylabel(['\epsilon _{2},\epsilon _{1}-', num2str(min_e1)],'FontSize',11)

    coeff=coeff_e2;
    coeff_index(:,1)=1:size(coeff,1);
    positionVector2 = [0.45, 0.1, 0.1, 0.8];
    %subplot(1,7,[3])
    subplot('Position', positionVector2)
    plot(coeff_index, coeff);
    axis([1 size(coeff,1) 0 1.5])
    xlabel('# coefficient','FontSize',11)
    
    positionVector3 = [0.65, 0.1, 0.3, 0.8];
    %subplot(1,7,[4,5])
    subplot('Position', positionVector3)
    drs_exp(:,1)=fit_parm.drs;
    plot(E_e1_e2(:,1), [DRS,drs_exp(:,1)])
    legend('DRS_{mod}','DRS_{exp}')
    xlabel('Energy [eV]','FontSize',11)
    ylabel('DRS','FontSize',11)
    title(['thickness ', num2str(model_parm.thickness), 'nm'])
    drawnow
    
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
    if sum(fit_parm.initial_parms.use(1:2)=='no')~=2
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
    
    E_e2x=[[E;max(E)+fit_parm.dEnergy], e2x, e10x];
    save (fullfile(fit_parm.pathname_fit_ml,['E_e2x', ...
        '.dat']), 'E_e2x', '-ascii');
    
    fitparms_xy=[E_e1_e2, [coeff_e2;zeros(numel(E)-numel(coeff_e2),1)]];
    save (fullfile(fit_parm.pathname_fit_ml, 'fitparms', ['finished','.xy.fitparms']), 'fitparms_xy', '-ascii');
end
end