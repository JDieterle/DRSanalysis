function ev=evaluate_e2(e2, coeff_e2,weight,fit_parm, DRS_fit)
pathname_fit_ml=fit_parm.pathname_fit_ml;
drs_exp(:,1)=fit_parm.drs;

e2_c=[e2;coeff_e2];
num_e2_c=numel(e2_c);
num_e2=numel(e2);
weight_coeff(:,1)=weight(numel(e2)+3:numel(weight));

% difference e2(n) e2(n+1)
d_e2_e2(:,1)=e2(2:num_e2,1)-e2(1:num_e2-1,1).^2;

dd_ev=0.02*d_e2_e2'*(d_e2_e2.*weight(2:num_e2));
%curvature:
nn(:,1)=(e2_c(3:num_e2_c,1)+e2_c(1:num_e2_c-2,1)-2*e2_c(2:num_e2_c-1,1)).^2;
nnw=nn.*weight(3:numel(weight));
ev_=fit_parm.weights.st*nn;

%correlate st to drs_exp:
drs_st=fit_parm.weights.st_del*abs(drs_exp-DRS_fit)'*abs(nnw(1:numel(e2)));

% couple first coeff_e2 to 2.last e2:
%difference coeff 1 and e2 end:
diff_c1_e2=abs(coeff_e2(1)-e2(numel(e2)-1)).^2+0.25*abs(coeff_e2(1)-e2(numel(e2))).^2+0.15*abs(coeff_e2(1)-e2(numel(e2)-2)).^2;
%difference coeff 2 and e2 end:
diff_c2_e2=0.3*abs(coeff_e2(2)-e2(numel(e2)-1)).^2+0.5*abs(coeff_e2(2)-e2(numel(e2))).^2+0.2*abs(coeff_e2(2)-e2(numel(e2)-2)).^2;
%difference coeff 2 and coeff 1:
diff_coeff_(:,1)= abs(coeff_e2(2:numel(coeff_e2)-1)-coeff_e2(1:numel(coeff_e2)-2));

diff_coeff_w(:,1)=diff_coeff_.*weight_coeff.^2;
diff_coeff=diff_coeff_'*diff_coeff_w;%+abs(coeff_e2(2)-coeff_e2(1));


difference=fit_parm.weights.st_diff_c*(20*diff_c1_e2+0.7*diff_c2_e2+diff_c2_e2^4+diff_coeff);
difference=difference*0;
%
num_points_without_absorption=numel(fit_parm.e2_eqal_zero(:,1));
diff_exp_fit_trans=(DRS_fit(1:num_points_without_absorption)-drs_exp(1:num_points_without_absorption));
e2_transparent_range=(diff_exp_fit_trans)*fit_parm.weights.diff_exp_mod_transparent_range;
%
lower_end=mean(diff_exp_fit_trans(1:max(1,round(num_points_without_absorption/3))));
upper_end=mean(diff_exp_fit_trans(max(1,round(num_points_without_absorption/2*1)):...
                                  max(1,round(num_points_without_absorption/6*5))));
e2_curv_transparent_range=(lower_end-upper_end)*fit_parm.weights.e2_curv_transparent_range*sum(abs(e2_transparent_range));

e2small=(e2(fit_parm.upper_limit_no_absorption_eV_index(2)+1:fit_parm.upper_limit_small_absorption_eV_index(2)).*fit_parm.weight_e2_small);

ev=struct('st' ,ev_, 'diff',difference, 'diff_ev', dd_ev, 'drs_st', drs_st,...
    'e2_transparent_range', e2_transparent_range, 'e2_curv_transparent_range',...
    e2_curv_transparent_range, 'e2small', e2small) ;

global laufindex 
if mod(laufindex+1,1000)==0
    if laufindex+1==10000000 % fit ready, save all data
        E_nn(2:numel(nn)+1,2)=nn;
        save (strcat(pathname_fit_ml,'E_st','.dat'), 'E_nn', '-ascii');
         ready_ev=1
    end
end


end