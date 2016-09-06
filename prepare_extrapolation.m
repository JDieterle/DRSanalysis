function extrapolation_data=prepare_extrapolation(Energy, pathname_fit, fit_parm)
%(E_lower_limit, E_upper_limit, dEnergy, pathname_fit)
E_lower_limit=min(Energy);
E_upper_limit=max(Energy);
dEnergy=abs(Energy(1)-Energy(2));
Energy2(:,1)=max(Energy):dEnergy:10*E_upper_limit-9*E_lower_limit;


quotient0=abs(kkt_(E_upper_limit+dEnergy,E_upper_limit)/kkt_(E_upper_limit+3*dEnergy,E_upper_limit))
%quotient0=abs(kkt_(E_upper_limit+dEnergy,E_upper_limit)-kkt_(E_upper_limit+2*dEnergy,E_upper_limit))


b=2;
element_list_E=[];
element_list_min=[];
element_list_max=[];

for a=1:1:numel(Energy2)
    quotient=abs(kkt_(Energy2(b),E_upper_limit)/kkt_(Energy2(a), E_upper_limit));
if quotient<quotient0
element_list_E(numel(element_list_E)+1,1)=a;

else
    b=a;
   element_list_min(numel(element_list_min)+1,:)=min(element_list_E);
   element_list_max(numel(element_list_max)+1,:)=max(element_list_E);
    element_list_E=[];
    element_list_E=a;
end
end

lmax= element_list_max+ numel(Energy);
lmin= element_list_min+ numel(Energy);

% 
ListE=[Energy2(element_list_min), Energy2(element_list_max), (Energy2(element_list_min)+Energy2(element_list_max))./2];


% add 1 column for constant offset (very high energies)
lmax(numel(lmax)+1)=lmax(numel(lmax));
lmin(numel(lmin)+1)=lmin(numel(lmin));

%list of all needed summands

%energy verctor for higher energies
num_en_experiment=numel(Energy);
maxE_exp=max(Energy);
maxE_extra=max(Energy2);
minE=E_lower_limit;
 
E_k_extra(:,1)=E_lower_limit:dEnergy:maxE_extra;
num_en_extrapolated=numel(E_k_extra(:,1));


% % extrapolation data:::
for a=1:numel(lmax)
E_k_extra(:,a+1)=zeros(numel(E_k_extra,1),1);
E_k_extra(lmin(a):lmax(a),a+1)=1;
end

E_k_extra(:,numel(lmax)+1)=0;


fit_parm_local=struct('extrapol', struct('E_e2', E_k_extra));

% % number of summands to extrapolate k:xxxxxxxx
%fit_parm.extrapol.
num_summands_k=size(E_k_extra,2)-1;
%fit_parm.extrapol.num_coeff_k=fit_parm.extrapol.num_summands_k;

%num_summands_k=fit_parm.extrapol.num_summands_k;
%num_coeff_k=fit_parm.extrapol.num_coeff_k;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%### aus k n ausrechnen...
for a=1:num_summands_k
n_extra(:,a)=KKTx(E_k_extra(:,1),...
    [zeros(num_en_experiment,1);...
    E_k_extra(num_en_experiment+1:num_en_extrapolated,a+1)],0);
end

%for very high energies where n ist a constant:
n_extra(:,num_summands_k)=1;


if max(n_extra(:,num_summands_k-1))<10^-7
    mmm=n_extra
    n_extra(:,num_summands_k-1)=[];
    E_k_extra(:,num_summands_k)=[];
    num_summands_k=size(E_k_extra,2)-1;
    remove=1
    lmax(numel(lmax)-1)=[];
    lmin(numel(lmin)-1)=[];
end



% n in fit_parm speichern:
extrapol.n=n_extra;
n_E=[E_k_extra(:,1),n_extra];
E_e2=[E_k_extra(:,1),...
    [zeros(num_en_experiment,num_summands_k);...
    E_k_extra(num_en_experiment+1:num_en_extrapolated,2:num_summands_k+1)]];

save (strcat(pathname_fit,'offset_e1.dat'), 'n_E', '-ascii');
save (strcat(pathname_fit,'offset_e2.dat'), 'E_e2', '-ascii');


%### extrapolation vorbereiten fertig

%prepare weigth data for smoothing:
Energy_range(:,1)=max(abs(lmax-lmin),1);
weight(:,1)=ones(numel(Energy),1);
weight(numel(Energy)+1:numel(Energy)+numel(lmin),1)=1./Energy_range;
weight(round(numel(Energy)/2):numel(Energy),1)=1;
weight(numel(Energy)-numel(lmin):numel(Energy),1)=1;
weight(numel(Energy),1)=1;
weight(numel(Energy)+numel(lmin),1)=0;
%weight(numel(Energy)-2*numel(lmin):numel(Energy)+round(numel(lmin)/2))=weight(numel(Energy)-2*numel(lmin):numel(Energy)+round(numel(lmin)/2))*5;
%weight(numel(Energy)-round(numel(lmin)/2):numel(Energy)+round(numel(lmin)/4))=weight(numel(Energy)-round(numel(lmin)/2):numel(Energy)+round(numel(lmin)/4))*5;
weight(1:round(numel(Energy)/7))=1;
%weight(numel(Energy)-24:numel(Energy)-21)=500;
%...
weight(1:200)=1;


E_weight(:,1)=E_lower_limit:dEnergy:2*E_upper_limit-E_lower_limit;
E_weight(2:numel(weight)+1,2)=weight;
save (strcat(pathname_fit,'E_weight.dat'), 'E_weight', '-ascii');

%########
extrapolation_data=struct('e1_extra',extrapol.n, 'num_summands_e2',num_summands_k, 'E_e2_extra',E_e2, 'weight', weight, 'ListE', ListE, 'E_upper_limit',E_upper_limit);

end


function n=kkt_(Eabs,E_n)

n=Eabs/(Eabs^2-E_n^2);

end