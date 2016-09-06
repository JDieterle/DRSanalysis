% Energy: Vector equally spaced

function extrapolation_data=prepare_extrapolation_old(Energy, pathname_fit, fit_parm)
E_lower_limit=min(Energy);
E_upper_limit=max(Energy);

gg=Energy;

dEnergy=abs(Energy(2)-Energy(1));


try
    extrapolation_data=evalin('base', 'extrapolation_data');
    if abs(E_lower_limit-extrapolation_data.E_lower_limit)<0.000001 && ...
       abs(E_upper_limit-extrapolation_data.E_upper_limit)<0.000001 && ...
       abs(dEnergy      -extrapolation_data.dEnergy)<0.000000001
       return
    end
catch
    disp('prepare extrapolation...')
end

Energy2(:,1)=E_upper_limit:dEnergy:10*E_upper_limit-9*E_lower_limit;

% calculate deviation in kkt:
quotient0=abs(kkt_(E_upper_limit+dEnergy,E_upper_limit)/kkt_(E_upper_limit+2*dEnergy,E_upper_limit))
%quotient0=abs(kkt_(E_upper_limit+dEnergy,E_upper_limit)-kkt_(E_upper_limit+2*dEnergy,E_upper_limit))

b=2;
element_list_E=[];
element_list_min=[];
element_list_max=[];

for a=1:1:numel(Energy2)
    quotient=abs(kkt_(Energy2(b),E_upper_limit)/kkt_(Energy2(a), E_upper_limit));
    if quotient<quotient0
        element_list_E=[element_list_E;a];
    else
        b=a;
        element_list_min=[element_list_min; min(element_list_E)];
        element_list_max=[element_list_max; max(element_list_E)];
        element_list_E=[];
        element_list_E=a;
    end
end

lmax= min(element_list_max+ numel(Energy),numel(Energy2));
lmin= element_list_min+ numel(Energy);
assignin('base', 'prepare_extrapolationlmax', lmax)
assignin('base', 'prepare_extrapolationlmin', lmin)
assignin('base', 'prepare_extrapolationEnergy2', Energy2)
ListE=[Energy2(lmin), Energy2(lmax), (Energy2(lmin)+Energy2(lmax))./2];
% add 1 column for constant offset (very high energies)
lmax=[lmax; numel(Energy2)];
lmin=[lmin; max(lmax)];

% list of all needed summands
% energy vector for higher energies
num_en_experiment=numel(Energy);
maxE_extra=max(Energy2);
 
E_e2_extra(:,1)=E_lower_limit:dEnergy:maxE_extra;
num_en_extrapolated=numel(E_e2_extra(:,1));

% extrapolation data:::
for a=1:numel(lmax)
    E_e2_extra(:,a+1)=zeros(numel(E_e2_extra,1),1);
    E_e2_extra(lmin(a):lmax(a),a+1)=1;
    E_e2_extra(lmin(a),a+1)=0.75;
    %if a>1
    E_e2_extra(lmin(a)-1,a+1)=0.25;
    %end
    E_e2_extra(lmax(a),a+1)=0.75;
    E_e2_extra(lmax(a)+1,a+1)=0.25;

end
%fit_parm_local=struct('extrapol', struct('E_e2', E_e2_extra));

% number of summands to extrapolate e2:
num_summands_e2=size(E_e2_extra,2)-1;

%### calculate effect of e2 on e1 by KKT:
kk=[];
for a=1:num_summands_e2
    k=E_e2_extra(:,a+1);
    %k(:,1)= [0.75*k(1)+0.25*k(2); ...
    %        k(1:numel(k)-2)*0.25+k(2:numel(k)-1)*0.5+k(3:numel(k))*0.25;...
    %        0.25*k(numel(k)-1)+0.75*k(numel(k))];
     e1_extra(:,a)=KKTx(E_e2_extra(:,1),k,0);
end

% for a=1:num_summands_e2
%         e1_extra(:,a)=KKTx(E_e2_extra(:,1),...
%     [zeros(num_en_experiment,1);...
%     E_e2_extra(num_en_experiment+1:num_en_extrapolated,a+1)],0);
% end



%for very high energies which result in a constant offset in e1:
% set e2=0 in column for constant offset (very high energies):
E_e2_extra(:,numel(lmax)+1)=0;
% set e1=1:
e1_extra(:,num_summands_e2)=1;

for a=1:num_summands_e2-1
    if max(e1_extra(1:numel(Energy),num_summands_e2-1))<10^-4
        e1_extra(:,num_summands_e2-1)=[];
        E_e2_extra(:,num_summands_e2)=[];
        num_summands_e2=size(E_e2_extra,2)-1;
        disp('remove unnecessary summands')
        lmax(numel(lmax)-1)=[];
        lmin(numel(lmin)-1)=[];
    end
end

% e1 in fit_parm speichern:
extrapol.e1=e1_extra;
E_e1=[E_e2_extra(:,1),e1_extra];
E_e2=[E_e2_extra(:,1),...
    [...
    E_e2_extra(1:num_en_extrapolated,2:num_summands_e2+1)]];

%    [zeros(num_en_experiment,num_summands_e2);...
%    E_e2_extra(num_en_experiment+1:num_en_extrapolated,2:num_summands_e2+1)]];


write_file(fullfile(pathname_fit,'offset_e1.dat'), E_e1, 'Energy\t e1');
write_file(fullfile(pathname_fit,'offset_e2.dat'), E_e2, 'Energy\t e2');
%E_e2_all=[E_e1(:,1),kk];
%write_file(fullfile(pathname_fit,'E_e2.dat'),      E_e2_all, 'Energy\t e2');


%### extrapolation vorbereiten fertig

% %prepare weigth data for smoothening:
% Energy_range(:,1)=max(abs(lmax-lmin),1);
% weight(:,1)=ones(numel(Energy),1);
% weight(numel(Energy)+1:numel(Energy)+numel(lmin),1)=1./Energy_range;
% weight(numel(Energy)-5*numel(lmin):numel(Energy),1)=15;
% weight(numel(Energy),1)=2;
% weight(numel(Energy)+numel(lmin),1)=0;
% %weight(numel(Energy)-2*numel(lmin):numel(Energy)+round(numel(lmin)/2))=weight(numel(Energy)-2*numel(lmin):numel(Energy)+round(numel(lmin)/2))*5;
% %weight(numel(Energy)-round(numel(lmin)/2):numel(Energy)+round(numel(lmin)/4))=weight(numel(Energy)-round(numel(lmin)/2):numel(Energy)+round(numel(lmin)/4))*5;
% weight(1:round(numel(Energy)/5))=1;
% weight(1:fit_parm.upper_limit_small_absorption_eV_index(2))=1;
% weight(1:fit_parm.upper_limit_no_absorption_eV_index(2))=1;
% 
% E_weight(:,1)=E_lower_limit:dEnergy:2*E_upper_limit-E_lower_limit;
% E_weight(2:numel(weight)+1,2)=weight;
% save (strcat(pathname_fit,'E_weight.dat'), 'E_weight', '-ascii');

%########
extrapolation_data=struct('e1_extra',extrapol.e1,...
                          'num_summands_e2',num_summands_e2,...
                          'E_e2_extra',E_e2, ...
                          'weight', [], 'ListE', ListE,...
                          'E_lower_limit',E_lower_limit,...
                          'E_upper_limit',E_upper_limit, ...
                          'lmax', lmax, 'lmin',lmin,...
                          'dEnergy', dEnergy);
                      
assignin('base', 'extrapolation_data', extrapolation_data);
end

function e1=kkt_(Eabs,Ee1)
    e1=Eabs/(Eabs^2-Ee1^2);
end