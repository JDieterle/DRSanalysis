function imported_data=import_E_e1_e2(filename,fit_parm,data,pathname_import)
if nargin==4
    [~,name,ext] = fileparts(filename);
    try
    copyfile(filename,fullfile(pathname_import,[name,ext]));
    catch
        disp(['import_E_e1_e2: unable to copy file: ', filename]);
    end
end


Energy(:,1)=data.Energy;

disp(strcat('importing_', filename, '...'));

imp_data=load_and_interpolate_file(filename, Energy);
E_e1_e2=imp_data.E_e1_e2;
E_min=imp_data.E_min;
index_E_min=energy_to_index(E_e1_e2(:,1), E_min);
fit_parm.index_E_min=index_E_min;

nf=numel(filename);
end_filename=filename(nf-7:nf);

e2extrapolation_coefficients=[];
if sum(end_filename=='fitparms')==8
    try
    [~,rawdata]=hdrload(filename);
    e2extrapolation_coefficients=rawdata(1:fit_parm.extrapol.num_coeff_e2,4)
    catch
        disp('import_E_e1_e2: error loading *.fitparms')
    end
end
nf=numel(filename);
end_filename=filename(nf-3:nf);
if sum(end_filename=='.dat')==4
    try
        disp(['   ', 'import extrapolation parameter directly form .dat file'])
        imp_data=load_and_interpolate_file(filename, Energy, [1,3,4]);
        E_e1_e2=imp_data.E_e1_e2;
        E_min=imp_data.E_min;
        index_E_min=energy_to_index(E_e1_e2(:,1), E_min);
        fit_parm.index_E_min=index_E_min;
        [~,rawdata]=hdrload(filename);
        e2extrapolation_coefficients=rawdata(1:fit_parm.extrapol.num_coeff_e2,5)
        disp('read *.dat file')
    catch
        disp('import_E_e1_e2: error loading *.dat')
    end
end
if isempty(e2extrapolation_coefficients)
        disp('import_E_e1_e2: fitting extrapolation coeff.')
    e2extrapolation_coefficients=tranform_e1_e2_to_e2_e2extrapolation_coefficients(E_e1_e2,fit_parm,data,filename,index_E_min,pathname_import);  
end 
imported_data=struct('E_e1_e2',E_e1_e2,'e2extrapolation_coefficients',e2extrapolation_coefficients);

if nargin==4
    [~,name,ext] = fileparts(filename);
    write_file(fullfile(pathname_import,['E_e1_e2_',name,ext, '_imported']),E_e1_e2, ...
    'Energy\t e1\t e2')
end
disp(strcat('...imported!'));
end


function e2extrapolation_coefficients=tranform_e1_e2_to_e2_e2extrapolation_coefficients(E_e1_e2,fit_parm,data,filename,index_E_min,pathname_import)
E_lower_limit=min(data.Energy);
E_upper_limit=max(data.Energy);
extrapolation_data=evalin('base', 'extrapolation_data');

% fit e2extrapolation_coefficients to difference:

e2extrapolation_coefficients(:,1)=ones(extrapolation_data.num_summands_e2,1)*1;

e2_end=E_e1_e2(size(E_e1_e2,1),3);

lb=e2extrapolation_coefficients.*0;
ub=e2extrapolation_coefficients+16;


%new
e2extrapolation_coefficients(:,1)=ones(extrapolation_data.num_summands_e2,1)*e2_end;
e2extrapolation_coefficients(numel(e2extrapolation_coefficients),1)=1;
e2extrapolation_coefficients(1:3,1)=e2_end;

lb(1)=e2_end*0.9;
ub(1)=max(0.2,e2_end*1.1);
lb(2)=e2_end*0.8;
ub(2)=max(0.2,e2_end*1.2);

options = optimset('Algorithm','trust-region-reflective','DiffMinChange', 0.0001,...
        'ScaleProblem','jacobian','FinDiffType','central',...
        'Display', 'iter','MaxFunEvals', 50000, 'TolFun', 1e-10);

    e2extrapolation_coefficients=...
    lsqcurvefit(@(e2extrapolation_coefficients, Energy)...
    calculate_e1...
    (Energy, E_e1_e2(index_E_min:numel(E_e1_e2(:,2)),:), fit_parm, data, e2extrapolation_coefficients, extrapolation_data, filename),...
    e2extrapolation_coefficients,...
    data.Energy(index_E_min:numel(data.Energy),1),...
    E_e1_e2(index_E_min:numel(E_e1_e2(:,2)),2),...
    lb,...
    ub,options);

fit_parm.index_E_min=0;

E_e1_e2_e1o_e2o(:,[1,4,5])=E_e1_e2;
E_e1_e2_e1o_e2o(:,3)=E_e1_e2(:,3);
E_e1_e2_e1o_e2o(:,2)=...
    calculate_e1...
    (data.Energy, E_e1_e2, fit_parm, data,...
    e2extrapolation_coefficients, extrapolation_data, filename);

% E_e1_e2_x(:,1)=E_e1_e2(:,1);
% E_e1_e2_x(:,3)=E_e1_e2(:,3)*0+0.1;
% E_e1_e2_e1o_e2o(:,2)=...
%     calculate_e1...
%     (data.Energy, E_e1_e2_x, fit_parm, data,...
%     e2extrapolation_coefficients.*0+0.1, extrapolation_data, filename);
length_E_long=size(extrapolation_data.E_e2_extra,1);
length_E_exp=size(E_e1_e2,1);
E_e2_long=[extrapolation_data.E_e2_extra(:,1),[E_e1_e2(:,3);zeros(length_E_long-length_E_exp,1)]+...
          sum(extrapolation_data.E_e2_extra(:,2:size(extrapolation_data.E_e2_extra,2))*e2extrapolation_coefficients,2)];
%.............................................E_e1_e2(:,3)
last_e2extrapolation_coefficient=e2extrapolation_coefficients(numel(e2extrapolation_coefficients));
E_e2_long(:,3)=KKTx(E_e2_long(:,1),E_e2_long(:,2),1+last_e2extrapolation_coefficient);
E_e2_=E_e2_long(1:length_E_exp,:);

%plot(E_e1_e2_e1o_e2o(:,1),[E_e1_e2_e1o_e2o(:,2:3),E_e2_(:,2)+0.01]);
plot(E_e2_long(:,1), E_e2_long(:,2:3));
hold on
plot(E_e1_e2_e1o_e2o(:,1),[E_e1_e2_e1o_e2o(:,2:3),E_e2_(:,2)+0.01]);
plot(E_e1_e2(:,1),E_e1_e2(:,2:3));


%save (strcat(filename,'_import'), 'E_e1_e2_e1o_e2o', '-ascii');
if nargin==6
    [~,name,ext] = fileparts(filename);
    write_file(fullfile(pathname_import,['E_e1_e2_e1o_e2o',name,ext, '_imported']),E_e1_e2_e1o_e2o, ...
    'Energy\t e1imp\t e2imp\t e1orig\t e2orig')
end

E_e2_extra(:,1)=extrapolation_data.E_e2_extra(:,1);
for i=1 : extrapolation_data.num_summands_e2
    e2_extra(:,i)=  [zeros(numel(data.Energy),1);...
                       extrapolation_data.E_e2_extra(numel(data.Energy)+1:...
                       numel(E_e2_extra(:,1)),1+i)].*...
                       e2extrapolation_coefficients(i);
end
E_e2_extra=[E_e2_extra,e2_extra];

%save (strcat(filename,'_import_E_e2_extra'), 'E_e2_extra', '-ascii');
if nargin==6
    [~,name,ext] = fileparts(filename);
    write_file(fullfile(pathname_import,['E_e2_extra',name,ext, '_imported']),E_e2_extra, ...
    'Energy\t e2extra')
end

%x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub)
end

function e1_kkt_e2=calculate_e1(Energy, E_e1_e2, fit_parm, data, e2extrapolation_coefficients, extrapolation_data, filename)
% kramers... transformation
e2=E_e1_e2(:,3);
e2x=e2;

e2_end=e2(numel(e2));%xxxx
e2(numel(e2))=e2_end*0.75;%xxxx

%e2x=[e2;0]; % append 0 for special stable KKT
%e2=e2alt;
%e2(numel(e2))=e2_end*0.75+e2extrapolation_coefficients(1)*0.25;%xxxx
%e2x(numel(e2x))=e2_end*0.25;%xxxx
%
%e2x(numel(e2x))=e2extrapolation_coefficients(1)*0.25;%xxxx
%e2x(numel(e2))=e2_end-e2extrapolation_coefficients(1)*0.25;


%e10=KKTx([Energy;fit_parm.dEnergy+max(Energy)],e2x,1); %n_offset=1
e10=KKTx(Energy,e2x,1); %n_offset=1

%e10(numel(e10))=[];

% add summand for higher energies xy:
e1_kkt_e2(:,1)=e10+ sum( extrapolation_data.e1_extra(fit_parm.index_E_min+1:fit_parm.index_E_min+numel(e10),:)* e2extrapolation_coefficients,2);

E_=[Energy, e10, e1_kkt_e2, E_e1_e2(:,3)];
%save (strcat(filename,'ccccc'), 'E_', '-ascii');
if nargin==4
    [~,name,ext] = fileparts(filename);
    write_file(fullfile(pathname_import,['E_e1',name,ext, '_imported']),E_, ...
    'Energy\t e1nooff\t e1wioff\t e1\t e2\t')
end
end

function imp_data=load_and_interpolate_file(filename, Energy, columns)
[~,E_e1_e2_raw_]=hdrload(filename);
if nargin>2
    columns=columns;
else
    columns=1:3;
end

E_e1_e2_raw=E_e1_e2_raw_(:,columns);

E_e1_e2_raw=make_E_ascending_vector(E_e1_e2_raw);

E_min=min(E_e1_e2_raw(:,1));

if min(E_e1_e2_raw(:,1))> min(Energy)
   E_e1_e2_raw=[[min(Energy),E_e1_e2_raw(1,2),E_e1_e2_raw(1,3)];E_e1_e2_raw];
end

if max(E_e1_e2_raw(:,1))< max(Energy)
   E_e1_e2_raw=[E_e1_e2_raw;[max(Energy),E_e1_e2_raw(size(E_e1_e2_raw,1),2),...
               E_e1_e2_raw(size(E_e1_e2_raw,1),3)]];
end

E_e1_e2(:,1)=Energy;
E_e1_e2(:,2:3)=...
        interp1(E_e1_e2_raw(:,1),E_e1_e2_raw(:,2:3),Energy,'pchip');
    
    
    imp_data=struct('E_e1_e2', E_e1_e2, 'E_min',E_min);
    assignin('base', 'import_E_e1_e2', E_e1_e2)
end

% wandelt eine Energieposition in den Index der nächsthöheren Energie im Energievektor um:
% Parameter:    Energy: vektor mit Energie
%               Energyposition: Position des Punktes im Spektrum (eV)
% Ergebnis:     index: Index des Energiewertes im Vektor Energy, der knapp
%               unterhalb von Energyposition liegt
function index=energy_to_index(Energy, Energyposition)
a=1;
while Energy(a)<Energyposition % falls Energie zu groß/klein funktionierts nicht
    a=a+1;
end
index=a;
end

function E_e1_e2_raw=make_E_ascending_vector(E_e1_e2_raw);
if E_e1_e2_raw(1,1)>E_e1_e2_raw(2,1)
    E_e1_e2_rawx=[];
    for a=1:size(E_e1_e2_raw,1)
    E_e1_e2_rawx=[E_e1_e2_raw(a,:);E_e1_e2_rawx];
    end
    E_e1_e2_raw=E_e1_e2_rawx;
else
    E_e1_e2_raw=E_e1_e2_raw;
end
end