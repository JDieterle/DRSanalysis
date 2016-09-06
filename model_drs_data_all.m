% transfer matrix algorithm for uniaxial anisotropic stratified samples
% simplified version of equations reported in:
% J. Hao and L. Zhou, PRB 77, 094201 (2008)
function data=model_drs_data_all(Energy_,ex_,d_, R_ss_0)
% aniso: vector with 1 for anisotropic and 1 for isotropic layers
if nargin==0
%Energy_(:,1)=1.4:0.01:6;
Energy_(:,1)=1:0.005:12;
R_ss_0=[];

% import material files:
% [~,si]=hdrload('Mat\si_jaw.mat');
% e1_e2_si(:,1)=interp1(si(:,1),si(:,2),Energy_);
% e1_e2_si(:,2)=interp1(si(:,1),si(:,3),Energy_);
% epsilon_si=e1_e2_si(:,1)-1i*e1_e2_si(:,2);
% 
% [~,sio2]=hdrload('Mat\sio2_jaw.mat');
% e1_e2_sio2(:,1)=interp1(sio2(:,1),sio2(:,2),Energy_);
% e1_e2_sio2(:,2)=interp1(sio2(:,1),sio2(:,3),Energy_);
% epsilon_sio2=e1_e2_sio2(:,1)-1i*e1_e2_sio2(:,2);
    h=6.6260755*10^(-34); %[Js]
    c=299792458; %[m/s]
    Joule_conversion=1.6*10^(-19);
    lambda_mu(:,1)=h*c./(Energy_(:,1)*Joule_conversion)*10^(6); %[µm]
    A_c=1.5147;
    B_c=0.0045386;
    C_c=-0.000022659;
    n_sub(:,1)=A_c+B_c./lambda_mu(:,1).^2+C_c./lambda_mu(:,1).^4; %real part
    ex_sub=n_sub.^2-1i*0;

epsilon_film_xy=Energy_.*0+1+1j*0;

%ex_=[epsilon_film_xy,epsilon_sio2,epsilon_si];
ex_=[epsilon_film_xy,ex_sub];
ex_=ex_sub;
%...%
num_layers=2;
num_layers=1;
  d_(:,1)=Energy_.*0+1000*10^-9;
% d_(:,1)=Energy_.*0+20*10^-9;
% d_(:,2)=Energy_.*0+130*10^-9;
% d_(:,3)=Energy_.*0+1000*10^-9;
% d_(:,2)=Energy_.*0+1000*10^-9;
assignin('base', 'ex_', ex_);
else
num_layers=size(ex_,2);
end

hbar=6.582119514*10^-16; % eV*s
c=299792458; % m/s
angle_of_incidence_=0;
Energy=zeros(numel(Energy_)*numel(angle_of_incidence_),1);
d=zeros(numel(Energy_)*numel(angle_of_incidence_),num_layers);
ex=complex(zeros(numel(Energy_)*numel(angle_of_incidence_),num_layers));

aoi=1;
Energy((aoi-1)*numel(Energy_)+1:aoi*numel(Energy_))=Energy_;
ex((aoi-1)*numel(Energy_)+1:aoi*numel(Energy_),:)=ex_;
d((aoi-1)*numel(Energy_)+1:aoi*numel(Energy_),:)=d_;


% data=[Energy,real(ex(:,1)),-imag(ex(:,1))];
% write_E_e1_e2_in_VASE_readable_file(data(size(data,1)/aoi,:), 'R:\xy.mat')
% data=[Energy,real(ez(:,1)),-imag(ez(:,1))];
% write_E_e1_e2_in_VASE_readable_file(data(size(data,1)/aoi,:), 'R:\z.mat')

%lE=numel(Energy);
omega_sq=(Energy./hbar).^2;
k_length=Energy/hbar/c;
kx=k_length.*0;
kz0_1(:,1)=k_length.*1;% kz in air, later replaced by kz(:,n_layer)
%kz0_2=-kz0_1;
%kz0_4=-kz0_1;

%%kl0_1=sqrt(kz0_1(:,1).^2+kx.^2);

for n_layer=1:num_layers
% calculate k_z
% for s-polarized light
k_z_1(:,1) = (sqrt(omega_sq.*ex(:,n_layer) - c^2*kx.^2)/c);
% for p-polarized light
% % if aniso(n_layer)==1
% % k_z_3(:,1) = sqrt(-(ex(:,n_layer).*( kx.^2-omega_sq/c^2.* ez(:,n_layer)))./ez(:,n_layer)); 
% % else
% % end
% length of k_vector
%%kl_1(:,1)=sqrt(k_z_1(:,1).^2+kx.^2);

Ms_1_1=(1+kz0_1./k_z_1)*0.5;
Ms_1_2=(1-kz0_1./k_z_1)*0.5;
Ms_2_1=(Ms_1_2);
Ms_2_2=(Ms_1_1);

if n_layer<num_layers
% propagation matrix:
P_1_1=exp(-k_z_1(:).*d(:,n_layer)*1i);
P_2_2=1./P_1_1;

% multiply P with M: P_n*M_n_(n-1)
pms_1_1=P_1_1.*Ms_1_1;
pms_1_2=P_1_1.*Ms_1_2;
pms_2_1=P_2_2.*Ms_2_1;
pms_2_2=P_2_2.*Ms_2_2;

else
pms_1_1=Ms_1_1;
pms_1_2=Ms_1_2;
pms_2_1=Ms_2_1;
pms_2_2=Ms_2_2;

end

% q: total transfer matrix
% start product with M_1_0 were layer #0 is air
if n_layer==1
qs_1_1=pms_1_1(:,1);
qs_1_2=pms_1_2(:,1);
qs_2_1=pms_2_1(:,1);
qs_2_2=pms_2_2(:,1);

else
qs_1_1_=qs_1_1;
qs_1_2_=qs_1_2;
qs_2_1_=qs_2_1;
qs_2_2_=qs_2_2;
    
qs_1_1=pms_1_1(:,1).*qs_1_1_+pms_1_2(:,1).*qs_2_1_;
qs_1_2=pms_1_1(:,1).*qs_1_2_+pms_1_2(:,1).*qs_2_2_;
qs_2_1=pms_2_1(:,1).*qs_1_1_+pms_2_2(:,1).*qs_2_1_;
qs_2_2=pms_2_1(:,1).*qs_1_2_+pms_2_2(:,1).*qs_2_2_;

end

% put k-values in variables for next M matrix
kz0_1=k_z_1;
end

r_ss=-qs_2_1./qs_2_2;
R_ss=r_ss.*conj(r_ss);

%t_ss=qs_1_1-qs_1_2.*qs_2_1./qs_2_2;
%T_ss=t_ss.*conj(t_ss);

%r_pp=-qp_2_1./qp_2_2;
%R_pp=r_pp.*conj(r_pp);

%t_pp=qp_1_1-qp_1_2.*qp_2_1./qp_2_2;
%T_pp=t_pp.*conj(t_pp);

% rho=r_pp./r_ss;
% absrho=abs(rho);
% psi=real(atan(absrho)*180/pi);
% del=real(-1i*log(rho./absrho)*180/pi);
% 
% data=[psi,del];
if isempty(R_ss_0)
    data=R_ss;
    assignin('base', 'R_ss_0', R_ss);
    else    
    data=(R_ss-R_ss_0)./R_ss_0;
end

if nargin==0
% plot(Energy, [R_pp, R_pp2-0.001, R_pp3+0.001])
try
[~,wvase]=hdrload('R:\r_p.dat');
wvase_e=reshape(wvase(:,1),[],aoi);
wvase_=reshape(wvase(:,4),[],aoi);
n=max(wvase(:,4));
catch
    disp('file R:\r_p.dat not found')
end
data_=reshape(r_ss,[],aoi);
m=max(data_);
figure
for index_=1:aoi
subplot(aoi,1,index_);
plot(Energy_, data_(:,index_));
hold on
subplot(aoi,1,index_);
try
m=max(data_(:,index_));
n=max(wvase_(:,index_));
plot(wvase_e(:,1), wvase_(:,index_),'--m');
catch
    disp('nothing else to plot')
end
end
end
%
% E_t_pp=[Energy, T_pp];
% E_t_ss=[Energy, T_ss];
% E_r_pp=[Energy, R_pp];
% E_r_ss=[Energy, R_ss];
% 
% save('R:\t_pp_m.dat', 'E_t_pp', '-ascii');
% save('R:\t_ss_m.dat', 'E_t_ss', '-ascii');
% save('R:\r_pp_m.dat', 'E_r_pp', '-ascii');
% save('R:\r_ss_m.dat', 'E_t_pp', '-ascii');