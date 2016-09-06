%Kramers-Kronig Transformation

%Energy(:,1) Spaltenvektor mit Energie
%k(:,1) Spaltenvektor mit k
%
%=> n(:,1) Spaltenvektor mit n


function n = KKTx(Energy,k, n_offset)
n(1:numel(Energy),1)=zeros;
energies_odd_number(:,1)=Energy(1:2:numel(Energy),1);
energies_odd_number_quadrat(:,1)=energies_odd_number.^2;
%k(:,1)=[0.75*k(1)+0.25*k(2); k(1:numel(k)-2)*0.25+k(2:numel(k)-1)*0.5+k(3:numel(k))*0.25; 0.25*k(numel(k)-1)+0.75*k(numel(k))];
k_odd_number(:,1)=k(1:2:numel(k),1);
E_k_j_odd(1,:)=energies_odd_number.*k_odd_number;
energies_even_number(:,1)=Energy(2:2:numel(Energy),1);
energies_even_number_quadrat(:,1)=energies_even_number.^2;
k_even_number(:,1)=k(2:2:numel(k),1);
E_k_j_even(1,:)=energies_even_number.*k_even_number;
Delta_E=(Energy(2,1)-Energy(1,1));


for i=1:2:numel(Energy)
    %E_ij=E_j^2-E_i^2
    E_ij=1./(energies_even_number_quadrat -(Energy(i))^2);
    % finish evaluation
    n(i,1)= E_k_j_even*E_ij;
    %end    
end

for i=2:2:numel(Energy)
    %E_ij=E_j^2-E_i^2
    E_ij=1./(energies_odd_number_quadrat -(Energy(i))^2);
    % finish evaluation
    n(i,1)= E_k_j_odd*E_ij;
    %end
end

n=n_offset+4 / pi *Delta_E*n;
%n=arrayfun(@(k)KKTelement(num, Energy,k,n_offset),k)
end

