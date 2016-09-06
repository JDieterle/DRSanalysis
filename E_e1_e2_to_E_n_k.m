function E_n_k=E_e1_e2_to_E_n_k(E_e1_e2)
E=E_e1_e2(:,1);
e1=E_e1_e2(:,2);
e2=E_e1_e2(:,3);


n=sqrt((sqrt(e1.^2+e2.^2)+e1)./2);
k=sqrt((sqrt(e1.^2+e2.^2)-e1)./2);

E_n_k=[E,n,k];
end