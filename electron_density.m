function rhoe=electron_density(hu,ct_MV)
%convert HU to ED using Thomas 1999 for kV, single line if MV
fhu=double(hu); 
if ct_MV
    rhoe=(fhu+945.1)/950.1;
    rhoe(rhoe<0)=0; % trap negative densities
else
    rhoe=1+fhu/1000;
    rhoe2=1+fhu/1950;
    rhoe(rhoe<0)=0; % trap negative densities
    rhoe(rhoe>1.1)=rhoe2(rhoe>1.1); %bone
end
return;