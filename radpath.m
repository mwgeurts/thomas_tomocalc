function  rpath=radpath(ctimages,xf,xp,zf,zp,MVCT)
global dose_dimensiony
%xf,zf single focus value (outside array)
%xp,zp multiple calc points in array
dx=xf-xp;
dz=zf-zp;
steps=max((max(abs(dx),abs(dz))));
r=sqrt(dx.*dx + dz.*dz);
dr=r/steps;
deltax=dx/steps;
deltaz=dz/steps;
x=xp;
z=zp;
rpath=0;
ctdim=max(size(ctimages));
ctsize=ctdim*ones(1,2);
%hutemp=zeros(3,15,225);
for i=1:steps
    iz=floor(z);
    iz(iz<1)=1; iz(iz>size(ctimages,3))=1;
    ix=floor(x);   
    ix(ix<1)=1; ix(ix>size(ctimages,2))=1; %use 1 as flag both sides
    if max(iz)==1, break, end
    if max(ix)==1, break, end  
    rhoe=ctimages(:,sub2ind([size(ctimages,2) size(ctimages,3)],ix,iz)); % electron_density() function was removed as the CT was already converted to density
    rpath=rpath+rhoe;  %will multiply by dr later

    %hutemp(1,:,i)=hu(:,113);
    %hutemp(2,:,i)=rhoe(:,113);
    %hutemp(3,:,i)=rpath(:,113);
    x=x+deltax;
    z=z+deltaz;
end
a=ones(dose_dimensiony,1);
dr2=a*dr;
rpath=rpath.*dr2;
return;
