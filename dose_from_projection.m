function dose=dose_from_projection(Yvaluemm,Xtheta,depthmm,dfromfocmm,gantry_period,ddepth)
global SP TPR OARXOPEN OARXLEAVES OARY segments;
global reference_doserate
%changed Sept 2013 to make all inputs in mm
%internal to this function only, will use cm
Yvalue=Yvaluemm/10;
depth=depthmm/10;
dfromfoc=dfromfocmm/10;
%
dose=zeros(1,size(Yvalue,2));

depthx=depth;  %depth for profiles - leave original for TPR
depthx(depthx>OARY.depths(OARY.ndepths))=OARY.depths(OARY.ndepths);
% this code assumes max depth same for all three sets of profiles
    depthx(depthx<OARY.depths(1))=OARY.depths(1);
% this code assumes max/min depth same for all three sets of profiles

%next lines added SJT 29/9/14 to trap out of range values
        maxopenindex=max(OARXOPEN.indices);
        maxleafindex=max(OARXLEAVES.indices);
        maxyindex=max(OARY.indices);

for seg=1:segments.nseg
    leaf1=segments.startval(seg);
    lastleaf=segments.endval(seg);
    weight=segments.value(seg);
    xnleaves=1+lastleaf-leaf1;
    nleaves=int16(xnleaves);
    width=(1+lastleaf-leaf1)*0.625; %0.625cm per leaf at 85cm
    midpos=(0.5*(leaf1+lastleaf) - 32.5)*0.625;  %mid between leaf 32 and 33
    eqsq=2*SP.length*width/(SP.length+width); 
    tprvalue=interp2(TPR.sizes,TPR.depths,TPR.tpr,eqsq,depth);
    tprvalue(depth<=0)=0; %set any outside phantom to zero
    spvalue=interp1(SP.sizes,SP.values,width);
%calculate index for 40-wide x-profile
        theta=abs(Xtheta);
        oarindex=theta/0.005; %5 milliradian steps
   

        oarindex(oarindex>maxopenindex)=maxopenindex;%added SJT 29/9/14 to trap out of range values
        
        oar40value=interp2(OARXOPEN.depths,OARXOPEN.indices,OARXOPEN.oar,depthx,oarindex);
    %calculate profile for field set by leaves
        theta=abs(Xtheta-atan(midpos/85));
        if nleaves>2
            thetaoffset=atan((xnleaves-4)*.5*.625/85.0);
            theta=theta-thetaoffset;
            theta(theta<0)=0;
            nleaves=3; %index to third leaf size
        end
        oarindex=theta/0.001; %milliradians
        oarindex(oarindex>maxleafindex)=maxleafindex;%added SJT 29/9/14 to trap out of range values
        oarxvalue=interp2(OARXLEAVES.depths,OARXLEAVES.indices,OARXLEAVES.oar(:,:,nleaves),depthx,oarindex,'linear',0);
        oarxv=oarxvalue.*oar40value;
    %calculate index for y-profile
        theta=atan(abs(Yvalue./(85-ddepth')));
  %      theta=atan(abs(Yvalue/85));
        if width<OARY.widths(1)
            width=OARY.widths(1); %this is OK, because width is a scalar
        end
        wvect=repmat(width,1,size(oarindex,1)); %gets it vectorised the same as oarindex
        oarindex=theta/0.001;  % milliradians
        oarindex(oarindex>maxyindex)=maxyindex;%added SJT 29/9/14 to trap out of range values
        oaryvalue=interp3(OARY.depths,OARY.indices,OARY.widths,OARY.oar,depthx',oarindex,wvect,'linear',0);
        
        dose=dose+ weight*spvalue*(tprvalue'.*oarxv'.*oaryvalue);

   
end
dfoc2=dfromfoc.*cos(Xtheta);
invsq=(85*85)./(dfoc2.*dfoc2);%inverse square from 85cm
%0.0125986 is (86.5/85)^2 x 1/1.37 x 1/60   (1.37 is TPR)
dose=dose.*invsq'*reference_doserate*0.0125986*gantry_period/51.0; %normalise to 5x40 field
return