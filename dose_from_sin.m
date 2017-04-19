function dose_from_sin(dvhonly,datasource,handles)
global ctimages;
global projection
global dosecube ct_ylist ct_yindex dose_dimensionxz dose_dimensiony ct_xindex ct_zindex
global matcen_x matcen_y matcen_z dose_gridxz dose_gridy
global Xvalue Zvalue YvalueRef start_y
global dcmpath dcmfile dose_tol dist_tol nfractions perctol
global pitch gantry_period field_width
global gantry_start sino_path sino_name ct_path ct_filename
global dosecubeinfo y_offset_gamma dosecubeindex plotlimit
global x_offset_gamma z_offset_gamma
global dicomdosecube n_of_proj TomoPositionPatient
global ct_rounding
global drt_loaded sinogram
global ctpix isoc_pos DefIndex mindepth
global amp_4d phase_4d min_index4d max_index4d
global rpm_data rpm_mark ctpatientposition
global comparetwo
global ctImPosPat num_of_subprojections original_gantry_start %GST change
global ct_sorted_names ct_sorted_couch ct_sorted_index
global TomoRoll
global sec_ct_path sec_ct_filename sec_ct_sorted_names sec_ct_sorted_couch sec_ct_sorted_index sec_ctpix
global sec_ct_patient_id CouchTranslation sec_ct_set
global sec_ctimages rtplaninfo ct_MV
global Edepth Dfoc Theta3 DdepthRef
global original_isoc_pos

%datasource=1 - read from text file
%datsource=2 - already in global variables
firstimagefound=0;
if dvhonly<2
    readTomodata(field_width);
    
    n_of_proj=size(sinogram,2);
    
    
    axes(handles.axes1);
    colormap(jet(128));
    imagesc(sinogram);
    wb=waitbar(0,'reading CT data');
    %maxsin=max(max(sinogram));
    %subimage((128/maxsin)*sinogram,jet(128));
    delta_y=(pitch*10*field_width)/(51.0*num_of_subprojections);  %couch move per projection %GST Change %10* SJT change
    %halflength=0.5*delta_y*(n_of_proj-1)*num_of_subprojections; %GST Change  -1 is SJT change
    
    ncpoints=rtplaninfo.BeamSequence.Item_1.NumberOfControlPoints; %should equal n_of_proj
    itemncp=ItemText(ncpoints);
    final_isoc_pos=rtplaninfo.BeamSequence.Item_1.ControlPointSequence.(itemncp).IsocenterPosition;
    halflength=(isoc_pos(3)-final_isoc_pos(3))/2;  %mm
    matcen_y=(isoc_pos(3)+final_isoc_pos(3))/2;  %added 11/9/13 to replace the previous C ofM
    
    %determine which slices we want
    ct_ylist=ones(dose_dimensiony,1);
    if (halflength>0)
        ct_ymin=start_y-dose_gridy*(dose_dimensiony-1)/2 - halflength; %head first
    else
        ct_ymin=start_y+dose_gridy*(dose_dimensiony-1)/2 - halflength; %feet first
    end

    %dir_struct = dir(fullfile(ct_path,'*.dcm'));
    %sorted_names = sortrows({dir_struct.name}');
    nimages=size(ct_sorted_names,1);
    % couchpos=zeros(nimages,1); %commented out since now done when reading CT
    % for i=1:nimages
    %     ctinfo = dicominfo([ct_path cell2mat(sorted_names(i))]);
    %     couchpos(i)=ctinfo.ImagePositionPatient(3); %mm
    %     waitbar(i/(2*nimages));
    % end
    % [sorted_couch,sorted_index]=sortrows(couchpos); %index is index to sorted_names
    pixelsize=ctpix;
    
    if datasource==2
        %round ct_ymin to nearest ctslice
        %ctinfo = dicominfo([ct_path ct_filename]);
        couch1=ct_sorted_couch(1); %mm
        ctstep=abs(couch1-ct_sorted_couch(2));
        couch2=ctstep*floor(0.5 + couch1/ctstep);
        dcouch=couch1-couch2; %mm if couch position not multiple of step
        
        ct_ratio=floor(0.5 + (ct_ymin-(dcouch))/(ctstep));
        ct_temp=ct_ratio*ctstep +(dcouch);
        ct_rounding=ct_ymin-ct_temp;
        ct_ymin=ct_temp;
        
    end
    
    
    for i=1:dose_dimensiony
        if halflength>0
            ct_ylist(i)=ct_ymin+(i-1)*dose_gridy; %head first
        else
            ct_ylist(i)=ct_ymin-(i-1)*dose_gridy; %feet first
        end
    end
    ct_yindex=(dose_dimensiony+1)/2; %index for display
    ct_zindex=(dose_dimensionxz+1)/2; %GST - these were originally set to the same as ct_yindex - not sure if these should be set to a different value?
    ct_xindex=(dose_dimensionxz+1)/2; %GST - these were originally set to the same as ct_yindex - not sure if these should be set to a different value?
    %identify and load ct slices
    
    gotCTimage=zeros(1,dose_dimensiony);
    for i=1:nimages
        couch_long=ct_sorted_couch(i);
        for j=1:dose_dimensiony
            if abs(couch_long-ct_ylist(j))<0.01
                ctinfo = dicominfo([ct_path cell2mat(ct_sorted_names(ct_sorted_index(i)))]);
                ctimage = int16(dicomread(ctinfo))+ctinfo.RescaleIntercept;
                if ~firstimagefound
                    ctimages=zeros(dose_dimensiony,ctinfo.Width,ctinfo.Height);
                    firstimagefound=1;
                end
                ctimages(j,:,:)=ctimage;
                gotCTimage(j)=1;
                break
            end
        end
        waitbar(i/nimages);
        if min(gotCTimage)==1
            break
        end
    end
    if ct_MV
        nfound=sum(gotCTimage);
        if (nfound<nimages)
            errortxt=sprintf('only using %d of %d MV images\n this may cause problems with DVH',nfound,nimages);
            choice = questdlg(errortxt,'too short a y dimension', ...
                'Use anyway','try larger','try larger');
            if ~strcmp(choice,'Use anyway')
                close(wb);
                return
            end
        end
    end
    %rotx=double(1+ctinfo.Width)/2.0;
    %rotz=double(1+ctinfo.Height)/2.0;
    %rotx=rotx+isoc_pos(1)/pixelsize;  %added 29/6/11
    %rotz=rotz+isoc_pos(2)/pixelsize;  %added to cope with isoc shifting patient in v4
    %new version of this added 25/1/13
    ctImPosPat=ctinfo.ImagePositionPatient;
    ctpatientposition=ctinfo.PatientPosition; %eg HFS
    %swap coords if not HFS
    %these have not been tested for FF patients
    %therefore almost certainly wrong for FFS and FFP
    if strncmp(ctpatientposition,'HFP',3)
        ctImPosPat(1)=-ctinfo.ImagePositionPatient(1);
        ctImPosPat(2)=-ctinfo.ImagePositionPatient(2);
        isoc_pos(2)=-isoc_pos(2); %added 23/4/14
    end
    if strncmp(ctpatientposition,'FFS',3)
        ctImPosPat(1)=-ctinfo.ImagePositionPatient(1);
        ctImPosPat(3)=-ctinfo.ImagePositionPatient(3);
    end
    if strncmp(ctpatientposition,'FFP',3)
        ctImPosPat(3)=-ctinfo.ImagePositionPatient(3);
        ctImPosPat(2)=-ctinfo.ImagePositionPatient(2);
        isoc_pos(2)=-isoc_pos(2); %added 23/4/14
    end
    rotx=1-(ctImPosPat(1)-isoc_pos(1))/pixelsize; %added 25/1/13
    rotz=1-(ctImPosPat(2)-isoc_pos(2))/pixelsize; %added 25/1/13
    close(wb);
%     if strncmpi(rtplaninfo.RTPlanName,'cheese',6)
%         %following fix for cheese DQA only
%         rotx=rotx-rtplaninfo.PatientSetupSequence.Item_1.SetupDeviceSequence.Item_1.SetupDeviceParameter/pixelsize;%IEC X red laser
%         rotz=rotz-rtplaninfo.PatientSetupSequence.Item_1.SetupDeviceSequence.Item_3.SetupDeviceParameter/pixelsize;;%IEC Z red laser
%     end
    
    if sec_ct_set
        wb=waitbar(0,'secondary CT data');
        sec_couch1=sec_ct_sorted_couch(1); %mm
        sec_ctstep=abs(sec_couch1-sec_ct_sorted_couch(2));
        %will accept a match within half of ct_step for secondary
        matchmax=0.5*(sec_ctstep) + 0.001;
        gotsecCTimage=zeros(1,dose_dimensiony);
        nsecimages=size(sec_ct_sorted_names,1);
        for i=1:nsecimages
            couch_long=sec_ct_sorted_couch(i);
            for j=1:dose_dimensiony
                if abs(couch_long-ct_ylist(j))<matchmax
                    ctinfo = dicominfo([sec_ct_path cell2mat(sec_ct_sorted_names(sec_ct_sorted_index(i)))]);
                    ctimage = int16(dicomread(ctinfo))+ctinfo.RescaleIntercept;
                    sec_ctimages(j,:,:)=imrotate((ctimage+1000),TomoRoll,'bilinear','crop')-1000;%rotate counterclockwise by TomoRoll
                    gotsecCTimage(j)=1;
                    break
                end
            end
            waitbar(i/nsecimages);
            if min(gotsecCTimage)==1
                break
            end
        end
        close(wb);
        if min(gotsecCTimage)==0
            errordlg('unable to find CT for each slice','secondary CT set');
            return
        end
        %mask to primary
        MVcent=ctImPosPat;
        
        ctarsize=size(ctimages);
        ct_imdim=ctarsize(2); %256+ or 512
        MVradius=ctpix*(ct_imdim-1)/2;
        shifted_MVcent= MVcent;
        shifted_MVcent(1)= shifted_MVcent(1)+ CouchTranslation(1);
        shifted_MVcent(2)= shifted_MVcent(2)+ CouchTranslation(3);


        
        shifted_MVcent=shifted_MVcent+MVradius;
        kvImPos=ctinfo.ImagePositionPatient;
        sec_ct_Rows=ctinfo.Rows; %kv not necessarily square if couch added
        sec_ct_Columns=ctinfo.Columns;
        kvpix=ctinfo.PixelSpacing(1);
        kvcent(2)=0.5+(shifted_MVcent(2)-kvImPos(2))/kvpix;
        kvcent(1)=0.5+(shifted_MVcent(1)-kvImPos(1))/kvpix; %corrds in kv pixel space of MV centre
        [Xgrid,Ygrid] = meshgrid(1:double(sec_ct_Columns),1:double(sec_ct_Rows));
        maskmat=zeros(sec_ct_Rows, sec_ct_Columns);
        radiusmat=kvpix*sqrt((Xgrid-kvcent(1)).^2 +(Ygrid-kvcent(2)).^2); %radius in mm
        maskmat(radiusmat<MVradius)=1;
        for j=1:dose_dimensiony
            if(gotCTimage(j)) %only mask on slices where MV image exists
                ct_temp=squeeze(sec_ctimages(j,:,:));
                ct_temp(maskmat==1)=-1000;
                sec_ctimages(j,:,:)=ct_temp;
            end
        end
        
        
    else
        gotsecCTimage=0;
    end
    
    if (min(gotCTimage)+min(gotsecCTimage))==0
        errordlg('unable to find CT for each slice. Please either reduce grid set y (slice spacing) or dimension y (number of slices). You may need to select LockYdim','CT set too short');
        return
    end
    
    
    if min(gotCTimage)==0 %missing images, but OK since secondary exist
        for j=1:dose_dimensiony
            if ~gotCTimage(j)
                ctimages(j,:,:)=-1000*ones(ct_imdim,ct_imdim);
            end
        end
    end %this is needed to keep dimensions OK
    
    
    
    %set up cube for dose calc
    totalpoints=dose_dimensionxz*dose_dimensionxz*dose_dimensiony;
    Xvalue=zeros(1,totalpoints);
    Yvalue=zeros(1,totalpoints);
    Zvalue=zeros(1,totalpoints);
    n_ones=ones(1,totalpoints);
    n_ones2D=ones(1,dose_dimensionxz*dose_dimensionxz);
    
    %create 2D list of points for CT
    Xvalue2D=zeros(1,dose_dimensionxz*dose_dimensionxz);
    Zvalue2D=zeros(1,dose_dimensionxz*dose_dimensionxz);
    n=0;
    for i=1:dose_dimensionxz
        xzlist(i)=(i-(dose_dimensionxz+1)/2)*dose_gridxz;
    end
    
    for i=1:dose_dimensionxz
        for k=1:dose_dimensionxz
            n=n+1;
            Xvalue2D(n)=matcen_x+xzlist(i);
            Zvalue2D(n)=matcen_z+xzlist(k);
        end
    end
    
    snapct=get(handles.radiobutton18,'Value');
    if sec_ct_set&&~snapct
        %snapping to dose cube to match Marina's calculation
        marinaMVcent=dosecubeinfo.ImagePositionPatient + MVradius;
        dosegridpix=dosecubeinfo.PixelSpacing(1);
        dosegrid_dimensionx=dosecubeinfo.Rows;
        x_step=double(dosegridpix);
        x_pos = marinaMVcent(1) - isoc_pos(1);
        %X = (ix - 1)*pix + imposx % MR: from 3.2devB
        % MR: assuming ix is the centre of the image, ix = dose_dimensionx/2
        x_offset = x_pos + dosegridpix/2 - dosegridpix*(double(dosegrid_dimensionx))/2;
        xlist = x_offset : x_step : (x_offset + x_step*double(dosegrid_dimensionx - 1));
        %
        y_step = double(dosegridpix);
        y_pos = marinaMVcent(2) - isoc_pos(2);
        %Z = (1 - iz)*pix - imposz % MR: from 3.2devB
        % MR: assuming iz is the centre of the image, ix = dose_dimensionz/2
        %y_pos = -1*y_pos % MR invert the Y-axis. Why?
        y_offset = y_pos + dosegridpix/2 - dosegridpix*(double(dosegrid_dimensionx))/2;
        ylist = y_offset : y_step : (y_offset + y_step*double(dosegrid_dimensionx - 1));
        ylist = -1*ylist;
        nlist=max(size(xlist));
        nvalue=max(size(xzlist));
        difmatx=zeros(nlist,nvalue);
        difmaty=zeros(nlist,nvalue);
        for ilist=1:nlist
            for ivalue=1:nvalue
                difmatx(ilist,ivalue)=xlist(ilist)-matcen_x-xzlist(ivalue);
                difmaty(ilist,ivalue)=ylist(ilist)-matcen_z-xzlist(ivalue);
            end
        end
        minxp=min(min(difmatx(difmatx>=0))); %nearest positive to zero
        minxn=max(max(difmatx(difmatx<=0))); %nearest negative to zero
        minyp=min(min(difmaty(difmaty>=0))); %nearest positive to zero
        minyn=max(max(difmaty(difmaty<=0))); %nearest negative to zero    
        if abs(minxp)<abs(minxn)
            Xvalue2D=Xvalue2D+minxp;
        else
            Xvalue2D=Xvalue2D+minxn;
        end
        if abs(minyp)<abs(minyn)
            Zvalue2D=Zvalue2D+minyp;
        else
            Zvalue2D=Zvalue2D+minyn;
        end
            
    end
        
    
    endalpha=0:(51*num_of_subprojections-1); %GST change
    fifty_one_ones=ones(1,51*num_of_subprojections);
    
    gantry_start=original_gantry_start+(180/(51*num_of_subprojections)); %GST
    gantry_start=gantry_start-TomoRoll; %added 1/8/13
    alphav=(gantry_start*pi/180)+endalpha*(2*pi/(51*num_of_subprojections)); %negative since Tomo goes clockwise
    alphav(alphav<0)=alphav(alphav<0)+2*pi;
    %calculate 51 focus points
    dxf=850*sin(alphav); %850mm
    dzf=850*cos(alphav);
    xfocus=(dxf)/pixelsize + rotx;
    zfocus=(-dzf)/pixelsize + rotz; %coords of all foci in pixel space - outside the CT image
    np2d=dose_dimensionxz*dose_dimensionxz;
    
    alphav=n_ones2D'*alphav; %turn into matrix
    [phi, rho]=cart2pol(Xvalue2D,Zvalue2D);
    
    phi=phi'*fifty_one_ones;
    rho=rho'*fifty_one_ones;
    delta_depth=rho.*cos((pi/2)-alphav-phi);
    ppp=rho.*sin((pi/2)-alphav-phi);
    theta=atan(ppp./(850-delta_depth));
    %calculate distance from focus and effective depth for each calcualtion
    %point at each of 51 gantry angle in each CT slice
    
    xpixel=(Xvalue2D)/pixelsize + rotx;
    zpixel=(-Zvalue2D)/pixelsize + rotz; %coords of all points in pixels
    if sec_ct_set
        sec_rotx=1-(kvImPos(1)-original_isoc_pos(1))/sec_ctpix; %using original isoc poistion since moves apply only in Mv, not kV
        sec_rotz=1-(kvImPos(2)-original_isoc_pos(2))/sec_ctpix;
        if strncmpi(rtplaninfo.RTPlanName,'cheese',6)
            %following fix for cheese DQA only
            sec_rotx=sec_rotx-rtplaninfo.PatientSetupSequence.Item_1.SetupDeviceSequence.Item_1.SetupDeviceParameter/pixelsize;%IEC X red laser
            sec_rotz=sec_rotz-rtplaninfo.PatientSetupSequence.Item_1.SetupDeviceSequence.Item_3.SetupDeviceParameter/pixelsize;;%IEC Z red laser
        end
        sec_xfocus=(dxf)/sec_ctpix + sec_rotx;
        sec_zfocus=(-dzf)/sec_ctpix + sec_rotz; %coords of all foci in pixel space - outside the CT image
        sec_xpixel=(Xvalue2D)/sec_ctpix + sec_rotx;
        sec_zpixel=(-Zvalue2D)/sec_ctpix + sec_rotz; %coords of all points in pixels
    end
    dfromfoc=zeros(dose_dimensionxz*dose_dimensionxz,51*num_of_subprojections); %preallocate
    effdepth=zeros(dose_dimensiony,dose_dimensionxz*dose_dimensionxz,51*num_of_subprojections);
    wb=waitbar(0,'ray tracing');
    ngctimages=ctimages;  %non-global copies
    ngsec_ct_set=sec_ct_set;
    ngsec_ctimages=sec_ctimages;
    ngnum_of_subprojections=num_of_subprojections;
    ngsec_ctpix=sec_ctpix;
    %matlabpool(2);
    for iang=1:51*num_of_subprojections
        dx=xpixel-xfocus(iang);
        dz=zpixel-zfocus(iang);
        dfromfoc(:,iang)=sqrt(dx.*dx + dz.*dz)*(pixelsize); %distance in mm
        effdepth(:,:,iang)=radpath(ngctimages,xfocus(iang),xpixel,zfocus(iang),zpixel,ct_MV)*(pixelsize);
        if ngsec_ct_set
            temp_kvdepth=radpath(ngsec_ctimages,sec_xfocus(iang),sec_xpixel,sec_zfocus(iang),sec_zpixel,0)*(ngsec_ctpix);
            effdepth(:,:,iang)=effdepth(:,:,iang)+temp_kvdepth;

        end
        waitbar(iang/51*ngnum_of_subprojections);
    end
    if (halflength<0) %works for one FFS patient - need to check if generally works
        matcen_y=-matcen_y;
        set(handles.edit16,'String',-2*matcen_y);
    end
    %if snapping to Ct slices, change matcen_y
    snapct=get(handles.radiobutton18,'Value');
    if snapct
        %round matcen_y to nearest ctslice
        couch1=ct_sorted_couch(1); %mm
        ctstep=abs(couch1-ct_sorted_couch(2));
        couch2=ctstep*floor(0.5 + couch1/ctstep);
        dcouch=couch1-couch2; %mm if couch position not multiple of step
        ct_ratio=floor(0.5 + (matcen_y-(dcouch))/(ctstep));
        ct_temp=ct_ratio*ctstep +(dcouch);
        matcen_rounding=ct_ymin-ct_temp;
        matcen_y=ct_temp;
    end
    
    
    
    
    %put dfromfoc and effdepth into 51 x n^3 vectors
    n=0;
    ndim=dose_dimensionxz*dose_dimensionxz*dose_dimensiony;
    Edepth=zeros(ndim,51*num_of_subprojections);
    Dfoc=zeros(ndim,51*num_of_subprojections);
    Theta3=zeros(ndim,51*num_of_subprojections);
    
    for j=1:dose_dimensiony
        n2d=0;
        for i=1:dose_dimensionxz
            for k=1:dose_dimensionxz
                n=n+1;
                n2d=n2d+1;
                Xvalue(n)=Xvalue2D(n2d);
                Yvalue(n)=ct_ylist(j)-start_y;
                Zvalue(n)=Zvalue2D(n2d);
                Edepth(n,:)=effdepth(j,n2d,:);
                Dfoc(n,:)=dfromfoc(n2d,:);
                Theta3(n,:)=theta(n2d,:);
            end
        end
    end
    
    YvalueRef=Yvalue;
    DdepthRef=repmat(delta_depth,dose_dimensiony,1);
    dosecube=0;
    close(wb);
else %skip to this point if second half of dvhcalc
    delta_y=(pitch*10*field_width)/(51.0*num_of_subprojections);  %couch move per projection %GST Change %10* SJT change
    %halflength=0.5*delta_y*(n_of_proj-1)*num_of_subprojections; %GST Change  -1 is SJT change
    
    ncpoints=rtplaninfo.BeamSequence.Item_1.NumberOfControlPoints; %should equal n_of_proj
    itemncp=ItemText(ncpoints);
    final_isoc_pos=rtplaninfo.BeamSequence.Item_1.ControlPointSequence.(itemncp).IsocenterPosition;
    halflength=(isoc_pos(3)-final_isoc_pos(3))/2;  %mm
    Yvalue=YvalueRef;
 
end 
if dvhonly==1
    return
end
gantryindex=0;
Ddepth=DdepthRef;
if amp_4d>0 %Varian RPM 4D analysis
    rpm_dat2=rpm_data(min_index4d:max_index4d,:);
    rpm_mark2=rpm_mark(min_index4d:max_index4d);
    selectP=rpm_dat2(strcmp(rpm_mark2,'P'),:);
    meanP=mean(selectP(:,1)); %mean peak
    selectZ=rpm_dat2(strcmp(rpm_mark2,'Z'),:);
    meanZ=mean(selectZ(:,1)); %mean zenith
    meanval=mean(rpm_dat2(:,1)); %average position
    oldamp=(meanP-meanZ)/2;
    rpm_dat2(:,4)=(rpm_dat2(:,1)-meanval)*(amp_4d/oldamp);  %rescale to amplitude
    %find first point with chosen phase
    iphase=1;
     while(1)
         if (rpm_dat2(iphase,2)<=phase_4d) && (rpm_dat2(iphase+1,2)>=phase_4d) 
             break;
         end
         iphase=iphase+1;
 
     end %iphase is index of first value in array with chosen phase
     time4d=rpm_dat2(iphase,3);
     start_time=time4d;
     max_time=rpm_data(max_index4d,3);
     proj_time=1000.0 * gantry_period/51; %time taken per projection
     subproj_time=proj_time/num_of_subprojections;
     sum4d=0; n4d=0;  %for monitoring of code - can remove when issues sorted
end   
    

wb=waitbar(0,'calculating dose');
%Yvalue=Yvalue+halflength;
for p=1:n_of_proj
    waitbar(p/n_of_proj);
    
    %Split Projection Into Subprojections; allows for odd values up to 11. GST
    whole_projection=sinogram(:,p);
        subprojections=zeros(64,num_of_subprojections);
        for MLCindex=1:64
            if whole_projection(MLCindex)<=(1/num_of_subprojections)
                subprojections(MLCindex,(num_of_subprojections+1)/2)=whole_projection(MLCindex);
            elseif whole_projection(MLCindex)<=(3/num_of_subprojections)
                subprojections(MLCindex,(num_of_subprojections+1)/2)=1/num_of_subprojections;
                subprojections(MLCindex,(num_of_subprojections+3)/2)=(whole_projection(MLCindex)-1/num_of_subprojections)/2;
                subprojections(MLCindex,(num_of_subprojections-1)/2)=(whole_projection(MLCindex)-1/num_of_subprojections)/2;
            elseif whole_projection(MLCindex)<=(5/num_of_subprojections)
                subprojections(MLCindex,((num_of_subprojections-1)/2):((num_of_subprojections+3)/2))=1/num_of_subprojections;
                subprojections(MLCindex,(num_of_subprojections+5)/2)=(whole_projection(MLCindex)-3/num_of_subprojections)/2;
                subprojections(MLCindex,(num_of_subprojections-3)/2)=(whole_projection(MLCindex)-3/num_of_subprojections)/2;
            elseif whole_projection(MLCindex)<=(7/num_of_subprojections)
                subprojections(MLCindex,((num_of_subprojections-3)/2):((num_of_subprojections+5)/2))=1/num_of_subprojections;
                subprojections(MLCindex,(num_of_subprojections+7)/2)=(whole_projection(MLCindex)-5/num_of_subprojections)/2;
                subprojections(MLCindex,(num_of_subprojections-5)/2)=(whole_projection(MLCindex)-5/num_of_subprojections)/2;
            elseif whole_projection(MLCindex)<=(9/num_of_subprojections)
                subprojections(MLCindex,((num_of_subprojections-5)/2):((num_of_subprojections+7)/2))=1/num_of_subprojections;
                subprojections(MLCindex,(num_of_subprojections+9)/2)=(whole_projection(MLCindex)-7/num_of_subprojections)/2;
                subprojections(MLCindex,(num_of_subprojections-7)/2)=(whole_projection(MLCindex)-7/num_of_subprojections)/2;    
            elseif whole_projection(MLCindex)<=(11/num_of_subprojections)
                subprojections(MLCindex,((num_of_subprojections-7)/2):((num_of_subprojections+9)/2))=1/num_of_subprojections;
                subprojections(MLCindex,(num_of_subprojections+11)/2)=(whole_projection(MLCindex)-9/num_of_subprojections)/2;
                subprojections(MLCindex,(num_of_subprojections-9)/2)=(whole_projection(MLCindex)-9/num_of_subprojections)/2; 
            end
        end

    for subprojindex=1:num_of_subprojections
        if amp_4d>0 % 4d correction required
            y4d=interp1(rpm_dat2(:,3),rpm_dat2(:,4),time4d); %this requires input file to use mm
            time4d=time4d+subproj_time; %GST change
            if (time4d>max_time)
                time4d=time4d+start_time-max_time; %restart trace
            end
            sum4d=sum4d+y4d; n4d=n4d+1;
        else
            y4d=0;
        end
    
    
        gantryindex=gantryindex+1;
        if gantryindex>51*num_of_subprojections
            gantryindex=1;
        end
        projection=subprojections(:,subprojindex);
        n=segmentprojection();
        if (n>0)  %negative y4d 28/9/11  need to decide on sign
            dosecube=dosecube+dose_from_projection(Yvalue-y4d,Theta3(:,gantryindex),Edepth(:,gantryindex),Dfoc(:,gantryindex),gantry_period,0.1.*Ddepth(:,gantryindex));
        end
        if (halflength>0)
            Yvalue=Yvalue+delta_y; %head first
        else
            Yvalue=Yvalue-delta_y; %feet first
        end

    end
end
close(wb);
%mean4d=sum4d/n4d;
mindepth=min(Edepth');
dosecube(mindepth<1.5)=-0.001;  %want small negative to flag outside patient
%set to 1.5mm here.  5mm (or other user-set value) to be applied later
if ~dvhonly
    %calculate limit for plots
    plotlim1=max(max(max(dosecube)));
    dose_scale=dosecubeinfo.DoseGridScaling/nfractions;
    plotlim2=dose_scale*double(max(max(max(dicomdosecube))));
    plotlimit=max(plotlim1,plotlim2);
    
    %display dose in planes
    plot_transverse(handles);
    ystart=TomoPositionPatient(3);
    yquote=(ct_yindex-(dose_dimensiony+1)/2)*dose_gridy + matcen_y;
    indextemp=1+ (ystart-(yquote+y_offset_gamma))/dosecubeinfo.SliceThickness;
    dosecubeindex(3)=int16(indextemp);
    
    xzdim=dosecubeinfo.Rows;
    pix=dosecubeinfo.PixelSpacing(1);
    %xzmid=double(xzdim)/2+0.5;  %eg 128.5 for 256
    %following lines correct for offsets of plan and dose cube
    imposx=TomoPositionPatient(1)-isoc_pos(1);
    imposz=TomoPositionPatient(2)-isoc_pos(2);
    xmid=1-imposx/pix;
    zmid=1-imposz/pix;
    
    dosecubeindex(2)=floor(0.5+ xmid+Xvalue(dose_dimensionxz*(ct_xindex-1)+1)*1/pix);
    dosecubeindex(1)=floor(0.5+ zmid-Zvalue(ct_zindex)*1/pix);  %z goes backwards
    

    plot_dicom(handles);
    perctol=dose_tol*100;  %for log file
    if comparetwo
        calc_ddtwo(handles);
    else
        %calc_gamma_index(handles);
        if strcmpi(DefIndex,'Gamma')
            calc_gamma_wendling(handles);
        elseif strcmpi(DefIndex,'Box')
            calc_boxindex(handles);
        elseif strcmp(Defindex,'Kappa')
            calc_kappa(handles);
        end
    end
end