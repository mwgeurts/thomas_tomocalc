function dose = CheckTomoDose(varargin)






% Start by setting globals. These are used by the CheckTomo specific
% functions below, and have been left as globals to maintain compatibility
% with future versions
global start_path SP TPR OARXOPEN OARXLEAVES OARY segments ...
    reference_doserate dose_dimensiony;

% Persistently store parallel pool
persistent pool;

% If no inputs provided, return calcdose flag
if nargin == 1
    
    % Store the plan variable
    plan = varargin{1};
    
% Otherwise, store image and plan
elseif nargin == 2
    
    % Store image and plan variables
    image = varargin{1};
    plan = varargin{2};
    
% Otherwise, if 3 inputs are provided, store image, plan, and parallel pool
elseif nargin == 2
    
    % Store image and plan variables
    image = varargin{1};
    plan = varargin{2};
    pool = varargin{3};

% If zero or more than three arguments passed, log error
else
    if exist('Event', 'file') == 2
        Event(['An incorrect number of input arguments were passed to', ...
            ' CheckTomoDose'], 'ERROR');
    else
        error(['An incorrect number of input arguments were passed to', ...
            ' CheckTomoDose']);
    end
end

%% Verify registration
% If no registration vector was provided, add an empty one
if ~isfield(plan, 'registration')
    plan.registration = [0 0 0 0 0 0];
end

% Throw an error if the image registration pitch or yaw values are non-zero
if plan.registration(1) ~= 0 || plan.registration(2) ~= 0
    if exist('Event', 'file') == 2
        Event(['Dose calculation cannot handle pitch or yaw ', ...
            'corrections at this time'], 'ERROR');
    else
        error(['Dose calculation cannot handle pitch or yaw ', ...
            'corrections at this time']);
    end
end

% Log beginning of dose calculation and start timer(s)
if exist('Event', 'file') == 2
    
    % If using a parallel pool
    if isobject(pool) && pool.Connected
        Event(sprintf(['Beginning dose calculation using parallel ', ...
            'pool with %i workers'], pool.NumWorkers));
        tic
        ticBytes(pool)
    else
        Event('Beginning dose calculation a single worker');
        tic
    end
end

% Execute in try/catch statement
try 
    
%% Initialize parallel pool
% Attach the necessary files, if they do not already exist
    

    
%% Define dose_from_sin() variables

reference_doserate=8.5;

field_width = sum(abs([plan.frontField plan.backField]));

[start_path, ~, ~] = fileparts(mfilename('fullpath'));

sinogram = plan.agnostic;

pitch = plan.pitch;

num_of_subprojections = 1;

% Loop through the events cell array
for i = 1:size(plan.events, 1)
    
    % If type is isoX, apply IECX registration adjustment
    if strcmp(plan.events{i,2}, 'isoX')
        
        % Define isoc_pos X position
        isoc_pos(1) = (plan.events{i,3} - plan.registration(4))*10;

    % Otherwise, if type is isoY, apply IECZ registration adjustment
    elseif strcmp(plan.events{i,2}, 'isoY')
        
        % Define isoc_pos Y position
        isoc_pos(2) = (plan.events{i,3} + plan.registration(6))*10;
        
    % Otherwise, if type is isoZ, apply IECY registration adjustment
    elseif strcmp(plan.events{i,2}, 'isoZ')
        
        % Define isoc_pos Z position
        isoc_pos(3) = (plan.events{i,3} - plan.registration(5))*10;
        start_y = isoc_pos(3);
        
    % Otherwise, if type is gantryAngle, apply roll registration adjustment
    elseif strcmp(plan.events{i,2}, 'gantryAngle')
        
        % Define original_gantry_start and TomoRoll
        original_gantry_start = plan.events{i,3};
        TomoRoll = plan.registration(3) * 180/pi;
    end
end

% Define dose dimensions
dose_dimensionxz=image.dimensions(2)/2;
dose_dimensiony=image.dimensions(3);
dose_gridxz=image.width(1)*2*10;
dose_gridy=image.width(3)*10;

ctpix=image.width(1)*10;

matcen_x = 0;
matcen_z = 0;

ctImPosPat = image.start*10;

ct_MV = false;

%% Run dose_from_sin code

readTomodata(field_width);

n_of_proj=size(sinogram,2);

delta_y=(pitch*10*field_width)/(51.0*num_of_subprojections);

ncpoints=size(sinogram,2);

final_isoc_pos=[isoc_pos(1) isoc_pos(2) (isoc_pos(3) - ncpoints/51 * ...
    pitch * field_width)];

halflength=(isoc_pos(3)-final_isoc_pos(3))/2;  %mm

matcen_y=(isoc_pos(3)+final_isoc_pos(3))/2;  %added 11/9/13 to replace the previous C ofM

%determine which slices we want
ct_ylist=ones(dose_dimensiony,1);
if (halflength>0)
    ct_ymin=start_y-dose_gridy*(dose_dimensiony-1)/2 - halflength; %head first
else
    ct_ymin=start_y+dose_gridy*(dose_dimensiony-1)/2 - halflength; %feet first
end

nimages=image.dimensions(3);

pixelsize=ctpix;

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

% ctimages is the CT image data converted to 
ctimages = interp1(image.ivdt(:,1), image.ivdt(:,2), ...
    permute(image.data, [3 1 2]), 'linear', 'extrap');

rotx=1-(ctImPosPat(1)-isoc_pos(1))/pixelsize; %added 25/1/13
rotz=1-(ctImPosPat(2)-isoc_pos(2))/pixelsize; %added 25/1/13
    

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

dfromfoc=zeros(dose_dimensionxz*dose_dimensionxz,51*num_of_subprojections); %preallocate
effdepth=zeros(dose_dimensiony,dose_dimensionxz*dose_dimensionxz,51*num_of_subprojections);

%% Ray Tracing
ngctimages=ctimages;  %non-global copies
ngnum_of_subprojections=num_of_subprojections;
%matlabpool(2);
for iang=1:51*num_of_subprojections
    dx=xpixel-xfocus(iang);
    dz=zpixel-zfocus(iang);
    dfromfoc(:,iang)=sqrt(dx.*dx + dz.*dz)*(pixelsize); %distance in mm
    effdepth(:,:,iang)=radpath(ngctimages,xfocus(iang),xpixel,zfocus(iang),zpixel,ct_MV)*(pixelsize);
end
if (halflength<0) %works for one FFS patient - need to check if generally works
    matcen_y=-matcen_y;
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
    
gantryindex=0;
Ddepth=DdepthRef;

for p=1:n_of_proj
    
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
        y4d=0;

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








    
    
    

% Catch errors, log, and rethrow
catch err
    if exist('Event', 'file') == 2
        Event(getReport(err, 'extended', 'hyperlinks', 'off'), 'ERROR');
    else
        rethrow(err);
    end
    
    % Return empty dose
    dose = [];
end
