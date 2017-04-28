function dose = CheckTomoDose(varargin)






% Start by setting globals. These are used by the CheckTomo specific
% functions below, and have been left as globals to maintain compatibility
% with future versions
global start_path SP TPR OARXOPEN OARXLEAVES OARY segments ...
    reference_doserate dose_dimensiony projection subproj lensubproj;

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

% Clear varargin
clear varargin;

%% Verify plan type
% Throw an error if the plan type is non-Helical
if ~isfield(plan, 'planType') || ~strcmp(plan.planType, 'Helical')
    if exist('Event', 'file') == 2
        Event('This function can only calculate dose for helical plans', ...
            'ERROR');
    else
        error('This function can only calculate dose for helical plans');
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
% The following variables are used in the original dose_from_sin() function
% that was adapted into this function. 

% Log event
if exist('Event', 'file') == 2
    Event('Defining runtime variables');
end

% Downsample the dose grid by this factor in the IEC X and Z directions
downsample = 8;

% Define the machine nominal dose rate in Gy/min as a global. This is used
% by dose_from_projection()
reference_doserate = 8.5;

% Store the field width from the plan structure
field_width = sum(abs([plan.frontField plan.backField]));

% Store the path of this function as a global. This is used by 
% readTomodata()
[start_path, ~, ~] = fileparts(mfilename('fullpath'));

% Store the machine agnostic sinogram (normalized to 1)
sinogram = plan.agnostic;

% Store the plan pitch (unitless)
pitch = plan.pitch;

% Define the number of subprojections to calculate. If set to 1, then one
% ray tracing/dose calculation is performed for each of the 51 projections.
% If set to 3, then three positions are modeled (similar to the 
% supersampling flag in the TPS). As discussed below, this value must be an 
% odd integer between 1 and 11. Note that increasing this variable
% significantly increases memory requirements and computation time!
num_of_subprojections = 1;

% Loop through the events cell array
for i = 1:size(plan.events, 1)
    
    % If type is isoX, apply IECX registration adjustment
    if strcmp(plan.events{i,2}, 'isoX')
        
        % Define isoc_pos IECX position, in mm. This corresponds to the X
        % position of the DICOM CT header, plus any IECX offset provided by
        % the registration array.
        isoc_pos(1) = (plan.events{i,3} + plan.registration(4)) * 10;

    % Otherwise, if type is isoY, apply IECZ registration adjustment
    elseif strcmp(plan.events{i,2}, 'isoY')
        
        % Define isoc_pos IECZ position, in mm. This corresponds to the Z
        % position of the DICOM header, which is inverted from the
        % TomoTherapy coordinate system (TomoTherapy specifies the position
        % of the lower left voxel)
        isoc_pos(2) = -(plan.events{i,3} + plan.registration(6)) * 10;
        
    % Otherwise, if type is isoZ, apply IECY registration adjustment
    elseif strcmp(plan.events{i,2}, 'isoZ')
        
        % Define isoc_pos IECY position, in mm. This differs from the plan 
        % isoZ position (which is the IECY direction in the TomoTherapy 
        % coordinate system) by the number of empty projections (defined by 
        % startTrim) converted into couch travel, in mm.
        isoc_pos(3) = (plan.events{i,3} + plan.registration(5) ...
            + (plan.startTrim(1) - 1) / 51 * field_width * pitch) * -10;
        
        % Define start_y as the isoc_pos IECY position
        start_y = isoc_pos(3);
        
    % Otherwise, if type is gantryAngle, apply roll registration adjustment
    elseif strcmp(plan.events{i,2}, 'gantryAngle')
        
        % Define original_gantry_start angle and TomoRoll, in degrees. The
        % plan start gantry angle must be adjusted by the number of gantry
        % rotations for all leading empty projections (defined by
        % startTrim).For TomoRoll, the registration vector is in radians, 
        % so must be converted to degrees.
        original_gantry_start = mod(plan.events{i,3} + ...
            (plan.startTrim(1) - 1) * 360/ 51, 360);
        TomoRoll = plan.registration(3) * 180/pi;
    end
end

% Define dose dimensions
% The axial (xz) dimension assumes the dose volume will be square, and is
% set to the current IECZ dimension of the CT (typically the smaller of the
% two because of couch insertion) divided by the downsampling factor
dose_dimensionxz = floor(image.dimensions(2) / downsample);

% The y (IECY) dimension is equal to the CT image dimension. This means
% dose will be calculated at the slice resolution
dose_dimensiony = image.dimensions(3);

% dose_gridxz defines the voxel width (in mm) of the dose grid in the axial 
% direction. It is set based on the CT voxel size multiplied by the
% downsample factor converted to mm.
dose_gridxz = image.width(1) * downsample * 10;

% dose_gridy defines the voxel width (in mm) of the dose grid in the IECY
% dimension. It is set to the CT slice width in mm.
dose_gridy = image.width(3) * 10;

% Define ctpix as the CT axial resolution in mm.
ctpix = image.width(1) * 10;

% Define the dose grid center, in mm. This is set to be the center of the
% CT images in the IECX and IECZ directions, in DICOM coordinates
matcen_x = (image.start(1) + image.width(1) * ...
    (image.dimensions(1)-1)/2) * 10;
matcen_z = -(image.start(2) + image.width(2) * ...
    (image.dimensions(2)-1)/2) * 10;

% Define ctImPosPat as the three element DICOM position of the CT image. A
% coordinate transformation is needed to go from image.start (TomoTherapy
% coordinate system) which defines the lower left voxel, to  DICOM, in mm.
ctImPosPat = [image.start(1) -(image.start(2) + image.width(2)...
    * (image.dimensions(2)-1)) -image.start(3)] * 10;

% Define the gantry period in seconds, defined by the plan scale (seconds 
% per projection) multiplied by 51 projections and divided by the number of
% fractions. Gantry periods truncate to 0.1 sec
gantry_period = single(ceil(plan.scale * 51 / plan.fractions * 10) / 10);

% Log the values defined above
if exist('Event', 'file') == 2
    Event(sprintf('Downsampling Factor: %i', downsample));
    Event(sprintf('Dose rate: %0.1f Gy/min', reference_doserate));
    Event(sprintf('Field Width: %0.1f cm', field_width));
    Event(['Start path: ', start_path]);
    Event(sprintf('Pitch: %0.3f', pitch));
    Event(sprintf('Gantry Period: %0.3f sec', gantry_period));
    Event(sprintf('Number of projections: %i', size(sinogram, 2)));
    Event(sprintf('Number of subprojections: %i', num_of_subprojections));
    Event(sprintf('Isocenter position: [%0.2f %0.2f %0.2f] mm', isoc_pos(1), ...
        isoc_pos(2), isoc_pos(3)));
    Event(sprintf('Start Y: %0.2f mm', start_y));
    Event(sprintf('Gantry start angle: %0.2f deg', original_gantry_start));
    Event(sprintf('Roll correction: %0.2f deg', TomoRoll));
    Event(sprintf('Dose grid dimensions (xz): %i', dose_dimensionxz));
    Event(sprintf('Dose grid dimensions (y): %i', dose_dimensiony));
    Event(sprintf('Dose grid resolution (xz): %0.2f mm', dose_gridxz));
    Event(sprintf('Dose grid resolution (y): %0.2f mm', dose_gridy));
    Event(sprintf('Dose grid center (xz): [%0.2f %0.2f] mm', ...
        matcen_x, matcen_z));
    Event(sprintf('CT resolution: %0.2f mm', ctpix));
    Event(sprintf('CT image position: [%0.2f %0.2f %0.2f] mm', ctImPosPat(1), ...
        ctImPosPat(2), ctImPosPat(3)));
end

%% Run dose_from_sin code
% Log the start of readTomodata
if exist('Event', 'file') == 2
    Event(sprintf(['Executing readTomodata(%0.1f) to load TPR, OAR, ', ...
        'and Sp factors'], field_width));
end

% Run readTomodata to populate TPR, OARXOPEN, OARXLEAVES, OARY, and SP
% global variables
readTomodata(field_width);

% Define the number of projections in the sinogram
n_of_proj = size(sinogram, 2);

% Define the amount of couch (IECY) motion per subprojection
delta_y = (pitch * 10 * field_width) / (51 * num_of_subprojections);

% Define the number of control points (this is the same as the number of
% projections?)
ncpoints = size(sinogram, 2);

% Define the final isocenter position as the DICOM position at the last
% projection, in mm. The amount of couch travel (applied to isoc_pos(3) is
% determined using the product of the couch rotations, pitch, and field
% width.
final_isoc_pos = [isoc_pos(1) isoc_pos(2) (isoc_pos(3) - ncpoints / 51 * ...
    pitch * field_width * 10)];

% Define half length as half the distance between the start and final
% isocenter position
halflength = (isoc_pos(3) - final_isoc_pos(3)) / 2;

% Determine which slices we want
% Initialize ct_ylist as a unit vector with the same length as the number 
% dose voxels in the IECY dimension
ct_ylist = ones(dose_dimensiony, 1);

% Define ct_ymin using the head first notation, as the TomoTherapy image
% always will assume HFS
ct_ymin = start_y - dose_gridy * (dose_dimensiony - 1)/2 - halflength;

% Loop through each dose slice
for i=1:dose_dimensiony

    % Compute the CT slice that corresponds to the dose slice, assuming
    % head first, as the TomoTherapy image always assumes HFS
    ct_ylist(i) = ct_ymin + (i - 1) * dose_gridy;

end

% Define ctimages as the 3D CT image volume converted to physical density
% using the IVDT associated with the plan. The image must also be flipped
% around to get put back into DICOM coordinate space (which is what the
% code below expects)
ctimages = flip(interp1(image.ivdt(:,1), image.ivdt(:,2), ...
    permute(image.data, [3 2 1]), 'linear', 'extrap'), 1);

% If the ct data is not square, pad it
if size(ctimages,2) > size(ctimages,3)
  
    % Log action
    Event(sprintf('Padding CT in Z dimension to %i x %i', ...
        size(ctimages,2), size(ctimages,2)));
    
    % Add empty voxels
    ctimages = cat(3, zeros(size(ctimages,1), size(ctimages,2), ...
        size(ctimages,2) - size(ctimages,3)), ctimages);
    
elseif size(ctimages,2) < size(ctimages,3)
    
    % Log action
    Event(sprintf('Padding CT in X dimension to %i x %i', ...
        size(ctimages,3), size(ctimages,3)));
    
    % Add empty voxels
    ctimages = cat(2, ctimages, zeros(size(ctimages,1), size(ctimages,3) - ...
        size(ctimages,2), size(ctimages,3)));
end

% Define rotx and rotz as the central voxel of the CT image given relative
% to isocenter
rotx = 1 - (ctImPosPat(1) - isoc_pos(1)) / ctpix;
rotz = 1 - (ctImPosPat(2) - isoc_pos(2)) / ctpix;

% Set up cube for dose calc
% Define totalpoints as the total number of dose voxels in the 3D volume
totalpoints = dose_dimensionxz * dose_dimensionxz * dose_dimensiony;

% Log the total number of dose voxels
if exist('Event', 'file') == 2
    Event(sprintf('Total number of calculation points: %i', totalpoints));
end

% Initialize empty vectors for the IECX, Y, and Z values of each dose voxel
Xvalue=zeros(1, totalpoints);
Yvalue=zeros(1, totalpoints);
Zvalue=zeros(1, totalpoints);

% Create 2D list of points for CT
% Initialize empty vectors for the IECX and Z values of each dose voxel
Xvalue2D = zeros(1, dose_dimensionxz * dose_dimensionxz);
Zvalue2D = zeros(1, dose_dimensionxz * dose_dimensionxz);

% Initialize counter
n = 0;

% Initialize empty vector of dose voxel coordinates, assuming grid is
% centered
xzlist = zeros(1, dose_dimensionxz);

% Loop through the number of dose voxels
for i = 1:dose_dimensionxz
    
    % Store the voxel position, centered on zero
    xzlist(i) = (i - (dose_dimensionxz + 1) / 2) * dose_gridxz;
end

% Loop through the number of dose voxels in the X direction
for i=1:dose_dimensionxz
    
    % Loop through the number of dose voxels in the Z direction
    for k=1:dose_dimensionxz
        
        % Increment counter
        n = n + 1;
        
        % Store the dose voxel X position as the centered position plus the
        % dose matrix X offset
        Xvalue2D(n) = matcen_x + xzlist(i);
        
        % Store the dose voxel Z position as the centered position plus the
        % dose matric Z offset
        Zvalue2D(n) = matcen_z + xzlist(k);
    end
end

% Define the gantry start angle as the original gantry angle plus the
% half number of degrees for a single subprojection. This basically centers
% the gantry start angle on the first subprojection.
gantry_start = original_gantry_start + (180 / (51 * num_of_subprojections));

% Subtract any roll corrections to the start position
gantry_start = gantry_start - TomoRoll;

% Define a vector of angles, converted to radians, of each subprojection
alphav = (gantry_start * pi / 180) + (0:(51 * num_of_subprojections - 1)) * ...
    (2 * pi / (51 * num_of_subprojections));

% Correct any angles less than zero to be positive
alphav(alphav < 0) = alphav(alphav < 0) + 2 * pi;

% Compute the pixel number, starting from the lower left pixel, of each
% subprojection focal spot. Note that the IEC Z positions are flipped as
% the DICOM images are from the top down.
xfocus = 850 * sin(alphav) / ctpix + rotx;
zfocus = -850 * cos(alphav) / ctpix + rotz;

% Compute the polar coordinates of the IEC X and Z values of each dose
% voxel
[phi, rho] = cart2pol(Xvalue2D, Zvalue2D);

% Compute several factors that are used to compute the distance from the
% focal point and effective depth of each voxel
delta_depth = repmat(rho', 1, 51 * num_of_subprojections) .* ...
    cos((pi/2) - repmat(alphav, dose_dimensionxz * dose_dimensionxz, 1) - ...
    repmat(phi', 1, 51 * num_of_subprojections));
ppp = repmat(rho', 1, 51 * num_of_subprojections) .* ...
    sin((pi/2) - repmat(alphav, dose_dimensionxz * dose_dimensionxz, 1) - ...
    repmat(phi', 1, 51 * num_of_subprojections));
theta = atan(ppp ./ (850 - delta_depth));

%% Garbage collection
clear alphav ppp rho xzlist ct_ymin;

%% Ray Tracing
if exist('Event', 'file') == 2
    Event('Starting ray tracing steps');
end

% Store the corresponding CT pixel value of each dose voxel
xpixel = (Xvalue2D) / ctpix + rotx;
zpixel = (-Zvalue2D) / ctpix + rotz;

% Initialize empty arrays to store focal distances and effective depth
% vectors. The distances are a 2D matrix, with one column for each 
% subprojection and one row for each XZ dose voxel in a single image. The 
% effective depth is a 3D matrix, with the first dimension storing each IEC 
% Y dose image, the second each XZ dose voxel, and the third for each 
% subprojection.  These are single precision arrays, to reduce overall 
% memory requirements
dfromfoc = single(zeros(dose_dimensionxz * dose_dimensionxz, ...
    51 * num_of_subprojections));
effdepth = single(zeros(dose_dimensiony, dose_dimensionxz * ...
    dose_dimensionxz, 51 * num_of_subprojections));

% Loop through each subprojection
for iang = 1:51 * num_of_subprojections
    
    % Log subprojection
    if exist('Event', 'file') == 2
        Event(sprintf('Calculating effective depths from angle %i of %i', ...
            iang, 51*num_of_subprojections));
    end
    
    % Compute the distance from each voxel to the focal spot of this
    % subprojection, in mm
    dfromfoc(:,iang) = single(sqrt((xpixel - xfocus(iang)).^2 + ...
        (zpixel - zfocus(iang)).^2) * ctpix);
    
    % Compute the effective depth by executing radpath()
    effdepth(:,:,iang) = single(radpath(ctimages, xfocus(iang), xpixel,...
        zfocus(iang), zpixel) * ctpix);
end

%put dfromfoc and effdepth into 51 x n^3 vectors
n=0;
ndim=dose_dimensionxz * dose_dimensionxz * dose_dimensiony;
Edepth=single(zeros(ndim, 51 * num_of_subprojections));
Dfoc=single(zeros(ndim, 51 * num_of_subprojections));
Theta3=single(zeros(ndim, 51 * num_of_subprojections));
    
for j=1:dose_dimensiony
    n2d=0;
    for i=1:dose_dimensionxz
        for k=1:dose_dimensionxz
            n=n+1;
            n2d=n2d+1;
            Xvalue(n)=single(Xvalue2D(n2d));
            Yvalue(n)=single(ct_ylist(j)-start_y);
            Zvalue(n)=single(Zvalue2D(n2d));
            Edepth(n,:)=single(effdepth(j,n2d,:));
            Dfoc(n,:)=single(dfromfoc(n2d,:));
            Theta3(n,:)=single(theta(n2d,:));
        end
    end
end

%% Garbage Collection
% Clear array variables that are no longer used
clear ct_ylist ctimages dfromfoc effdepth theta xfocus xpixel ...
    Xvalue2D zfocus zpixel Zvalue2D ct_MV n n2d;

%% Dose Calculation
if exist('Event', 'file') == 2
    Event('Starting dose calculation steps');
end

DdepthRef=single(repmat(delta_depth,dose_dimensiony,1));
dosecube=0;
    
gantryindex=0;
Ddepth=DdepthRef;

for p=1:n_of_proj
    
    if exist('Event', 'file') == 2
        Event(sprintf('Calculating dose from projection %i of %i', ...
            p, n_of_proj));
    end
    
    %Split Projection Into Subprojections; allows for odd values up to 11. GST
    whole_projection=sinogram(:,p);
    subprojections=zeros(64,num_of_subprojections);
    for MLCindex=1:64
        if whole_projection(MLCindex)<=(1/num_of_subprojections)
            subprojections(MLCindex,(num_of_subprojections+1)/2)=...
                whole_projection(MLCindex);
        elseif whole_projection(MLCindex)<=(3/num_of_subprojections)
            subprojections(MLCindex,(num_of_subprojections+1)/2)=...
                1/num_of_subprojections;
            subprojections(MLCindex,(num_of_subprojections+3)/2)=...
                (whole_projection(MLCindex)-1/num_of_subprojections)/2;
            subprojections(MLCindex,(num_of_subprojections-1)/2)=...
                (whole_projection(MLCindex)-1/num_of_subprojections)/2;
        elseif whole_projection(MLCindex)<=(5/num_of_subprojections)
            subprojections(MLCindex,((num_of_subprojections-1)/2):...
                ((num_of_subprojections+3)/2))=1/num_of_subprojections;
            subprojections(MLCindex,(num_of_subprojections+5)/2)=...
                (whole_projection(MLCindex)-3/num_of_subprojections)/2;
            subprojections(MLCindex,(num_of_subprojections-3)/2)=...
                (whole_projection(MLCindex)-3/num_of_subprojections)/2;
        elseif whole_projection(MLCindex)<=(7/num_of_subprojections)
            subprojections(MLCindex,((num_of_subprojections-3)/2):...
                ((num_of_subprojections+5)/2))=1/num_of_subprojections;
            subprojections(MLCindex,(num_of_subprojections+7)/2)=...
                (whole_projection(MLCindex)-5/num_of_subprojections)/2;
            subprojections(MLCindex,(num_of_subprojections-5)/2)=...
                (whole_projection(MLCindex)-5/num_of_subprojections)/2;
        elseif whole_projection(MLCindex)<=(9/num_of_subprojections)
            subprojections(MLCindex,((num_of_subprojections-5)/2):...
                ((num_of_subprojections+7)/2))=1/num_of_subprojections;
            subprojections(MLCindex,(num_of_subprojections+9)/2)=...
                (whole_projection(MLCindex)-7/num_of_subprojections)/2;
            subprojections(MLCindex,(num_of_subprojections-7)/2)=...
                (whole_projection(MLCindex)-7/num_of_subprojections)/2;    
        elseif whole_projection(MLCindex)<=(11/num_of_subprojections)
            subprojections(MLCindex,((num_of_subprojections-7)/2):...
                ((num_of_subprojections+9)/2))=1/num_of_subprojections;
            subprojections(MLCindex,(num_of_subprojections+11)/2)=...
                (whole_projection(MLCindex)-9/num_of_subprojections)/2;
            subprojections(MLCindex,(num_of_subprojections-9)/2)=...
                (whole_projection(MLCindex)-9/num_of_subprojections)/2; 
        end
    end

    for subprojindex=1:num_of_subprojections

        gantryindex=gantryindex+1;
        if gantryindex>51*num_of_subprojections
            gantryindex=1;
        end
        projection=subprojections(:,subprojindex);
        n=segmentprojection();
        if (n>0)  %negative y4d 28/9/11  need to decide on sign
            dosecube=dosecube+dose_from_projection(Yvalue, ...
                Theta3(:,gantryindex), Edepth(:,gantryindex), ...
                Dfoc(:,gantryindex), gantry_period, ...
                0.1.*Ddepth(:,gantryindex));
        end
        Yvalue=Yvalue+delta_y; %head first


    end
end

mindepth=min(Edepth, [], 2);
dosecube(mindepth<1.5) = 0;

%% Garbage collection
clear Dfoc Edepth sinogram segments projection subproj SP OARXLEAVES ...
    OARXOPEN OARY Theta3 TPR mindepth;

%% Interpolate back to CT resolution
% Log interpolation step
if exist('Event', 'file') == 2
    Event(sprintf(['Upsampling calculated dose image by %i using ', ...
        'nearest neighbor interpolation'], downsample));
end
        
% Multiply by number of fractions to get total dose
dosecube = dosecube * plan.fractions;

% Initialize dose.data array
dose.data = zeros(image.dimensions);  

% Copy dose image start, width, and dimensions from CT image
dose.start = image.start;
dose.width = image.width;
dose.dimensions = image.dimensions;

x = reshape(Xvalue, dose_dimensionxz, dose_dimensionxz, dose_dimensiony) ...
    + dose_gridxz/2;
z = reshape(Zvalue, dose_dimensionxz, dose_dimensionxz, dose_dimensiony) ...
    + dose_gridxz/2;
dose_small = reshape(dosecube, dose_dimensionxz, dose_dimensionxz, ...
    dose_dimensiony);

% Compute target meshgrids
[my, mx] = meshgrid((dose.start(2)*10+dose.width*10*(dose.dimensions(2)-1))...
    :-dose.width(2)*10:dose.start(2)*10, ...
    dose.start(1)*10:dose.width(1)*10:...
    (dose.start(1)*10+dose.width*10*(dose.dimensions(1)-1)));

% Loop through slices
for i = 1:dose.dimensions(3)
   
    % Interpolate back to image dimensions
    dose.data(:,:,i) = interp2(x(:,:,i), z(:,:,i), ...
        dose_small(:,:,i), mx, my, 'nearest');
end

% Clear temporary variables
clear x y z i mx my dose_small dosecube;

% Log conclusion and stop timer
if exist('Event', 'file') == 2
    
    % If using a parallel pool
    if isobject(pool) && pool.Connected
        Event(sprintf(['Dose calculation completed successfully in %0.3f ', ...
            'seconds, sending %0.1f kB to %i workers'], toc, ...
            max(max(tocBytes(pool), pool.NumWorkers))));
    else
        Event(sprintf(['Dose calculation completed successfully in %0.3f ', ...
            'seconds'], toc));
    end
end

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
