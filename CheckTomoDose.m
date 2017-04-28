function dose = CheckTomoDose2(varargin)






% Start by setting globals. These are used by the CheckTomo specific
% functions below, and have been left as globals to maintain compatibility
% with future versions
global SP TPR OARXOPEN OARXLEAVES OARY segments ...
    reference_doserate projection subproj lensubproj;

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

%% Initialize parallel pool
% Attach the necessary files, if they do not already exist
if isobject(pool) && pool.Connected
    
    
    
end


%% Begin computations
% Execute in try/catch statement
try 
    
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
        
    % Otherwise, if type is gantryAngle, apply roll registration adjustment
    elseif strcmp(plan.events{i,2}, 'gantryAngle')
        
        % Define original_gantry_start angle adjusting for roll, in degrees. The
        % plan start gantry angle must be adjusted by the number of gantry
        % rotations for all leading empty projections (defined by
        % startTrim).
        original_gantry_start = mod(plan.events{i,3} + (plan.startTrim(1) ...
            - 1) * 360/51 + plan.registration(3) * 180/pi, 360);
    end
end

% Clear temporary variables
clear i;

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
    Event(sprintf('Pitch: %0.3f', pitch));
    Event(sprintf('Gantry Period: %0.3f sec', gantry_period));
    Event(sprintf('Number of projections: %i', size(sinogram, 2)));
    Event(sprintf('Number of subprojections: %i', num_of_subprojections));
    Event(sprintf('Isocenter position: [%0.2f %0.2f %0.2f] mm', isoc_pos(1), ...
        isoc_pos(2), isoc_pos(3)));
    Event(sprintf('Gantry start angle: %0.2f deg', original_gantry_start));
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
[TPR, SP, OARXOPEN, OARXLEAVES, OARY] = readTomodata(field_width);

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
% Log the total number of dose voxels
if exist('Event', 'file') == 2
    Event(sprintf('Total number of calculation points: %i', ...
        dose_dimensionxz * dose_dimensionxz * dose_dimensiony));
end

% Store the position of each dose voxel in the X and Z locations, centered
% around zero
xzlist = ((1:dose_dimensionxz) - (dose_dimensionxz + 1) / 2) * dose_gridxz;

% Store the dose voxel X position as the centered position plus the dose
% matrix X offset
Xvalue2D = reshape(repmat(matcen_x + xzlist, dose_dimensionxz, 1), 1, []);

% Store the dose voxel Z position as the centered position plus the dose
% matrix Z offset
Zvalue2D = reshape(repmat(matcen_z + xzlist, dose_dimensionxz, 1)', 1, []);

% Define the gantry start angle as the original gantry angle plus the
% half number of degrees for a single subprojection. This basically centers
% the gantry start angle on the first subprojection.
gantry_start = original_gantry_start + (180 / (51 * num_of_subprojections));

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
clear alphav phi ppp rho xzlist matcen_x matcen_z;

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

% Define the final isocenter position as the DICOM position at the last
% projection, in mm. The amount of couch travel (applied to isoc_pos(3) is
% determined using the product of the couch rotations, pitch, and field
% width.
final_isoc_pos = [isoc_pos(1) isoc_pos(2) (isoc_pos(3) - ...
    size(sinogram, 2) / 51 * pitch * field_width * 10)];

% Define ct_ymin using the head first notation, as the TomoTherapy image
% always will assume HFS
ct_ymin = isoc_pos(3) - dose_gridy * (dose_dimensiony - 1)/2 - ...
    (isoc_pos(3) - final_isoc_pos(3)) / 2;

% Compute the CT slice that corresponds to the dose slice, assuming
% head first, as the TomoTherapy image always assumes HFS
ct_ylist = ct_ymin + ((1:dose_dimensiony) - 1) * dose_gridy;

% Store the IECX, Y, and Z values of each dose voxel in a vector
Xvalue = single(repmat(Xvalue2D, 1, dose_dimensiony));
Yvalue = single(reshape(repmat((ct_ylist-isoc_pos(3)), dose_dimensionxz * ...
    dose_dimensionxz, 1), 1, []));
Zvalue = single(repmat(Zvalue2D, 1, dose_dimensiony));

% Reshape the radiation path lengths, focal distance, angle and depths to a 
% 2D array of voxel x subprojection. The length of Xvalue, Yvalue, and 
% Zvalue will equal the size of the first dimension.
Edepth = single(reshape(permute(effdepth, [2 1 3]), dose_dimensionxz * ...
    dose_dimensionxz * dose_dimensiony, 51 * num_of_subprojections));
Dfoc = single(repmat(dfromfoc, dose_dimensiony, 1));
Theta3 = single(repmat(theta, dose_dimensiony, 1));
Ddepth = single(repmat(delta_depth, dose_dimensiony, 1));

%% Garbage Collection
% Clear array variables that are no longer used
clear ct_ylist ctimages dfromfoc effdepth theta xfocus xpixel delta_depth...
   Xvalue2D zfocus zpixel Zvalue2D n2d ndim ctpix isoc_pos rotx rotz;

%% Dose Calculation
% Log start of dose calculation
if exist('Event', 'file') == 2
    Event('Starting dose calculation steps');
end

% Initialize empty dose cube array
dosecube = single(zeros(1, dose_dimensionxz * dose_dimensionxz * ...
    dose_dimensiony));

% Loop through each projection in the sinogram
for p = 1:size(sinogram, 2)
    
    % Log current projection number and total projections
    if exist('Event', 'file') == 2
        Event(sprintf('Calculating dose from projection %i of %i', ...
            p, size(sinogram, 2)));
    end
    
    %% Split Projection Into Subprojections
    % This code block will separate a single leaf open time for each MLC
    % leaf in a sinogram projection into an array of leaf open times for
    % each subprojection.
    
    % Initialize subprojection leaf open time array
    subprojections = zeros(64, num_of_subprojections);
    
    % Loop through each leaf
    for MLCindex = 1:64
        
        % If the leaf open time fits within one subprojection
        if sinogram(MLCindex, p) <= (1 / num_of_subprojections)
            
            % Store the entire leaf open time in the central subprojection
            subprojections(MLCindex, (num_of_subprojections + 1) / 2) = ...
                sinogram(MLCindex, p);
            
        % Otherwise, if the leaf open time is less than three
        % subprojections
        elseif sinogram(MLCindex, p) <= (3 / num_of_subprojections)
            
            % Store a full leaf open time in the central subprojection
            subprojections(MLCindex, (num_of_subprojections + 1) / 2) = ...
                1 / num_of_subprojections;
            
            % Store half of the remaining time into the upper and lower
            % subprojections
            subprojections(MLCindex, (num_of_subprojections + 3) / 2) = ...
                (sinogram(MLCindex, p) - 1 / num_of_subprojections) / 2;
            subprojections(MLCindex, (num_of_subprojections - 1) / 2) = ...
                (sinogram(MLCindex, p) - 1 / num_of_subprojections) / 2;
            
        % Otherwise, if the leaf open time is less than five subprojections
        elseif sinogram(MLCindex, p) <= (5 / num_of_subprojections)
            
            % Store a full leaf open time in the central subprojections
            subprojections(MLCindex, ((num_of_subprojections - 1) / 2):...
                ((num_of_subprojections + 3) / 2)) = ...
                1/num_of_subprojections;
            
            % Store half of the remaining time into the adjacent 
            % upper and lower subprojections
            subprojections(MLCindex, (num_of_subprojections + 5) / 2) = ...
                (sinogram(MLCindex, p) - 3 / num_of_subprojections) / 2;
            subprojections(MLCindex , (num_of_subprojections - 3) / 2) = ...
                (sinogram(MLCindex, p) - 3 / num_of_subprojections) / 2;
            
        % Otherwise, if the leaf open time is less than seven 
        % subprojections 
        elseif sinogram(MLCindex, p) <= (7 / num_of_subprojections)
            
            % Store a full leaf open time in the central subprojections
            subprojections(MLCindex, ((num_of_subprojections - 3) / 2):...
                ((num_of_subprojections + 5) / 2)) = ...
                1 / num_of_subprojections;
            
            % Store half of the remaining time into the adjacent 
            % upper and lower subprojections
            subprojections(MLCindex , (num_of_subprojections + 7) / 2) = ...
                (sinogram(MLCindex, p) - 5 / num_of_subprojections) / 2;
            subprojections(MLCindex , (num_of_subprojections - 5) / 2) = ...
                (sinogram(MLCindex, p) - 5 / num_of_subprojections) / 2;
        
        % Otherwise, if the leaf open time is less than nine subprojections 
        elseif sinogram(MLCindex, p) <= (9 / num_of_subprojections)
            
            % Store a full leaf open time in the central subprojections
            subprojections(MLCindex, ((num_of_subprojections - 5) / 2):...
                ((num_of_subprojections + 7) / 2)) = ...
                1 / num_of_subprojections;
            
            % Store half of the remaining time into the adjacent 
            % upper and lower subprojections
            subprojections(MLCindex , (num_of_subprojections + 9) / 2) = ...
                (sinogram(MLCindex, p) - 7 / num_of_subprojections) / 2;
            subprojections(MLCindex , (num_of_subprojections-7) / 2) = ...
                (sinogram(MLCindex, p) - 7 / num_of_subprojections) / 2;    
        
        % Otherwise, if the leaf open time is less than eleven 
        % subprojections 
        elseif sinogram(MLCindex, p) <= (11 / num_of_subprojections)
            
            % Store a full leaf open time in the central subprojections
            subprojections(MLCindex, ((num_of_subprojections - 7) / 2):...
                ((num_of_subprojections + 9) / 2)) = ...
                1 / num_of_subprojections;
            
            % Store half of the remaining time into the adjacent 
            % upper and lower subprojections
            subprojections(MLCindex , (num_of_subprojections + 11) / 2) = ...
                (sinogram(MLCindex, p) - 9 / num_of_subprojections) / 2;
            subprojections(MLCindex , (num_of_subprojections - 9) / 2) = ...
                (sinogram(MLCindex, p) - 9 / num_of_subprojections) / 2; 
        end
    end

    %% Loop through each subprojection, calculating dose
    for i = 1:num_of_subprojections
        
        % Compute the angular index corresponding to this subprojection
        gantryindex = mod((p-1) * num_of_subprojections + i - 1, ...
            51 * num_of_subprojections) + 1;
        
        projection=subprojections(:,i);

        if (segmentprojection()>0)  %negative y4d 28/9/11  need to decide on sign
            dosecube=dosecube+dose_from_projection(Yvalue, ...
                Theta3(:,gantryindex), Edepth(:,gantryindex), ...
                Dfoc(:,gantryindex), gantry_period, ...
                0.1.*Ddepth(:,gantryindex));
        end
        Yvalue=Yvalue+(pitch * 10 * field_width) / ...
            (51 * num_of_subprojections); %head first


    end
end

mindepth=min(Edepth, [], 2);
dosecube(mindepth<1.5) = 0;

%% Garbage collection
clear Dfoc Edepth sinogram segments projection subproj SP OARXLEAVES ...
    OARXOPEN OARY Theta3 TPR mindepth gantry_period pitch Ddepth Yvalue ...
    subprojections;

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
[my, mx] = meshgrid((dose.start(2)*10+dose.width*10*(dose.dimensions(2)))...
    :-dose.width(2)*10:(dose.start(2)+dose.width)*10, ...
    dose.start(1)*10:dose.width(1)*10:...
    (dose.start(1)*10+dose.width*10*(dose.dimensions(1)-1)));

% Loop through slices
for i = 1:dose.dimensions(3)
   
    % Interpolate back to image dimensions
    dose.data(:,:,i) = interp2(x(:,:,i), z(:,:,i), ...
        dose_small(:,:,i), mx, my, 'nearest');
end

% Clear temporary variables
clear x y z i mx my dose_small dosecube downsample Xvalue Zvalue;

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
