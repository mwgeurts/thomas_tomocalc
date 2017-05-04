function dose = CheckTomoDose(varargin)
% CheckTomoDose is a modified form of the original dose_from_sin function
% in the CheckTomo tool developed by Simon Thomas. It was adapted to
% calculate dose given input data provided in the format used by the 
% tomo_extract and dicom_tools GitHub repositories developed by Mark
% Geurts (https://github.com/mwgeurts). The function was also expanded to
% allow computation on a MATLAB cluster to improve computation time.
% Various other minor changes have been made to different lines of the
% code; refer to the comments in this and the other functions of this
% repository for details.
% 
% For more information on the methods employed in this tool, see Thomas et 
% al. Independent dose calculation software for tomotherapy, Med Phys 2015; 
% 39: 160-167. The original tool, CheckTomo, can be obtained through the
% GPL license by contacting Simon Thomas (refer to the correspondence
% address in the journal article referenced above for contact information).
%
% This function will report progress information by calling the function
% Event(message, flag) if it exists in the application path, where the
% message is a string and the flag is one of the following strings: 'INFO',
% 'WARN', or 'ERROR'. See the file at the following address for an example:
% https://github.com/mwgeurts/exit_detector/blob/master/Event.m
%
% This function can be executed with various combinations of input
% arguments. Upon first execution, at least two arguments are required, as
% shown in the following example, where image and plan are structures
% returned by the tomo_extract functions LoadImage and LoadPlan:
%
%   dose = CheckTomoDose(image, plan);
%
% The image is stored persistently, so after the first call, a second plan
% may be calculated with only one input argument:
%
%   dose = CheckTomoDose(plan);
%
% A parallel pool can also be passed in the third input argument (the first
% two must be image and plan):
%
%   dose = CheckTomoDose(image, plan, pool);
%
% Finally, additional configuration options can be passed as name/value
% pairs for input arguments 4 and on. The available options are 
% 'downsample', 'reference_doserate', and 'num_of_subprojections'. See the
% code below for detail on each option:
%
%   dose = CheckTomoDose(image, plan, pool, 'reference_doserate', 8.2);
%   dose = CheckTomoDose(image, plan, [], 'num_of_subprojections', 3);
%
% Upon successful completion, this function returns the structure dose,
% which contains the following fields:
%   start: 1 x 3 vector of position cooordinates, in cm, for the lower left
%       voxel. These coordinates are identical to image.start.
%   width: 1 x 3 vector of voxel dimensions, in cm
%   dimensions: 1 x 3 vector of the size of dose.data
%   data: 3D array of dose values, in Gy, of dimensions provided above
%   
% Below is complete example of how this function is used:
%
%   % Load image and plan data for the plan specified by the UID below
%   image = LoadImage('./path/', 'TomoPhant^^^^_patient.xml', ...
%       '1.2.826.0.1.3680043.2.200.1828118229.362.96568.1276');
%   plan = LoadPlan('./path/', 'TomoPhant^^^^_patient.xml', ...
%       '1.2.826.0.1.3680043.2.200.1828118229.362.96568.1276');
%
%   % Execute CheckTomoDose with 2 workers and apply a downsampling factor
%   pool = parpool(2);
%   dose = CheckTomoDose(image, plan, pool, 'downsample', 4);
%
% Author: Simon Thomas, adapted by Mark Geurts, mark.w.geurts@gmail.com
% Original work Copyright (C) 2011-15  Simon Thomas 
% Adapted work Copyright (C) 2017 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

% Persistently store parallel pool and image
persistent pool image;

% Downsample the dose grid by this factor in the IEC X and Z directions.
% This value can be overridden in the input arguments
downsample = 8;

% Define the machine nominal dose rate in Gy/min as a global. This is used
% by dose_from_projection(). This value can be overriden in the input
% arguments
reference_doserate = 8.5;

% Define the number of subprojections to calculate. If set to 1, then one
% ray tracing/dose calculation is performed for each of the 51 projections.
% If set to 3, then three positions are modeled (similar to the 
% supersampling flag in the TPS). As discussed below, this value must be an 
% odd integer between 1 and 11. Note that increasing this variable
% significantly increases memory requirements and computation time. This 
% value can be overriden in the input arguments
num_of_subprojections = 1;

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
elseif nargin >= 3
    
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

% Store remaining input arguments
for i = 4:2:nargin
    
    % If the user passed a downsample input argument, update it
    if strcmpi(varargin{i}, 'downsample')
        downsample = varargin{i+1};
    
    % Else, if the user passed a reference_doserate input argument
    elseif strcmpi(varargin{i}, 'reference_doserate')
        reference_doserate = varargin{i+1};
    
    % Else, if the user passed a num_of_subprojections input argument
    elseif strcmpi(varargin{i}, 'num_of_subprojections')
        num_of_subprojections = varargin{i+1};
    end   
end

% Clear varargin
clear i varargin;

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
    
    % Define necessary files for parallel computation
    necessaryFiles = {
        'dose_from_projection.m'
        'segmentprojection.m'
        'split_projection.m'
    };

    % Loop through files
    for i = 1:length(necessaryFiles)
        
        % If the file does not exist
        if max(strcmp(which(necessaryFiles{i}), pool.AttachedFiles)) == 0
            
            % Log attachment
            if exist('Event', 'file') == 2
                Event(['Attaching runtime file ', necessaryFiles{i}, ...
                    ' to pool']);
            end
                
            % Attach it
            addAttachedFiles(pool,necessaryFiles{i});
        else
            % Log attachment
            if exist('Event', 'file') == 2
                Event(['Runtime file ', necessaryFiles{i}, ...
                    ' already attached to pool']);
            end
        end
    end
    
    % Clear temporary variables
    clear i necessaryFiles;
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

% Store the field width from the plan structure
field_width = sum(abs([plan.frontField plan.backField]));

% Store the machine agnostic sinogram (normalized to 1)
sinogram = plan.agnostic;

% Store the plan pitch (unitless)
pitch = plan.pitch;

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

% If a parallel pool exists
if isobject(pool) && pool.Connected
    
    % Log start of dose calculation
    if exist('Event', 'file') == 2
        Event(['Executing calculations on parallel pool ', ...
            '(progress not displayed)']);
    end
    
    % Loop through each projection in the sinogram in parallel
    parfor p = 1:size(sinogram, 2)

        % Split Projection Into Subprojections
        subprojections = ...
            split_projection(sinogram, p, num_of_subprojections);

        %% Loop through each subprojection, calculating dose
        for i = 1:num_of_subprojections

            % Compute the angular index corresponding to this subprojection
            gantryindex = mod((p-1) * num_of_subprojections + i - 1, ...
                51 * num_of_subprojections) + 1;

            % Separate subprojection into segments using the function
            % segmentprojection. If all leaves are closed, n will be zero
            [n, segments] = segmentprojection(subprojections(:, i));

            % If at least one segment is found, continue to calculate dose for
            % this subprojection
            if (n > 0)

                % Compute local Y position in cm, considering couch motion
                localY = (Yvalue + ((p-1) * num_of_subprojections + i) * ...
                    (pitch * 10 * field_width) / ...
                    (51 * num_of_subprojections)) / 10;
                
                % Execute dose_from_projection, adding result to existing
                % dosecube array. Note that dose_from_projection uses cm, so
                % the depths are converted by diving by 10
                dosecube = dosecube + dose_from_projection(localY, ...
                    Theta3(:,gantryindex), Edepth(:,gantryindex)/10, ...
                    Dfoc(:,gantryindex)/10, gantry_period, ...
                    Ddepth(:,gantryindex)/10, reference_doserate, SP, TPR, ...
                    OARXOPEN, OARXLEAVES, OARY, segments); %#ok<PFBNS,PFRUS>
            end
        end
    end 
else
    % Loop through each projection in the sinogram
    for p = 1:size(sinogram, 2)

        % Log current projection number and total projections
        if exist('Event', 'file') == 2
            Event(sprintf('Calculating dose from projection %i of %i', ...
                p, size(sinogram, 2)));
        end

        % Split Projection Into Subprojections
        subprojections = ...
            split_projection(sinogram, p, num_of_subprojections);

        %% Loop through each subprojection, calculating dose
        for i = 1:num_of_subprojections

            % Compute the angular index corresponding to this subprojection
            gantryindex = mod((p-1) * num_of_subprojections + i - 1, ...
                51 * num_of_subprojections) + 1;

            % Separate subprojection into segments using the function
            % segmentprojection. If all leaves are closed, n will be zero
            [n, segments] = segmentprojection(subprojections(:, i));

            % If at least one segment is found, continue to calculate dose for
            % this subprojection
            if (n > 0)

                % Compute local Y position in cm, considering couch motion
                localY = (Yvalue + ((p-1) * num_of_subprojections + i) * ...
                    (pitch * 10 * field_width) / ...
                    (51 * num_of_subprojections)) / 10;
                
                % Execute dose_from_projection, adding result to existing
                % dosecube array. Note that dose_from_projection uses cm, so
                % the depths are converted by diving by 10
                dosecube = dosecube + dose_from_projection(localY, ...
                    Theta3(:,gantryindex), Edepth(:,gantryindex)/10, ...
                    Dfoc(:,gantryindex)/10, gantry_period, ...
                    Ddepth(:,gantryindex)/10, reference_doserate, SP, TPR, ...
                    OARXOPEN, OARXLEAVES, OARY, segments);
            end
        end
    end
end

% Compute the minimum effective depth (from any source angle) of each dose
% voxel. This is used to truncate voxels that are near or outside the
% patient; see the next line of code below.
mindepth = min(Edepth, [], 2);

% Truncate all dose voxels where the effective depth is less than 1.5 mm.
% This effectively removes dose outside the patient/couch
dosecube(mindepth < 1.5) = 0;

%% Garbage collection
% Clean up unused variables
clear Dfoc Edepth sinogram segments projection subproj SP OARXLEAVES ...
    OARXOPEN OARY Theta3 TPR mindepth gantry_period pitch Ddepth Yvalue ...
    subprojections localY;

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

% Copy dose image start, width, and dimensions from CT image. After
% interpolation, the dose array will have the same coordinate system.
dose.start = image.start;
dose.width = image.width;
dose.dimensions = image.dimensions;

% Reshape the dose voxel IEC X and Z position vectors into meshgrid format.
x = reshape(Xvalue, dose_dimensionxz, dose_dimensionxz, dose_dimensiony);
z = reshape(Zvalue + dose_gridxz, dose_dimensionxz, dose_dimensionxz, ...
    dose_dimensiony);
dose_small = reshape(dosecube, dose_dimensionxz, dose_dimensionxz, ...
    dose_dimensiony);

% Compute target meshgrids
[my, mx] = meshgrid((dose.start(2) + dose.width(2) * dose.dimensions(2))...
    :-dose.width(2):(dose.start(2) + dose.width(2)), dose.start(1):...
    dose.width(1):(dose.start(1) + dose.width(2) * (dose.dimensions(1) - 1)));

% Loop through slices
for i = 1:dose.dimensions(3)
   
    % Interpolate back to image dimensions
    dose.data(:,:,i) = interp2(x(:,:,i)/10, z(:,:,i)/10, ...
        dose_small(:,:,i), mx, my, 'nearest', 0);
end

% Clear temporary variables
clear x y z i mx my dose_small dosecube downsample Xvalue Zvalue;

% Log conclusion and stop timer
if exist('Event', 'file') == 2
    
    % If using a parallel pool
    if isobject(pool) && pool.Connected
        t = tocBytes(pool);
        Event(sprintf(['Dose calculation completed successfully in %0.3f ', ...
            'seconds, sending %0.1f MB to %i workers'], toc, ...
            sum(t(:,1))/1024^2, pool.NumWorkers));
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
