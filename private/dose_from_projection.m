function dose = dose_from_projection(Yvalue, Xtheta, depth, dfromfoc, ...
    gantry_period, ddepth, reference_doserate, SP, TPR, OARXOPEN, ...
    OARXLEAVES, OARY, segments)
% dose_from_projection computes the 3D dose from a TomoTherapy projection.
% It is called by CheckTomoDose for each subprojection in a sinogram. This 
% function has been slightly modified from the original CheckTomo function 
% by Simon Thomas to remove its dependency on global variables by passing
% the global variables as additional input arguments.
%
% The following inputs are required for execution:
%   Yvalue: n x 1 vector of IEC Y positions, in cm, for each of n dose 
%       voxels to be computed. The length of Yvalue must equal Xtheta,
%       depth, dfromfoc, and ddepth
%   Xtheta: n x 1 vector of off axis angles, in radians, for each dose
%       voxel to be computed.
%   depth: n x 1 vector of effective depths, in g/cm^2, for each dose voxel
%       to be computed.
%   dfromfoc: n x 1 vector of distances from the focal source to each dose 
%       voxel to be computed. 
%   gantry_period: helical plan gantry period, in seconds.
%   ddepth: n x 1 vector of depths, in cm, for each dose voxel to be
%       computed.
%   reference_doserate: treatment system dose rate, in Gy/min
%   SP: structure with fields length, nsizes, sizes, and values. See the
%       function readTomodata for more details on each field.
%   TPR: structure with fields ndepths, depths, nsizes, sizes, and tpr. See
%       the function readTomodata for more details on each field.
%   OARXOPEN: structure with fields ndepths, depths, nvalues, indices, and
%       oar. See the function readTomodata for more details on each field.
%   OARXLEAVES: structure with fields nwidths, widths, ndepths, depths,
%       nvalues, indices, and oar. See the function readTomodata for more 
%       details on each field.
%   OARY: structure with fields nwidths, widths, ndepths, depths, nvalues, 
%       indices, and oar. See the function readTomodata for more details on 
%       each field.
%   segments: structure containing segment details, with fields nseg, 
%       startval, endval, and value. See segmentprojection for details on
%       each value.
%
% The following variables are returned upon successful completion:
%   dose: 1 x n vector of dose values, in Gy, for each voxel defined in 
%       Yvalue, depth, dfromfoc, and Xtheta for the provided segments. Note
%       that the dose vector is transposed from its respective input
%       vectors.
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

% Initialize return vector
dose = zeros(1, size(Yvalue, 2));

% Cap upper and lower limits of OARY depth values. Note that this assumes
% that the OARXOPEN and OARXLEAVES depth arrays have the same max depth.
depth = max(OARY.depths(1), min(OARY.depths(OARY.ndepths), depth));

% Loop through each segment in the subprojection
for seg = 1:segments.nseg
    
    % Store the number of leaves for this segment
    nleaves = 1 + segments.endval(seg) - segments.startval(seg);
    
    % Compute the MLC field width of this segment, using an effective leaf 
    % width at 85 cm SAD of 0.625 cm
    width = single(nleaves * 0.625);
    
    %% Calculate the TPR, Sp, and open field MLC and longitudinal OARs
    % Compute the TPR values of each voxel depth, computing the equivalent
    % square as 2 * SP.length * width / (SP.length + width). The method was 
    % modified from nearest in the original CheckTomo code to nearest to 
    % improve computation time
    tprvalue = interp2(TPR.sizes, TPR.depths, TPR.tpr, ...
        2 * SP.length * width / (SP.length + width), depth, 'nearest');
    
    % Calculate the output factor using the MLC field width
    spvalue = interp1(SP.sizes, SP.values, width, 'linear', 'extrap');
    
    % Compute the open field MLC OAR values for each voxel looking up the
    % radial theta (indexed to the nearest 0.005 radian). The interpolation
    % method was modified from nearest in the original CheckTomo code to
    % improve computation time
    oar40value = interp2(OARXOPEN.depths,OARXOPEN.indices,OARXOPEN.oar, ...
        depth, min(abs(Xtheta) / 0.005, max(OARXOPEN.indices)), 'nearest');
    
    % Compute the longitudinal OAR value for each voxel, indexing the 
    % radial theta to the nearest 0.001 radian and limiting the maximum 
    % value to the maximum index available in OARY. The interpolation 
    % method was modified from nearest in the original CheckTomo code 
    % to improve computation time
    oaryvalue = interp3(OARY.depths, OARY.indices, OARY.widths, ...
        OARY.oar, depth', min(max(OARY.indices), atan(abs(Yvalue ./ ...
        (85 - ddepth'))) / 0.001), single(repmat(max(width, ...
        OARY.widths(1)), 1, length(Xtheta))), 'nearest');

    %% Calculate the MLC defined OAR
    % Store the radial off-axis angle for each voxel relative to the center
    % of the MLC defined field
    theta = abs(Xtheta - atan((0.5 * (segments.startval(seg) + ...
        segments.endval(seg)) - 32.5) * 0.625 / 85));
       
    % If the number of open leaves exceeds the second to last index in
    % OARXLEAVES
    if nleaves > OARXLEAVES.widths(OARXLEAVES.nwidths-1)
        
        % Adjust the off axis angle to be relative to the maxium leaf width
        % stored in OARXLEAVES (4).
        theta = max(0, theta - atan((nleaves - max(OARXLEAVES.widths)) * ...
            0.5 * 0.625 / 85.0));
        
        % Set the OARXLEAVES widths index to its maximum value 
        nleaves = OARXLEAVES.nwidths;
    end
    
    % Compute the MLC defined field OAR values for each voxel, indexing the 
    % radial distance to the nearest 0.001 radian and limiting the maximum 
    % value to the maximum index available in OARXLEAVES. The interpolation 
    % method was modified from nearest in the original CheckTomo code 
    % to improve computation time
    oarxvalue = interp2(OARXLEAVES.depths, OARXLEAVES.indices, ...
        OARXLEAVES.oar(:, :, nleaves), depth, min(max(OARXLEAVES.indices), ...
        theta / 0.001), 'nearest');
    
    %% Compute total output factor
    % Compute the total output factor for each voxel from this segment as 
    % the product of the leaf open fraction, output factor, and three OAR 
    % factors, summing with all other segments. This is converted to dose
    % later.
    dose = dose + segments.value(seg) * spvalue * tprvalue' .* ...
        (oarxvalue .* oar40value)' .* oaryvalue;
end

% Convert to dose by multiplying the output factors by an inverse square 
% factor (85 cm/depth)^2, the machine reference dose rate (Gy/min),
% reference dose measurement inverse square correction (86.5/85)^2, 
% reference measurement TPR (1.37), minute to second conversion (1/60), and
% total projection time (gantry period/51)
dose = dose .* (85 ./ (dfromfoc .* cos(Xtheta)))'.^2 * reference_doserate ...
    * (86.5/85)^2 / (1.37 * 60) * gantry_period / 51;

% Clear temporary variables
clear oarxvalue nleaves theta oaryvalue oar40value spvalue tprvalue width ...
    Yvalue Xtheta depth dfromfoc gantry_period ddepth reference_doserate ...
    SP TPR OARXOPEN OARXLEAVES OARY segments seg;