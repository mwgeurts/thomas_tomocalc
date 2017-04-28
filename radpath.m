function  rpath = radpath(ctimages, xf, xp, zf, zp)
% radpath computes the radiation path length (in units g-mm/cc) for each
% voxel in a dose volume to a focal point specified in the input variables 
% xp and zp. This function is called once for each subprojection by
% CheckTomoDose to compute the radiation path length for every dose voxel
% to every subprojection. This function has been slightly modified from the
% original radpath CheckTomo function by Simon Thomas to remove the
% electron density conversion and set the dose IEC Y dimension to the CT 
% dimension; see below documentation for more information.
%
% The following inputs are required for execution:
%   ctimages: a 3D array of CT image densities (in units g/cc), where
%       the first dimension is the IEC Y, the second dimension is IEC Z,
%       and the third dimension is IEC X.
%   xf: vector of IEC X positions of the center of each dose voxel, in mm
%   xp: IEC X position of the focal source, in mm
%   zf: vector of IEC Z positions of the center of each dose voxel, in mm
%   zp: IEC Z position of the focal source, in mm
%
% The following variables are returned upon successful completion:
%   rpath: a 2D array of radiation path lengths, in g-mm/cc, where the
%       first dimension is the number of CT images and the second dimension
%       is each dose voxel (corresponding to each voxel defined in xf and 
%       zf)
%
% Author: Simon Thomas, adapted by Mark Geurts, mark.w.geurts@gmail.com 
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

% Compute the difference between the focal position and dose voxels in the
% IEC X and Z directions
dx = xf - xp;
dz = zf - zp;

% Calculate the maximum number of steps needed to trace from the most 
% distant dose voxel to the focal position, one voxel at a time
steps = max([abs(dx) abs(dz)]);

% Initialize empty matrix for return variable
rpath = zeros(size(ctimages,1), length(dx));

% Store the largest CT dimension in the IEC X and Z directions as a two
% element vector. This is necessary because the sub2ind function used in
% the for loop below requires a square indexing system. The CT does not
% necessarily need to be square (there are checks below)
ctsize = repmat(max(size(ctimages,2), size(ctimages,3)), 1, 2);

% Loop through each step
for i=1:steps
    
    % Convert the distances into voxel indices for IEC X and Z
    ix = floor(xp);   
    iz = floor(zp);
    
    % Flag all indices less than 1 by setting the index to 1
    ix(ix < 1) = 1;
    iz(iz < 1) = 1;
    
    % Flag all indices larger than the CT dimensions by setting to 1
    ix(ix > size(ctimages, 2)) = 1;
    iz(iz > size(ctimages, 3)) = 1;
    
    % If a flag has been reached (on either side of the CT), stop the loop
    if max(iz)==1 || max(ix)==1
        break;
    end
    
    % Otherwise, retrieve the CT density value at this step using the ix
    % and iz indices and add it to the dose array. Note that this assumes 
    % that the ctimages array already contains density values and not 
    % Hounsfield units. This is different from the original version of this 
    % function.
    rpath = rpath + ctimages(:, sub2ind(ctsize, iz, ix));  

    % Increment the IEC X and Z position arrays by one step to prepare for
    % the next loop
    xp = xp + dx / steps;
    zp = zp + dz / steps;
end

% Multiply the density array by the length of each step (dr). This converts
% the units of the return variable to g-mm/cc
rpath = rpath .* repmat(sqrt(dx.^2 + dz.^2)/steps, size(ctimages,1), 1);
