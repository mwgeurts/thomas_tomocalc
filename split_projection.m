function subprojections = split_projection(sinogram, p, n)
% split_projection will separate a single leaf open time for each MLC leaf
% in a sinogram projection into an array of leaf open times for each 
% subprojection. This function has been created from a code block
% originally part of the CheckTomo function dose_from_sin. It was separated
% as it is used in both the parallel (parfor) and single-threaded options
% in the function CheckTomoDose.
%
% The following inputs are required for execution:
%   sinogram: 64 x # array containing leaf open time fractions (0 to 1) for
%       each leaf and projection in the helical plan, where # is the number
%       of projections.
%   p: current projection to split (index 1 to #)
%   n: number of subprojections to split the projection into (odd integer
%       between 1 and 11).
%
% The following variables are returned upon successful completion:
%   subprojections: 64 x n array of leaf open fractions for each
%       sub-projection of projection p
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
% You should have received a copy of the GNU General Public License a

% Initialize subprojection leaf open time array
subprojections = zeros(64, n);

% Loop through each leaf
for MLCindex = 1:64

    % If the leaf open time fits within one subprojection
    if sinogram(MLCindex, p) <= (1 / n)

        % Store the entire leaf open time in the central subprojection
        subprojections(MLCindex, (n + 1) / 2) = sinogram(MLCindex, p);

    % Otherwise, if the leaf open time is less than three
    % subprojections
    elseif sinogram(MLCindex, p) <= (3 / n)

        % Store a full leaf open time in the central subprojection
        subprojections(MLCindex, (n + 1) / 2) = 1 / n;

        % Store half of the remaining time into the upper and lower
        % subprojections
        subprojections(MLCindex, (n + 3) / 2) = ...
            (sinogram(MLCindex, p) - 1 / n) / 2;
        subprojections(MLCindex, (n - 1) / 2) = ...
            (sinogram(MLCindex, p) - 1 / n) / 2;

    % Otherwise, if the leaf open time is less than five subprojections
    elseif sinogram(MLCindex, p) <= (5 / n)

        % Store a full leaf open time in the central subprojections
        subprojections(MLCindex, ((n - 1) / 2):((n + 3) / 2)) = 1/n;

        % Store half of the remaining time into the adjacent 
        % upper and lower subprojections
        subprojections(MLCindex, (n + 5) / 2) = ...
            (sinogram(MLCindex, p) - 3 / n) / 2;
        subprojections(MLCindex , (n - 3) / 2) = ...
            (sinogram(MLCindex, p) - 3 / n) / 2;

    % Otherwise, if the leaf open time is less than seven 
    % subprojections 
    elseif sinogram(MLCindex, p) <= (7 / n)

        % Store a full leaf open time in the central subprojections
        subprojections(MLCindex, ((n - 3) / 2):((n + 5) / 2)) = 1 / n;

        % Store half of the remaining time into the adjacent 
        % upper and lower subprojections
        subprojections(MLCindex , (n + 7) / 2) = ...
            (sinogram(MLCindex, p) - 5 / n) / 2;
        subprojections(MLCindex , (n - 5) / 2) = ...
            (sinogram(MLCindex, p) - 5 / n) / 2;

    % Otherwise, if the leaf open time is less than nine subprojections 
    elseif sinogram(MLCindex, p) <= (9 / n)

        % Store a full leaf open time in the central subprojections
        subprojections(MLCindex, ((n - 5) / 2):((n + 7) / 2)) = 1 / n;

        % Store half of the remaining time into the adjacent 
        % upper and lower subprojections
        subprojections(MLCindex , (n + 9) / 2) = ...
            (sinogram(MLCindex, p) - 7 / n) / 2;
        subprojections(MLCindex , (n-7) / 2) = ...
            (sinogram(MLCindex, p) - 7 / n) / 2;    

    % Otherwise, if the leaf open time is less than eleven 
    % subprojections 
    elseif sinogram(MLCindex, p) <= (11 / n)

        % Store a full leaf open time in the central subprojections
        subprojections(MLCindex, ((n - 7) / 2):((n + 9) / 2)) = 1 / n;

        % Store half of the remaining time into the adjacent 
        % upper and lower subprojections
        subprojections(MLCindex , (n + 11) / 2) = ...
            (sinogram(MLCindex, p) - 9 / n) / 2;
        subprojections(MLCindex , (n - 9) / 2) = ...
            (sinogram(MLCindex, p) - 9 / n) / 2; 
    end
end

% Clear temporary variables 
clear n p sinogram MLCindex;