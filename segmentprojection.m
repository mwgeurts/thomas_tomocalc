function [nproj, segments] = segmentprojection(projection)
% segmentprojection converts a single projection into a set of MLC segments 
% as described in Thomas et al. Independent dose calculation software for 
% tomotherapy, Med Phys 2015; 39: 160-167. This function contains two
% nested functions, split_projection and segmentsubproj, that are called
% recursively. These functions have been slightly modified from the  
% original CheckTomo function by Simon Thomas to remove its dependency on 
% global variables.
%
% The following inputs are required for execution:
%   projection: 64 x 1 vector of leaf open fractions for each leaf in the
%       projection to be segmented
%
% The following variables are returned upon successful completion:
%   nproj: number of segments in the segments structure
%   segments: structure containing the fields nseg (integer number of
%       segments), startval (1 x nseg vector containing starting leaf of 
%       each segment), endval (1 x nseg vector  contianing ending leaf of 
%       each segment), and value (1 x nseg vector containing leaf open time 
%       fraction of each segment)
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

% Initialize nested function shared variables
segments.nseg = 0;
subproj = [];
lensubproj = [];

% Execute split_projection to identify how many segments are needed
nproj = split_projection(projection, 64, 0, 0);

% If at least one projection is needed
if nproj > 0

    % Initialize segment counter
    projval = 1;

    % While the counter is less than the total number of segments needed
    while projval <= nproj
        
        % Execute segmentsubproj to add segment details to the return
        % structure. If more segments are found to be needed, nproj will
        % increase and extend the while loop
        nproj = segmentsubproj(projval, nproj);
        
        % Invemement the segment counter
        projval = projval + 1;
    end
end

% Clear temporary variables
clear subproj lensubproj projval projection;

% Nested function split_projection. This returns the total number of neeed
% projections as well as storing the subprojection details in the shared
% variables lensubproj and subproj. This function essentially identifies
% adjacent leaves
function nproj = split_projection(proj, nleaf, nproj, offsetval)
    
    % Initialize flag
    inproj = false;
    
    % Loop through each leaf
    for leaf = 1:nleaf
        
        % If the leaf open time is nonzero
        if proj(leaf) > 0
            
            % If the flag is set (by a previous adjacent open leaf)
            if inproj
                
                % Increment the length of this segment
                lensubproj(nproj) = lensubproj(nproj) + 1;
                
                % Add a new endval for this segment
                subproj(nproj, lensubproj(nproj), 1) = leaf + offsetval;
                
                % Add a new open fraction based on this leaf
                subproj(nproj, lensubproj(nproj), 2) = proj(leaf);  
                
            % Otherwise, all previous leaves have been closed
            else 
                
                % Increment the projection count
                nproj = nproj + 1;
                
                % Store a new endval for this segment
                subproj(nproj, 1, 1) = leaf + offsetval;
                
                % Store an intitial endval
                subproj(nproj, 1, 2) = proj(leaf);
                
                % Initialize the length of this segment
                lensubproj(nproj) = 1;
                
                % Update the flag to indicate that an open leaf was found
                inproj = true;
            end
            
        % Otherwise, if this leaf's open time is zero
        else
            
            % Update the flag to indicate that the previous segment ended
            inproj = false;
        end
    end
    
    % Clear temporary variables
    clear inproj proj nleaf offsetval;
end 

% Nested function segmentsubproj. This function updates the segments return
% structure and returns a new number of segments by either re-running
% split_projection (if additional leaf open time exists) or ending the
% while loop above
function nproj = segmentsubproj(projval, nproj)

    % Increment the index of segments for the return structure field
    % vectors
    segments.nseg = segments.nseg + 1;
    
    % Add the starting leaf for this segment to the return structure
    segments.startval(segments.nseg) = subproj(projval, 1, 1);
    
    % Add the ending leaf for this segment to the return structure
    segments.endval(segments.nseg) = ...
        subproj(projval, lensubproj(projval), 1);
    
    % Add the smallest leaf open time fraction to the return structure. 
    % This is the open time for the base segment
    segments.value(segments.nseg) = ...
        min(subproj(projval, 1:lensubproj(projval), 2));
    
    % Subtract the base leaf open time from the remaining leaf open times
    subproj(projval, :, 2) = subproj(projval, :, 2) - ...
        segments.value(segments.nseg);
    
    % If the remaining value is greater than 0.01%, continue by
    % re-splitting and repeating this function
    if max(subproj(projval, :, 2)) > 0.0001 
        
        % Re-run split_projection on the remaining leaf open fractions 
        nproj = split_projection(subproj(projval, :, 2), 1 + ...
            segments.endval(segments.nseg) - ...
            segments.startval(segments.nseg), nproj, ...
            subproj(projval, 1, 1) - 1);
    end
    
    % Clear temporary variables
    clear projval;
end

% End the parent function
end