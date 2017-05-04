function [TPR, SP, OARXOPEN, OARXLEAVES, OARY] = readTomodata(field_width)
% readTomodata loads the dose calculation parameters (TPR, SP, and OARs)
% for a specified field width. This function has been slightly modified 
% from the original radpath CheckTomo function by Simon Thomas to remove
% the use of global variables. The directory of the .txt files, designated 
% by start_path, is now assumed to be the same directory as this function.
% See the below documentation for more information.
% 
% The following inputs are required for execution:
%   field_width: double representing the plan field width in cm. This
%       function selects the nominal field width that most closely matches
%       the provided value (1, 2.5, or 5)
%
% The following variables are returned upon successful completion:
%   TPR: structure with fields ndepths, depths, nsizes, sizes, and tpr. See
%       the documentation in the code below for more details on each field.
%   SP: structure with fields length, nsizes, sizes, and values. See the
%       documentation in the code below for more details on each field.
%   OARXOPEN: structure with fields ndepths, depths, nvalues, indices, and
%       oar. See the documentation in the code below for more details on 
%       each field.
%   OARXLEAVES: structure with fields nwidths, widths, ndepths, depths,
%       nvalues, indices, and oar. See the documentation in the code below 
%       for more details on each field.
%   OARY: structure with fields nwidths, widths, ndepths, depths, nvalues, 
%       indices, and oar. See the documentation in the code below for more 
%       details on each field.
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

% Store the path of this function. The code below will search in this same
% path for the field width specific input files.
[start_path, ~, ~] = fileparts(mfilename('fullpath'));

% Initialize empty return variables
TPR = struct;
SP = struct;
OARXOPEN = struct;
OARXLEAVES = struct;
OARY = struct;

% Store the input variable as the length field in Sp
SP.length = single(field_width); 

% Store the file names and nominal field widths for each jaw position
files = {
    5.0 'LA4data50.txt'
    2.5 'LA4data25.txt'
    1.0 'LA4data10.txt'
};

% Match the field width to the closest input file nominal field width
[~, i] = min(abs(cell2mat(files(:,1)) - field_width));

% Open a file handle to the file coresponding to the matching field width
fid = fopen(fullfile(start_path, files{i,2}));

% Loop through the file contents
while ~feof(fid)   

    % Retrieve next line, and search for a header tag
    tagvalue = regexp(fgetl(fid), '\[.+\]', 'match');
    
    % If no tag was found, continue
    if isempty(tagvalue)
        continue;
    
    % If the header is [TPRdata]
    elseif strcmp(tagvalue{1}, '[TPRdata]')

        % Store the number of TPR depths (integer)
        D = textscan(fid, '%d', 1);
        TPR.ndepths = D{1};

        % Store the depths as single (vector of length ndepths)
        D = textscan(fid, '%f', TPR.ndepths);
        TPR.depths = single(D{1});

        % Store the number of TPR column sizes (integer)
        D = textscan(fid, '%d', 1);
        TPR.nsizes = D{1};

        % Store the TPR column sizes (vector of length nsizes)
        D = textscan(fid, '%f', TPR.nsizes);
        TPR.sizes = single(D{1});
        
        % Store the TPR data (2D array of ndepths x nsizes)
        D = textscan(fid, ...
            strjoin(repmat({'%f'}, 1, TPR.nsizes), ' '), TPR.ndepths);
        TPR.tpr = single(cell2mat(D));

    % If the header is [SPdata]
    elseif strcmp(tagvalue{1},'[SPdata]')
        
        % Store the number of output factors (integer)
        D = textscan(fid, '%d', 1);
        SP.nsizes = D{1};
        
        % Store the field width values (vector of length nsizes)
        D = textscan(fid, '%f', SP.nsizes);
        SP.sizes = single(D{1});
        
        % Store the output factors (vector of length nsizes)
        D = textscan(fid, '%f', SP.nsizes);
        SP.values = single(D{1});

    % If the header is [OARXOPENdata]
    elseif strcmp(tagvalue{1},'[OARXOPENdata]')
        
        % Store the number of depths (integer)
        D = textscan(fid, '%d', 1);
        OARXOPEN.ndepths = D{1};
        
        % Store the OAR depths (vector of length ndepths)
        D = textscan(fid, '%f', OARXOPEN.ndepths);
        OARXOPEN.depths = single(D{1});
        
        % Store the number of values (integer)
        D = textscan(fid, '%d', 1);
        OARXOPEN.nvalues = D{1};
        
        % Store the OAR projection indices
        OARXOPEN.indices = single(0:50);

        % Store the OAR open field values (2D array of ndepths x nvalues)
        D = textscan(fid, strjoin(repmat({'%f'}, 1, ...
            OARXOPEN.ndepths), ' '), OARXOPEN.nvalues);
        OARXOPEN.oar = single(cell2mat(D));

    % If the header is [OARXLEAVESdata]
    elseif strcmp(tagvalue{1},'[OARXLEAVESdata]')
        
        % Store the number of widths (integer)
        D = textscan(fid, '%d', 1);    
        OARXLEAVES.nwidths = D{1};
        
        % Store the leaf widths (vector of length nwidths)
        D = textscan(fid, '%d', OARXLEAVES.nwidths);
        OARXLEAVES.widths = single(D{1});
        
        % Store the number of depths (integer)
        D = textscan(fid, '%d', 1);
        OARXLEAVES.ndepths = D{1};
        
        % Store the depths (vector of length ndepths)
        D = textscan(fid, '%f', OARXLEAVES.ndepths);
        OARXLEAVES.depths = single(D{1}); 
        
        % Store the number of leaf indices (integer)
        D = textscan(fid, '%d', 1);
        OARXLEAVES.nvalues = D{1};
        
        % Store the leaf indices (vector of length nvalues)
        D = textscan(fid, '%f', OARXLEAVES.nvalues);
        OARXLEAVES.indices = single(D{1});
        
        % Store the OAR leaf values (3D array of nvalues x ndepths x nwidths)
        D = textscan(fid, strjoin(repmat({'%f'}, 1, OARXLEAVES.ndepths), ...
            ' '), OARXLEAVES.nvalues * OARXLEAVES.nwidths);
        OARXLEAVES.oar=single(permute(reshape(cell2mat(D), ...
            OARXLEAVES.nvalues, OARXLEAVES.nwidths, OARXLEAVES.ndepths), ...
            [1 3 2])); 

    % If the header is [OARYdata]
    elseif strcmp(tagvalue{1}, '[OARYdata]')
        
        % Store the number of widths (integer)
        D = textscan(fid, '%d', 1);    
        OARY.nwidths = D{1};
        
        % Store the widths (vector of length nwidths)
        D = textscan(fid, '%f', OARY.nwidths);
        OARY.widths = single(D{1});

        % Store the number of depths (integer)
        D = textscan(fid, '%d', 1);
        OARY.ndepths = D{1};
        
        % Store the depths (vector of length ndepths)
        D = textscan(fid, '%f', OARY.ndepths);
        OARY.depths = single(D{1}); 
        
        % Store the number of indices (integer)
        D = textscan(fid, '%d', 1);
        OARY.nvalues = D{1};
        
        % Store the indices (vector of length nvalues)
        D = textscan(fid, '%f', OARY.nvalues);
        OARY.indices = single(D{1});
        
        % Store the OAR Y values (3D array of nvalues x ndepths x nwidths)
        D = textscan(fid, strjoin(repmat({'%f'}, 1, OARY.ndepths), ...
            ' '), OARY.nvalues * OARY.nwidths);
        OARY.oar=single(permute(reshape(cell2mat(D), OARY.nvalues, ...
            OARY.nwidths, OARY.ndepths), [1 3 2])); 
    end
end

% Clear temporary variables
clear i D fid tagvalue start_path files field_width;