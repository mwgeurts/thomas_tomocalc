## TomoTherapy MATLAB Dose Calculator

by Simon Thomas, adapted by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2017, University of Wisconsin Board of Regents
<br>Original work copyright &copy; 2011-15, Simon Thomas 

The TomoTherapy&reg; MATLAB Dose Calculator is a modified form of the original `dose_from_sin` function in the [CheckTomo tool](http://onlinelibrary.wiley.com/doi/10.1118/1.3668061/full) developed by Simon Thomas. It was adapted to calculate dose given input data provided in the format used by the [tomo_extract](https://github.com/mwgeurts/tomo_extract) and [dicom_tools](https://github.com/mwgeurts/dicom_tools) repositories. The function was also expanded to allow computation on a MATLAB cluster to improve computation time. Various other minor changes have been made to different lines of the code; refer to the comments in this and the other functions of this repository for details.

TomoTherapy is a registered trademark of Accuray Incorporated.

## Contents

* [Installation and Use](README.md#installation-and-use)
* [Compatibility and Requirements](README.md#compatibility-and-requirements)
* [Methods](README.md#methods)
* [Troubleshooting](README.md#troubleshooting)
* [License](README.md#license) 

## Installation and Use

To install this function, copy all MATLAB .m files from this repository into your MATLAB path. If installing as a submodule into another git repository, execute `git submodule add https://github.com/mwgeurts/thomas_tomocalc`. To execute a dose calculation, call `CheckTomoDose`.

The function `CheckTomoDose` can be executed with various combinations of input arguments. Upon first execution, at least two arguments are required, as shown in the following example, where image and plan are structures returned by the [tomo_extract](https://github.com/mwgeurts/tomo_extract) functions `LoadImage` and `LoadPlan`:

```matlab
dose = CheckTomoDose(image, plan);
```

The image is stored persistently, so after the first call, a second plan may be calculated with only one input argument:

```matlab
dose = CheckTomoDose(plan);
```

A parallel pool can also be passed in the third input argument (the first two must be image and plan):

```matlab
dose = CheckTomoDose(image, plan, pool);
```

Finally, additional configuration options can be passed as name/value pairs for input arguments 4 and on. The available options are 'downsample', 'reference_doserate', and 'num_of_subprojections'. See the code below for detail on each option:

```matlab
dose = CheckTomoDose(image, plan, pool, 'reference_doserate', 8.2);
dose = CheckTomoDose(image, plan, [], 'num_of_subprojections', 3);
```

Upon successful completion, this function returns the structure dose, which contains the following fields:

* start: 1 x 3 vector of position cooordinates, in cm, for the lower leftvoxel. These coordinates are identical to image.start.
* width: 1 x 3 vector of voxel dimensions, in cm
* dimensions: 1 x 3 vector of the size of dose.data
* data: 3D array of dose values, in Gy, of dimensions provided above
  
Below is complete example of how this function is used:

```matlab
% Load image and plan data for the plan specified by the UID below
image = LoadImage('./path/', 'TomoPhant^^^^_patient.xml', ...
    '1.2.826.0.1.3680043.2.200.1828118229.362.96568.1276');
plan = LoadPlan('./path/', 'TomoPhant^^^^_patient.xml', ...
    '1.2.826.0.1.3680043.2.200.1828118229.362.96568.1276');

% Execute CheckTomoDose with 2 workers and apply a downsampling factor 
pool = parpool(2);
dose = CheckTomoDose(image, plan, pool, 'downsample', 4);
```

## Compatibility and Requirements

These functions have been developed in MATLAB R2016b (9.1) on Mac OS X version 10.12.4 and Parallel Computing Toolbox version 6.9. The Parallel Computing Toolbox is only required if you wish to run the calculation against a local or MATLAB cluster parallel pool.

## Methods

For more information on the methods employed in this tool, see [Thomas et al. Independent dose calculation software for tomotherapy, Med Phys 2015; 39: 160-167](http://onlinelibrary.wiley.com/doi/10.1118/1.3668061/full). The original tool, CheckTomo, can be obtained through the GPL license by contacting Simon Thomas (refer to the correspondence address in the journal article referenced above for contact information).

## Troubleshooting

This function will report progress information by calling the function
% Event(message, flag) if it exists in the application path, where the
% message is a string and the flag is one of the following strings: 'INFO',
% 'WARN', or 'ERROR'. See the file at the following address for an example:
% https://github.com/mwgeurts/exit_detector/blob/master/Event.m

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
