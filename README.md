## TomoTherapy MATLAB Dose Calculator

by Simon Thomas, adapted by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2017, University of Wisconsin Board of Regents
<br>Original work copyright &copy; 2011-15, Simon Thomas 

The TomoTherapy&reg; MATLAB Dose Calculator is a modified form of the original `dose_from_sin` function in the [CheckTomo tool](http://onlinelibrary.wiley.com/doi/10.1118/1.3668061/full) developed by Simon Thomas. It was adapted to calculate dose given input data provided in the format used by the [tomo_extract](https://github.com/mwgeurts/tomo_extract) and [dicom_tools](https://github.com/mwgeurts/dicom_tools) repositories.

TomoTherapy is a registered trademark of Accuray Incorporated.

## Installation

To install this function, copy all MATLAB .m and .txt files from this repository into your MATLAB path. If installing as a submodule into another git repository, execute `git submodule add https://github.com/mwgeurts/thomas_tomocalc`. To execute a dose calculation, call `CheckTomoDose`.

## Usage and Documentation

The function `CheckTomoDose` can be executed with various combinations of input arguments. Upon first execution, at least two arguments are required, as shown in the following example. See the [wiki documentation](https://github.com/mwgeurts/checktomo/wiki/Dose-Calculator-Runtime-Requirements) for more information on the input and return variable formats.

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

Finally, additional configuration options can be passed as name/value pairs for input arguments 4 and on. The available options are `downsample`, `reference_doserate`, `outside_body`, `density_threshold`, and `num_of_subprojections`.

```matlab
dose = CheckTomoDose(image, plan, pool, 'reference_doserate', 8.2);
dose = CheckTomoDose(image, plan, [], 'num_of_subprojections', 3);
```

For examples of how this function is used within larger applications, see the [checktomo](https://github.com/mwgeurts/checktomo) and [exit_detector](https://github.com/mwgeurts/exit_detector) repositories.

## License

Released under the GNU GPL v3.0 License.  See the [LICENSE](LICENSE) file for further details.
