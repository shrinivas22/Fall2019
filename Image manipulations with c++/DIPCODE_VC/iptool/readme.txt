*** This project is built under Visual Studio 2015


This software is architectured as follows

iptool/iptool.cpp- This file includes the main function.

iptool/iptools -This folder hosts the files that are compiled into a static library. 
	image - This folder hosts the files that define an image.
	utility- this folder hosts the files that students store their implemented algorithms.



*** INSTALATION ***

On Windows

Open the project by double click iptool.sln.

*** FUNCTIONS ***

1. Add intensity: add
Increase the intensity for a gray-level image.

2. Binarization: binarize
Binarize the pixels with the threshold.

3. Scaling: Scale
Reduce or expand the heigh and width with twp scale factors.
Scaling factor = 2: double height and width of the input image.
Scaling factor = 0.5: half height and width of the input image.



*** PARAMETERS ***

The first parameter of the parameters.txt is the number of operations (lines).
There are four parameters for each operation (line):
1. the input file name;
2. the output file name;
3. the name of the filter. Use "add", "binarize", and "scale" for your filters;

4.HistogramStretching:
Here, the histogram is stretched into given range so that the pixel intensites are evenly distributed.
The peaks are reduced and redistributed.
5.OptimalThresholding:
In this method,we find an optimum threshold for the image by looping several times until the
difference of avearage of the foreground and background is less than the user defined value 'x'.
6.Combination:
This process is the combination of the above 2 processes, that also involves highlighting the 
foreground and background images , stretching them separately and combining them.

*** PARAMETERS ***
HistogramStretching:
We have a total of 15 +6 =21 parameters for 3 ROIs
every 5th (in the first 15) parameter is dummy **Please Ignore**
every first 4 parameters until 15 are to determine the ROI as x,y,Sx,Sy
**(every Sx,Sy is taken as endpoint rather than length)**
the last 6 values correspond to the C,D values for 3 ROI(2*3=6)
where C,D are fromthe formula[ p'=p-C(B-A/D-C)+A]

Optimal Threshold:
We have a total of 15 +6 =21 parameters for 3 ROIs
every 5th (in the first 15) parameter 
is taken as the threshold value to be given ('x' a user given specific value)
used in the stopping condition
every first 4 parameters until 15 are to determine the ROI as x,y,Sx,Sy
**(every Sx,Sy is taken as endpoint rather than length)**

Combination:
The 4 values here correspond to the ROIS x,y,Sx,Sy.

SobelEdgeDetector:
Takes in  a total of 12+3+2=17 parameters for 3 ROIs
every 5th (in the first 15) parameter 
is taken as the threshold value to be given ('x' a user given specific value)
every first 4 parameters until 15 are to determine the ROI as x,y,Sx,Sy
then we have 2 values the first one being angle, Sobel kernel size.

Gives 3 outputs
1) SobelED gives the output after thresholding based on the gradient amplitude formula
2) SobelAngle applies direction on the detected edges based on the amplitude of the gradient
3) SobelAmplitude replaces pixel values with the normalised gradient value

OpenCV functions
it includes 
Sobel 5X5, 3X3 in X and Y directions
Canny edge detection
Then 
Histogram Stretching
OTSU algorithm
Taking the Foreground based on OTSU and doing histogram stretching on this 

***  DFTFunctions 
The first 4 values here correspond to the ROIS x,y,Sx,Sy. 
the next 2 sets of 4 values also correspond to ROIs. 
the next three values correspond to the radius thresholds R1,R2, R3 used for various filters.

*********Please NOte******************
http://breckon.eu/toby/teaching/dip/opencv/lecture_demos/
Took help from various codes in the above repository for DFT, IDFT and filtering
Open CV functions 
DFT(input Mat, output Mat)
gives the output of the spatial domain image in frequency domain
Need to perform normalization and rearranging of the planes to give a comprehendable input done using ShiftDFT function and magnitude_Display function
IDFT(input Mat, output Mat)

create_filter(input image, radius threshold1, filter name, radius threshold2, radius threshold3 )
Kindly note that RT1>RT3>RT2 in the above implemenation

Tried for extra credits- image sharpening 
butterworth(image, order, radius threshold)
applied high pass using the butterWorth filter formula


gives the output of the frequency domain image in spatial  domain
***All the image names are self-explanatory kindly let me know if any issues***

*** Run the program: For this 

Directly debug in Visual Studio.
You can find the output image in output folder.


*** Important information ***

Application assumes the next format of input image (ppm/pgm) file:
line1: <version>
line2: <#columns> <#rows>
line3: <max_value>
line4-end-of-file:<pix1><pix2>...<pix_n>

if it is a grayscale image then every pixel is a char value. If it is a RGB image then every pixel is a set of 3 char values: <red><green><blue>

Thus, if you have a problem with reading image, the first thing you should check is input file format.
