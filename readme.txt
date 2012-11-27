This code implements an algorithm to automatically orient a 3D CT heart scan with contrast agent such that the principle axes of the left ventricle are parallel to the coordinate axes.
The algorithm requires a 3D CT advanced DICOM file as input which shows a cardiac phase where the left ventricle has a high contrast, while the right ventricle does not.
The code relies on the ITK framework. It has been tested with ITK 4.2.0. Before you use the code please make sure, the ITK libraries are available in a compiled form on your system. ITK is available at http://www.itk.org 
You can then use cmake to create the software project.

For information about the usage please read the section USAGE below.


The algorithm is structured as follows. Please refer to the source code for details.

1. 
Read the image data from file

2. 
Downsample the input image to speed up the following computations

3. 
Apply a threshold to eliminate low contrast structures and only keep high contrast structures

4. 
Perform a morphological opening to eliminate small structures

5. 
Conduct a connected component analysis and only keep the largest component i.e. the left ventricle and attached vessels. Features as e.g. size and principal axes for each component are calculated during this.

6.
Rotate the image about its centre to make the blob's principle axes parallel to the coordinate axes. The long ventricular axis will then be parallel to the z-axis

7.
Examine all the slices along the long ventricular axis and count the number of blobs in each slice. In the region of the ventricle there should only be one blob i.e. the ventricle. In the upper part there will be two blobs i.e. the left atrium and the aorta. The distance of their centers will get larger the farther up the examined slice is. This Information can be used to delineate the ventricle from the atrium and the aorta. Save the distances of the blobs in an array to obtain a distance distribution along the slices

8. 
Median filter the array of obtained blob distances to make the procedure more robust

8.1 
Calculate the regression line of the distances to identify the direction of the bifurcation of aorta and ventricle with reference to the maximum calculated distance

9.
Go from the slice with the largest blob distance along the decending regression line until the blob distance gets zero. Now keep only slices below that point (these contain the ventricle)

10.
Recalculate the principle axes of the remaining structure (hopefully more or less the left ventricle)

11.
Calculate the rotation of the image which is needed to make the structure's principle axes parallel to the coordinate axes and with it calculate the combined rotation of the rotation in step 5 and the now calculated rotation

12.
Apply the combined rotation to the original input image 

13.
Save the result


USAGE

The command line expects 6 parameters:

InputFilename:			Absolute filename of the input image

OutputFilename:			Absolute filename of the output image including the suffix

lowerThreshold: 		The lower threshold of the initial thresholding filter in Hounsfield Units

upperThreshold: 		The upperthreshold of the initial thresholding filter in Hounsfield Units

openingKernelSize:		The Kernelsize of the opening filter in voxel

selectOutput:			0 = Output transformed original image
				1 = Output transformed threshold image of left ventricle
				2 = Output transformed threshold image

Example: alignLeftVentricle.exe d:\imageData\inputImage.dcm d:\results\outputImage.dcm 200 1000 5 0