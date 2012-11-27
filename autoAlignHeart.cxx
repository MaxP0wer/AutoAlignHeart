/*===========================================================================
 
   This is a program to automatically detect the left ventricular myocardial
   long axis in cardiac CT scans and align it parallel to the coordinate axes
	 The program requires ITK which is available at <http://www.itk.org>
	 It has been tested with ITK version 4.2.0.
 
   Copyright (C) 2012 Till Huelnhagen, Marc Dewey, Henning Meyer
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 ===========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkTimeProbe.h"
#include "itkGDCMImageIO.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryShapeKeepNObjectsImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

//******************************
// Global type definitions
//******************************
// Image Types
const unsigned int InputDimension = 3;

typedef itk::Image< signed short, 2 > SignedShortImageType2D;
typedef itk::Image< signed short, 3 > SignedShortImageType3D;
typedef itk::Image< unsigned char, 2 > UnsignedCharImageType2D;
typedef itk::Image< unsigned char, 3 > UnsignedCharImageType3D;
typedef itk::Image< itk::Vector< unsigned char , 3>, 2> ColorImageType2D;
typedef itk::Image< itk::Vector< unsigned char , 3>, 3> ColorImageType3D;
	
// Further type definitions
typedef itk::LinearInterpolateImageFunction< SignedShortImageType3D, double > InterpolatorType;
typedef itk::LinearInterpolateImageFunction< UnsignedCharImageType3D, double > InterpolatorBWType;
typedef itk::AffineTransform< double, InputDimension > AffineTransformType;

typedef itk::ResampleImageFilter< SignedShortImageType3D, SignedShortImageType3D > ResampleFilterType;
typedef itk::ResampleImageFilter< UnsignedCharImageType3D, UnsignedCharImageType3D > ResampleFilterBWType;

typedef itk::BinaryBallStructuringElement< unsigned char, 3 > StructuringElementType;
typedef itk::BinaryErodeImageFilter< UnsignedCharImageType3D, UnsignedCharImageType3D, StructuringElementType > ErodeFilterType;
typedef itk::BinaryDilateImageFilter< UnsignedCharImageType3D,	UnsignedCharImageType3D, StructuringElementType > DilateFilterType;

typedef itk::LabelMap< itk::ShapeLabelObject< itk::SizeValueType, 2> > LabelMapType2D;
typedef itk::LabelMap< itk::ShapeLabelObject< itk::SizeValueType, 3> > LabelMapType3D;

// function prototypes
SignedShortImageType3D::Pointer resampleImage( SignedShortImageType3D::Pointer image, AffineTransformType* transform, InterpolatorType* interpolator, SignedShortImageType3D::SpacingType& outputSpacing, bool enlargeImage = false);
UnsignedCharImageType3D::Pointer resampleImage( UnsignedCharImageType3D::Pointer image, AffineTransformType* transform, InterpolatorBWType* interpolator, UnsignedCharImageType3D::SpacingType& outputSpacing, bool enlargeImage = false);

UnsignedCharImageType3D::Pointer binaryThreshold( SignedShortImageType3D::Pointer image, SignedShortImageType3D::PixelType lowerThreshold, SignedShortImageType3D::PixelType upperThreshold, unsigned char insideValue, unsigned char outsideValue); 
UnsignedCharImageType2D::Pointer binaryThreshold2D( SignedShortImageType2D::Pointer image, SignedShortImageType2D::PixelType lowerThreshold, SignedShortImageType2D::PixelType upperThreshold, unsigned char insideValue, unsigned char outsideValue);

LabelMapType2D::Pointer createLabelMap2D( UnsignedCharImageType3D::Pointer BWImage2D );
LabelMapType3D::Pointer createLabelMap( UnsignedCharImageType3D::Pointer BWImage );

SignedShortImageType3D::Pointer createImageFromLabelMap( LabelMapType3D::Pointer labelMap );
UnsignedCharImageType3D::Pointer createBWImageFromLabelMap2D( LabelMapType2D::Pointer labelMap );

double getDistanceOfLargestBlobs ( SignedShortImageType3D::Pointer image, int sliceNumber );
LabelMapType2D::Pointer getBlobs( SignedShortImageType3D::Pointer image, int sliceNumber );

double median( std::vector<double> vec );
std::vector<double> medianFilterVector( std::vector<double> vec, int radius);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//******************************
// Main program
//******************************
int main( int argc, char **argv )
{
	// Validate input parameters
	if( argc < 7 )
	{
		std::cerr << "\nInvalid Number of Arguments. \n\n Usage: " 
			<< argv[0]
		<< " InputFilename OutputFilename lowerThreshold upperThreshold openingKernelSize selectOutput"
			<< std::endl;
		return EXIT_FAILURE;
	}

	itk::TimeProbe clock;
	clock.Start();

	//***********************************  
	// 1. Read the input image
	//***********************************
	typedef itk::ImageFileReader< SignedShortImageType3D > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

	typedef itk::GDCMImageIO           ImageIOType;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  
	reader->SetImageIO( gdcmImageIO );
	reader->SetFileName( argv[1] );
  
	std::cout << "\nReading image... ";
	try
	{
		reader->Update();
		clock.Stop();
		std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
		clock.Start();
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while reading the input image" << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	//**************************************
	// 2. Downsample the image to save time
	//**************************************
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
 
  const SignedShortImageType3D::SpacingType& inputSpacing = reader->GetOutput()->GetSpacing();
  const SignedShortImageType3D::RegionType& inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
  const SignedShortImageType3D::SizeType& inputSize = inputRegion.GetSize();
 
	std::cout << "Input image Information: " << std::endl;
	std::cout << "Origin:  " << reader->GetOutput()->GetOrigin() << std::endl;
	std::cout << "Size: " << inputSize << std::endl;
	std::cout << "Spacing: " << inputSpacing << std::endl;

	// Change image spacing to 1.5mm for downsampling
  SignedShortImageType3D::SpacingType outputSpacing;
	outputSpacing[0] = 1.5 * sgn(inputSpacing[0]);
	outputSpacing[1] = 1.5 * sgn(inputSpacing[1]);
	outputSpacing[2] = 1.5 * sgn(inputSpacing[2]);

 
	std::cout << "\nDownsampling image... ";

	AffineTransformType::Pointer identityTransform = AffineTransformType::New();
  identityTransform->SetIdentity();

	SignedShortImageType3D::Pointer image = resampleImage(reader->GetOutput(), identityTransform, interpolator, outputSpacing, true);

	clock.Stop();
	std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	std::cout << "New size: " << image->GetLargestPossibleRegion().GetSize() << "\n" << std::endl;

	// calculate center of the image
	SignedShortImageType3D::PointType imageCenter;
	SignedShortImageType3D::IndexType imageCenterIndex;
	SignedShortImageType3D::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();
	imageCenterIndex[0] = imageSize[0] / 2 + 0.5;
	imageCenterIndex[1] = imageSize[1] / 2 + 0.5;
	imageCenterIndex[2] = imageSize[2] / 2 + 0.5;
	image->TransformIndexToPhysicalPoint( imageCenterIndex, imageCenter );

	//***************************************
	// 3. Thresholding
	//***************************************
	std::cout << "Thresholding... ";

	UnsignedCharImageType3D::PixelType outsideValue, insideValue;
	SignedShortImageType3D::PixelType lowerThreshold, upperThreshold;
	outsideValue = 0;
	insideValue = 255;
	lowerThreshold = ::atoi(argv[3]);
	upperThreshold = ::atoi(argv[4]);

	UnsignedCharImageType3D::Pointer thresholdImage = binaryThreshold(image, lowerThreshold, upperThreshold, insideValue, outsideValue);

	clock.Stop();
	std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	//***************************************
	// 4. Strong morphological opening
	//***************************************
	std::cout << "Morphological opening... ";

	ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();
	DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

	StructuringElementType structuringElement;
	structuringElement.SetRadius( atoi(argv[5]) ); // 3x3 structuring element for radius 1
	structuringElement.CreateStructuringElement();

	binaryErode->SetKernel( structuringElement );
	binaryDilate->SetKernel( structuringElement );

	binaryErode->SetInput( thresholdImage );
	binaryDilate->SetInput( binaryErode->GetOutput());

	binaryErode->Update();

	clock.Stop();
	std::cout << "\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	//***************************************
	// 5. Convert binary image to label map
	//***************************************
	std::cout << "Creating label map... ";

	LabelMapType3D::Pointer labelMap = createLabelMap( binaryDilate->GetOutput());

	clock.Stop();
	std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();
	std::cout << "There are " << labelMap->GetNumberOfLabelObjects() << " objects." << std::endl;

	if( labelMap->GetNumberOfLabelObjects() == 0 ) {
		std::cout << "ERROR: No objects could be found in the image. Please check thresholds and morphology kernel size. \n\nAbort" << std::endl;
		return EXIT_FAILURE;
	}

	//***************************************
	// 5.1 Remove all but the largest object
	//***************************************
	std::cout << "\nRemoving all but the largest struture... ";

	//find largest object
	int objectNumber = labelMap->GetNumberOfLabelObjects();
	double maxSize = 0;
	int maxLabel = 0;
	for( int j = 0; j < objectNumber; j++) {
		LabelMapType3D::LabelObjectType* currentLabelObject = labelMap->GetNthLabelObject(j);
		double currentSize = currentLabelObject->GetNumberOfPixels();
		if( currentSize > maxSize) {
			maxSize = currentSize;
			maxLabel = j;
		}
	}
	//std::cout << "Label with maximum size: " << maxLabel << std::endl;
	//delete the remaining objects
	while( objectNumber > 0) {
		objectNumber--;
		if( objectNumber != maxLabel) {
			labelMap->RemoveLabelObject(labelMap->GetNthLabelObject(objectNumber));
		}
	}
	maxLabel = 0;

	clock.Stop();
	std::cout << "\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();
	
	// print out info
	std::cout << "There is(are) " << labelMap->GetNumberOfLabelObjects() << " object(s) remaining." << std::endl;
	LabelMapType3D::LabelObjectType* labelObject = labelMap->GetNthLabelObject(maxLabel);

	//************************************
	// 5.2 Create Label Image from Label Map
	//************************************
	std::cout << "Creating label map... ";

	SignedShortImageType3D::Pointer labelImage = createImageFromLabelMap( labelMap );

	clock.Stop();
	std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	//******************************************************************************************************
	// 6. Rotate the image to make the principal axes of the main component parallel to the coordinate axes
	//******************************************************************************************************
	const SignedShortImageType3D::SizeType& inputImageSize = image->GetBufferedRegion().GetSize();

	AffineTransformType::Pointer trafo = AffineTransformType::New();
	trafo->SetIdentity();

	trafo->SetCenter( imageCenter );
	trafo->SetMatrix(labelObject->GetPrincipalAxesToPhysicalAxesTransform()->GetMatrix());

	std::cout << "\nResampling Image... ";

	// Do the resampling
	ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
		resampler2->SetTransform(trafo);
		resampler2->SetInput( labelImage );  // Set the Label Image as output to be written
		resampler2->SetInterpolator( interpolator );
		resampler2->SetDefaultPixelValue( 0 );
		resampler2->SetOutputOrigin( labelImage->GetOrigin());
		resampler2->SetOutputSpacing( outputSpacing );
		resampler2->SetOutputDirection( labelImage->GetDirection());
		resampler2->SetSize( inputImageSize );
		resampler2->Update();

	clock.Stop();
	std::cout << "\t\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	//************************************
	// 7. Analyse blobs in all slices
	//************************************
	std::cout << "Analyzing blobs in each slice... ";

	// create stringstream to save data for statisical analysis
	std::ostringstream stringStream;
	stringStream << "slice; distance; roundness0; roundness1; size0; size1; ellipsoidDiameter1_0; ellipsoidDiameter1_1; ellipsoidDiameter2_0; ellipsoidDiameter2_1; elongation0; elongation1" << std::endl; 

	double distance = 0;
	std::vector<double> distances;
	std::vector<int> numberOfBlobs;
	int maxDistanceSlice = 0;
	double maxDistance = 0;
	UnsignedCharImageType2D::PointType c0_old, c1_old, c0, c1;
	c0_old.Fill( 0 );
	c1_old.Fill( 0 );

	for( int sliceNumber = 0; sliceNumber < inputImageSize[2]; sliceNumber++ ) {
		distance = 0;

		// conduct blob analysis for current slice
		LabelMapType2D::Pointer labelMap = getBlobs( resampler2->GetOutput(), sliceNumber);
		int currentNumberOfBlobs = labelMap->GetNumberOfLabelObjects();
		numberOfBlobs.push_back( currentNumberOfBlobs );

		if( currentNumberOfBlobs == 1 ) {
			LabelMapType2D::LabelObjectType::Pointer object0 = labelMap->GetNthLabelObject( 0 );
			stringStream << sliceNumber << "; 0; " << object0->GetRoundness() << "; 0; " 
				<< object0->GetNumberOfPixels() << "; 0; " 
				<< object0->GetEquivalentEllipsoidDiameter()[0] << "; 0; " 
				<< object0->GetEquivalentEllipsoidDiameter()[1] << "; 0; "
				<< object0->GetElongation() << "; 0"
				<< std::endl;
			c0_old = object0->GetCentroid();
			c1_old.Fill ( 0 );
		}

		else if( currentNumberOfBlobs > 1 ) {
			LabelMapType2D::LabelObjectType::Pointer object0 = labelMap->GetNthLabelObject( 0 );
			LabelMapType2D::LabelObjectType::Pointer object1 = labelMap->GetNthLabelObject( 1 );

			// check assignment of labels
			UnsignedCharImageType2D::PointType c0new = object0->GetCentroid();
			UnsignedCharImageType2D::PointType c1new = object1->GetCentroid();
			double d0 = sqrt( (c0new[0] - c0_old[0]) * (c0new[0] - c0_old[0]) + (c0new[1] - c0_old[1]) * (c0new[1] - c0_old[1]) );
			double d1 = sqrt( (c1new[0] - c0_old[0]) * (c1new[0] - c0_old[0]) + (c1new[1] - c0_old[1]) * (c1new[1] - c0_old[1]) );
			if( d0 > d1 ) {
				object1 = labelMap->GetNthLabelObject( 0 );
				object0 = labelMap->GetNthLabelObject( 1 );
			}

			c0 = object0->GetCentroid();
			c1 = object1->GetCentroid();
			c0_old = c0;
			c1_old = c1;
			double sizeRatio = double(object0->GetNumberOfPixels()) / double(object1->GetNumberOfPixels());
			sizeRatio = (sizeRatio>1)? 1/sizeRatio:sizeRatio;
			if( sizeRatio > 0.1 ) {
				distance = sqrt( (c0[0] - c1[0]) * (c0[0] - c1[0]) + (c0[1] - c1[1]) * (c0[1] - c1[1]) );
				if (distance > maxDistance) {
					maxDistance = distance;
					maxDistanceSlice = sliceNumber;
				}
			}
			stringStream << sliceNumber << "; " << distance << "; " 
				<< object0->GetRoundness() << "; " << object1->GetRoundness() << "; " 
				<< object0->GetNumberOfPixels() << "; " << object1->GetNumberOfPixels() << "; " 
				<< object0->GetEquivalentEllipsoidDiameter()[0] << "; " << object1->GetEquivalentEllipsoidDiameter()[0] << "; "
				<< object0->GetEquivalentEllipsoidDiameter()[1] << "; " << object1->GetEquivalentEllipsoidDiameter()[1] << "; "
				<< object0->GetElongation() << "; " << object1->GetElongation() << std::endl;
		}

		else
			stringStream << sliceNumber << "; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0" << std::endl;
		distances.push_back( distance );
	}

	clock.Stop();
	std::cout << "\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	clock.Start();

	// write resulting data to file
	std::ofstream myfile;
  myfile.open ("data.txt");
	myfile << stringStream.str();
  myfile.close();
	std::cout << "\nOutput written to data.txt and distancesMedian.txt" << std::endl;

	//****************************************************
	// 8. Calculate regression line of the blob distances
	//****************************************************
	// calculate mean in x and y
	double sumY = 0, meanX, meanY;
	for(int x = 0; x < distances.size(); x++) {
		sumY += distances[x];
	}
	meanX = distances.size() / 2;
	meanY = sumY / distances.size();

	double numerator = 0, denominator = 0;

	for(int x = 0; x < distances.size(); x++) {
		numerator += ( x - meanX ) * ( distances[x] - meanY );
		denominator += ( x - meanX ) * ( x - meanX );
	}
	double slope = numerator / denominator;
	std::cout << "Slope of regression line: " << slope << std::endl;

	std::vector<double> filteredVector = medianFilterVector( distances, 2 );
  myfile.open ("distancesMedian.txt");
	for( int i = 0; i < filteredVector.size(); i++)
		myfile << filteredVector[i] << std::endl;
  myfile.close();
	std::cout << "\nOutput written to distancesMedian.txt" << std::endl;

	//************************************************
	// 9. Find slice to cut to delinate the ventricle
	//************************************************
	int cutSlice = maxDistanceSlice;

	while( filteredVector[cutSlice] > 0) {
		cutSlice = cutSlice - sgn(slope) * 1;
	}
	std::cout << "\nCut slice: " << cutSlice << "\n" << std::endl;

	//extract image region with ventricle
	typedef itk::ExtractImageFilter< SignedShortImageType3D, SignedShortImageType3D > ExtractImageFilterType;

	UnsignedCharImageType3D::IndexType desiredStart;
	desiredStart.Fill( 0 );
	UnsignedCharImageType3D::SizeType desiredSize;

	// check which side of the image to keep
	if (slope > 0) {
		desiredSize = resampler2->GetOutput()->GetLargestPossibleRegion().GetSize();
		desiredSize[2] = cutSlice;
	}
	else {
		desiredStart[2] = cutSlice;
		desiredSize = resampler2->GetOutput()->GetLargestPossibleRegion().GetSize();
		desiredSize[2] = desiredSize[2] - cutSlice;
	}

	SignedShortImageType3D::RegionType desiredRegion(desiredStart, desiredSize);

	ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
		extractImageFilter->SetInput( resampler2->GetOutput() );
		extractImageFilter->SetExtractionRegion( desiredRegion );
		extractImageFilter->SetDirectionCollapseToIdentity();
		extractImageFilter->Update();

	//****************************************************************
	//	10. Recalculate the principle axes of the remaining structure
	//****************************************************************
	AffineTransformType::Pointer trafo2 = AffineTransformType::New();
	trafo2->SetIdentity();
	// thresholding
	UnsignedCharImageType3D::Pointer croppedImage = binaryThreshold( extractImageFilter->GetOutput(), 1, 255, 255, 0);
	// create label map
	LabelMapType3D::Pointer labelMap2 = createLabelMap( croppedImage );

	// set Transform center to center of the image
	trafo2->SetCenter( imageCenter );
	
	//*****************************************************************
	//	11. Get the rotation and calculate the combined roation matrix
	//*****************************************************************
	itk::Matrix<double,3,3> rotationMatrix;
	rotationMatrix = trafo->GetMatrix() * labelMap2->GetNthLabelObject( 0 )->GetPrincipalAxesToPhysicalAxesTransform()->GetMatrix();
	//std::cout << "Rotation 1 :" << trafo->GetMatrix() << std::endl;
	//std::cout << "Rotation 2 :" << labelMap2->GetNthLabelObject( 0 )->GetPrincipalAxesToPhysicalAxesTransform()->GetMatrix() << std::endl;
	//std::cout << "Combined Rotation :" << rotationMatrix << std::endl;

	trafo2->SetMatrix( rotationMatrix );
	
	itk::Vector< double, 3 > axis;
	axis.Fill( 0 );
	axis[1] = 1;
	trafo2->Rotate3D( axis, 1.57, true );


	//********************************************************************
	//	12. Transform the input image using the calculated transformation
	//********************************************************************
	SignedShortImageType3D::Pointer resultingImage = resampleImage( image, trafo2, interpolator, outputSpacing);

  //**********************
	// 13. Write output image
	//**********************
	std::cout << "\nWriting output image... ";
	typedef itk::ImageFileWriter< SignedShortImageType3D > WriterType;
	WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( argv[2] );
		writer->SetImageIO( gdcmImageIO );
	if( std::atoi( argv[6] ) == 1 )
		writer->SetInput( extractImageFilter->GetOutput() );
	else if( std::atoi( argv[6] ) == 2 )
		writer->SetInput( resampler2->GetOutput() );
	else
		writer->SetInput( resultingImage );
  try
  {
	  writer->Update();
		clock.Stop();
		std::cout << "\t\t\tdone\t( " << clock.GetMean() << "s )\n" << std::endl;
	  std::cout << "Transformation successful" << std::endl;
		std::cout << "Total time: " << clock.GetTotal() << "s" << std::endl;
  }
  catch( itk::ExceptionObject & err )
  {
	  std::cerr << "ExceptionObject caught !" << std::endl;
	  std::cerr << err << std::endl;
	  return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


UnsignedCharImageType3D::Pointer resampleImage( UnsignedCharImageType3D::Pointer image, AffineTransformType* transform, InterpolatorBWType* interpolator, UnsignedCharImageType3D::SpacingType& outputSpacing, bool enlargeImage )
{
	//calculate the output size according to the desired output spacing
	const UnsignedCharImageType3D::SizeType& inputSize = image->GetLargestPossibleRegion().GetSize();
	const UnsignedCharImageType3D::SpacingType& inputSpacing = image->GetSpacing();

	UnsignedCharImageType3D::SizeType   outputSize;
  typedef UnsignedCharImageType3D::SizeType::SizeValueType SizeValueType;
	if( !enlargeImage ) {
		outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
		outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
		outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);
	}
	else {
		outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
		outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
		outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);
		SizeValueType diagonal = static_cast<SizeValueType>( sqrt( double(outputSize[0]*outputSize[0] + outputSize[1]*outputSize[1] + outputSize[2]*outputSize[2]) ) );
		
		AffineTransformType::OutputVectorType translation;
		translation[0] = -outputSpacing[0] * ( diagonal - outputSize[0] ) / 2;
		translation[1] = -outputSpacing[1] * ( diagonal - outputSize[1] ) / 2;
		translation[2] = -outputSpacing[2] * ( diagonal - outputSize[2] ) / 2;
		transform->Translate( translation );

		outputSize[0] = static_cast<SizeValueType>( diagonal );
		outputSize[1] = static_cast<SizeValueType>( diagonal );
		outputSize[2] = static_cast<SizeValueType>( diagonal );
	}

	ResampleFilterBWType::Pointer resampler = ResampleFilterBWType::New();
	resampler->SetInput(image);
	resampler->SetTransform(transform);
	resampler->SetInterpolator(interpolator);
	resampler->SetOutputOrigin ( image->GetOrigin());
	resampler->SetOutputSpacing ( outputSpacing );
	resampler->SetOutputDirection ( image->GetDirection());
	resampler->SetSize ( outputSize );
	resampler->SetDefaultPixelValue( 0 );
	resampler->Update ();

	return resampler->GetOutput();
}

// Image resampling function
SignedShortImageType3D::Pointer resampleImage( SignedShortImageType3D::Pointer image, AffineTransformType* transform, InterpolatorType* interpolator, SignedShortImageType3D::SpacingType& outputSpacing, bool enlargeImage ) 
{
	//calculate the output size according to the desired output spacing
	const SignedShortImageType3D::SizeType& inputSize = image->GetLargestPossibleRegion().GetSize();
	const SignedShortImageType3D::SpacingType& inputSpacing = image->GetSpacing();

	SignedShortImageType3D::SizeType outputSize;
  typedef SignedShortImageType3D::SizeType::SizeValueType SizeValueType;
	if( !enlargeImage ) {
		outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
		outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
		outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);
	}
	else {
		outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
		outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
		outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);
		SizeValueType diagonal = static_cast<SizeValueType>( sqrt( double(outputSize[0]*outputSize[0] + outputSize[1]*outputSize[1] + outputSize[2]*outputSize[2]) ) );
		
		AffineTransformType::OutputVectorType translation;
		translation[0] = -outputSpacing[0] * ( diagonal - outputSize[0] ) / 2;
		translation[1] = -outputSpacing[1] * ( diagonal - outputSize[1] ) / 2;
		translation[2] = -outputSpacing[2] * ( diagonal - outputSize[2] ) / 2;
		transform->Translate( translation );

		outputSize[0] = static_cast<SizeValueType>( diagonal );
		outputSize[1] = static_cast<SizeValueType>( diagonal );
		outputSize[2] = static_cast<SizeValueType>( diagonal );
	}

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetInput(image);
	resampler->SetTransform(transform);
	resampler->SetInterpolator(interpolator);
	resampler->SetOutputOrigin ( image->GetOrigin());
	resampler->SetOutputSpacing ( outputSpacing );
	resampler->SetOutputDirection ( image->GetDirection());
	resampler->SetSize ( outputSize );
	resampler->SetDefaultPixelValue( 100 );
	resampler->Update ();

	return resampler->GetOutput();
}

// Thresholding function
UnsignedCharImageType3D::Pointer binaryThreshold( SignedShortImageType3D::Pointer image, SignedShortImageType3D::PixelType lowerThreshold, SignedShortImageType3D::PixelType upperThreshold, unsigned char insideValue, unsigned char outsideValue)
{
	typedef itk::BinaryThresholdImageFilter<SignedShortImageType3D, UnsignedCharImageType3D> BinaryThresholdFilterType;

	BinaryThresholdFilterType::Pointer binaryThresholdFilter = BinaryThresholdFilterType::New();
	binaryThresholdFilter->SetInput(image);
	binaryThresholdFilter->SetOutsideValue( outsideValue );
	binaryThresholdFilter->SetInsideValue( insideValue );
	binaryThresholdFilter->SetUpperThreshold(upperThreshold);
	binaryThresholdFilter->SetLowerThreshold(lowerThreshold);
	binaryThresholdFilter->Update();

	return binaryThresholdFilter->GetOutput();
}

UnsignedCharImageType2D::Pointer binaryThreshold2D( SignedShortImageType2D::Pointer image, SignedShortImageType2D::PixelType lowerThreshold, SignedShortImageType2D::PixelType upperThreshold, unsigned char insideValue, unsigned char outsideValue)
{
	typedef itk::BinaryThresholdImageFilter<SignedShortImageType2D, UnsignedCharImageType2D> BinaryThresholdFilterType;

	BinaryThresholdFilterType::Pointer binaryThresholdFilter = BinaryThresholdFilterType::New();
	binaryThresholdFilter->SetInput(image);
	binaryThresholdFilter->SetOutsideValue( outsideValue );
	binaryThresholdFilter->SetInsideValue( insideValue );
	binaryThresholdFilter->SetUpperThreshold(upperThreshold);
	binaryThresholdFilter->SetLowerThreshold(lowerThreshold);
	binaryThresholdFilter->Update();

	return binaryThresholdFilter->GetOutput();
}

// Function to create a label map from a binary image
LabelMapType3D::Pointer createLabelMap( UnsignedCharImageType3D::Pointer BWImage )
{
	typedef itk::BinaryImageToShapeLabelMapFilter< UnsignedCharImageType3D > BinaryImageToShapeLabelMapFilterType;
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilter->SetInput(BWImage);
	binaryImageToShapeLabelMapFilter->Update();

	return binaryImageToShapeLabelMapFilter->GetOutput();
}

// Function to create a label map from a binary image
LabelMapType2D::Pointer createLabelMap2D( UnsignedCharImageType2D::Pointer BWImage )
{
	typedef itk::BinaryImageToShapeLabelMapFilter< UnsignedCharImageType2D > BinaryImageToShapeLabelMapFilterType;
	BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilter->SetInput(BWImage);
	binaryImageToShapeLabelMapFilter->Update();

	return binaryImageToShapeLabelMapFilter->GetOutput();
}

SignedShortImageType3D::Pointer createImageFromLabelMap( LabelMapType3D::Pointer labelMap )
{
	typedef itk::LabelMapToLabelImageFilter<LabelMapType3D, SignedShortImageType3D> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter->SetInput( labelMap );
	labelMapToLabelImageFilter->Update();

	return labelMapToLabelImageFilter->GetOutput();
}

UnsignedCharImageType3D::Pointer createBWImageFromLabelMap( LabelMapType3D::Pointer labelMap )
{
	typedef itk::LabelMapToLabelImageFilter<LabelMapType3D, UnsignedCharImageType3D> LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter->SetInput( labelMap );
	labelMapToLabelImageFilter->Update();

	return labelMapToLabelImageFilter->GetOutput();
}

double getDistanceOfLargestBlobs ( SignedShortImageType3D::Pointer image, int sliceNumber )
{
	typedef itk::ExtractImageFilter< SignedShortImageType3D, SignedShortImageType2D > ExtractImageFilterType;
	UnsignedCharImageType3D::IndexType desiredStart;
	desiredStart[0] = 0;
	desiredStart[1] = 0;
	desiredStart[2] = sliceNumber;

	UnsignedCharImageType3D::SizeType desiredSize;
	desiredSize = image->GetLargestPossibleRegion().GetSize();
	desiredSize[2] = 0;

	SignedShortImageType3D::RegionType desiredRegion(desiredStart, desiredSize);

	ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
	
	extractImageFilter->SetInput( image );
	extractImageFilter->SetExtractionRegion( desiredRegion );
	extractImageFilter->SetDirectionCollapseToIdentity();
	extractImageFilter->Update();

	// Thresholding
	UnsignedCharImageType2D::Pointer bwImage = binaryThreshold2D( extractImageFilter->GetOutput(), 1, 255, 255, 0);

	//keep only the largest two blobs
	typedef itk::BinaryShapeKeepNObjectsImageFilter<UnsignedCharImageType2D> BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryShapeKeepNObjectsImageFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryShapeKeepNObjectsImageFilter->SetInput(bwImage);
	binaryShapeKeepNObjectsImageFilter->SetNumberOfObjects( 2 );
	binaryShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
	binaryShapeKeepNObjectsImageFilter->SetReverseOrdering( false );
	binaryShapeKeepNObjectsImageFilter->SetAttribute( "NumberOfPixels" );
	binaryShapeKeepNObjectsImageFilter->Update();

	// Create label map
	LabelMapType2D::Pointer labelMap = createLabelMap2D( binaryShapeKeepNObjectsImageFilter->GetOutput() );

	std::cout << "Slice number: " << sliceNumber << std::endl;
	double distance = 0;
	if( labelMap->GetNumberOfLabelObjects() > 1) {
		// check size ratio of the blobs > 10%
		double sizeRatio = double(labelMap->GetNthLabelObject(0)->GetNumberOfPixels()) / double(labelMap->GetNthLabelObject(1)->GetNumberOfPixels());
		sizeRatio = (sizeRatio>1)? 1/sizeRatio:sizeRatio;
		std::cout << "Size ratio: " << sizeRatio << std::endl;
		if( sizeRatio > 0.1 ) {
			UnsignedCharImageType2D::PointType c0 = labelMap->GetNthLabelObject(0)->GetCentroid();
			UnsignedCharImageType2D::PointType c1 = labelMap->GetNthLabelObject(1)->GetCentroid();
			distance = sqrt( (c0[0] - c1[0]) * (c0[0] - c1[0]) + (c0[1] - c1[1]) * (c0[1] - c1[1]) );
			std::cout << "Distance between centroids: " << distance << std::endl;
			std::cout << "Size0: " << labelMap->GetNthLabelObject(0)->GetNumberOfPixels() << "\nSize1: " << labelMap->GetNthLabelObject(1)->GetNumberOfPixels() << std::endl;
		}
	}

	return distance;
}

LabelMapType2D::Pointer getBlobs( SignedShortImageType3D::Pointer image, int sliceNumber ) 
{
	typedef itk::ExtractImageFilter< SignedShortImageType3D, SignedShortImageType2D > ExtractImageFilterType;
	UnsignedCharImageType3D::IndexType desiredStart;
	desiredStart[0] = 0;
	desiredStart[1] = 0;
	desiredStart[2] = sliceNumber;

	UnsignedCharImageType3D::SizeType desiredSize;
	desiredSize = image->GetLargestPossibleRegion().GetSize();
	desiredSize[2] = 0;

	SignedShortImageType3D::RegionType desiredRegion(desiredStart, desiredSize);

	ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
	
	extractImageFilter->SetInput( image );
	extractImageFilter->SetExtractionRegion( desiredRegion );
	extractImageFilter->SetDirectionCollapseToIdentity();
	extractImageFilter->Update();

	// Thresholding
	UnsignedCharImageType2D::Pointer bwImage = binaryThreshold2D( extractImageFilter->GetOutput(), 1, 255, 255, 0);

	//keep only the largest two blobs
	typedef itk::BinaryShapeKeepNObjectsImageFilter<UnsignedCharImageType2D> BinaryShapeKeepNObjectsImageFilterType;
	BinaryShapeKeepNObjectsImageFilterType::Pointer binaryShapeKeepNObjectsImageFilter = BinaryShapeKeepNObjectsImageFilterType::New();
	binaryShapeKeepNObjectsImageFilter->SetInput(bwImage);
	binaryShapeKeepNObjectsImageFilter->SetNumberOfObjects( 2 );
	binaryShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
	binaryShapeKeepNObjectsImageFilter->SetReverseOrdering( false );
	binaryShapeKeepNObjectsImageFilter->SetAttribute( "NumberOfPixels" );
	binaryShapeKeepNObjectsImageFilter->Update();

	// Create label map
	LabelMapType2D::Pointer labelMap = createLabelMap2D( binaryShapeKeepNObjectsImageFilter->GetOutput() );
	
	return labelMap;
}

// returns the median of a vector
double median( std::vector<double> vec )
{
	typedef std::vector<double>::size_type vec_sz;

	vec_sz size = vec.size();
	
	std::nth_element(vec.begin(), vec.begin()+size, vec.end());

	vec_sz mid = size/2;

	if( size % 2 == 0) {
		std::nth_element(vec.begin(), vec.begin() + mid - 1, vec.end());
		double mid1 = vec[mid - 1];
		std::nth_element(vec.begin(), vec.begin() + mid, vec.end());
		double mid2 = vec[mid]; 
		return ( mid1 + mid2 ) / 2;
	}
	else {
		std::nth_element(vec.begin(), vec.begin() + mid, vec.end());
		return vec[ mid ];
	}
}

std::vector<double> medianFilterVector( std::vector<double> vec, int radius)
{
	std::vector<double> filteredVector( vec.size() );
	// copy entries at end and beginning which are not filtered
	for( int i = 0; i < radius; i++)
		filteredVector[i] = vec[i];
	for( int i = vec.size() - 1; i > vec.size() - radius; i--)
		filteredVector[i] = vec[i];
	// fill middle values
	for( int i = radius; i < vec.size() - radius; i++) {
		std::vector<double> tempVector( 2*radius + 1 );
		for( int j = 0; j < 2*radius + 1; j++) {
			tempVector[j] = vec[ j - radius + i];
			filteredVector[i] = median( tempVector );
		}
	}
	return filteredVector;
}