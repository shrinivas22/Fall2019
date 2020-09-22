/************************************************************
 *															*
 * This sample project include three functions:				*
 * 1. Add intensity for gray-level image.					*
 *    Input: source image, output image name, value			*
 *															*
 * 2. Image thresholding: pixels will become black if the	*
 *    intensity is below the threshold, and white if above	*
 *    or equal the threhold.								*
 *    Input: source image, output image name, threshold		*
 *															*
 * 3. Image scaling: reduction/expansion of 2 for 			*
 *    the width and length. This project uses averaging 	*
 *    technique for reduction and pixel replication			*
 *    technique for expansion.								*
 *    Input: source image, output image name, scale factor	*
 *															*
 ************************************************************/

#include "./iptools/core.h"
#include <opencv2/opencv.hpp>
#include <string.h>
#include <fstream>


using namespace std;
using  namespace cv;

#define MAXLEN 256
/* Histogram equalization with mask only for polygonal foreground regions*/
void equalizeHistWithMask(const Mat1b& src, Mat1b& dst, Mat1b mask = Mat1b())
{
	int cnz = countNonZero(mask);
	if (mask.empty() || (cnz == src.rows * src.cols))
	{
		equalizeHist(src, dst);
		return;
	}

	dst = src.clone();

	// Histogram
	vector<int> hist(256, 0);
	for (int r = 0; r < src.rows; ++r) {
		for (int c = 0; c < src.cols; ++c) {
			if (mask(r, c)) {
				hist[src(r, c)]++;
			}
		}
	}

	// Cumulative histogram
	float scale = 255.f / float(cnz);
	vector<uchar> lut(256);
	int sum = 0;
	for (int i = 0; i < hist.size(); ++i) {
		sum += hist[i];
		lut[i] = saturate_cast<uchar>(sum * scale);
	}

	// Apply equalization
	for (int r = 0; r < src.rows; ++r) {
		for (int c = 0; c < src.cols; ++c) {
			if (mask(r, c)) {
				dst(r, c) = lut[src(r, c)];
			}
		}
	}
}

/**/
void shiftDFT(Mat& fImage)
{
	Mat tmp, q0, q1, q2, q3;

	// first crop the image, if it has an odd number of rows or columns

	fImage = fImage(Rect(0, 0, fImage.cols & -2, fImage.rows & -2));

	int cx = fImage.cols / 2;
	int cy = fImage.rows / 2;

	// rearrange the quadrants of Fourier image
	// so that the origin is at the image center

	q0 = fImage(Rect(0, 0, cx, cy));
	q1 = fImage(Rect(cx, 0, cx, cy));
	q2 = fImage(Rect(0, cy, cx, cy));
	q3 = fImage(Rect(cx, cy, cx, cy));

	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);
}

/******************************************************************************/
// return a floating point spectrum magnitude image scaled for user viewing
// complexImg- input dft (2 channel floating point, Real + Imaginary fourier image)
// rearrange - perform rearrangement of DFT quadrants if true

// return value - pointer to output spectrum magnitude image scaled for user viewing

Mat magnitude_display(Mat & complexImg)
{
	Mat planes[2];


	split(complexImg, planes);
	magnitude(planes[0], planes[1], planes[0]);

	Mat mag = (planes[0]).clone();
	mag += Scalar::all(1);
	log(mag, mag);
		// re-arrange the quaderants
		shiftDFT(mag);
	normalize(mag, mag, 0, 1, NORM_MINMAX);

	return mag;

}
/******************************************************************************/

// create various filters based on inputs of filter, radius D, order n


void create_filter(Mat& dft_Filter, int D,string filter,int D2=0, int D3=0)
{
	Mat tmp = Mat(dft_Filter.rows, dft_Filter.cols, CV_32F);

	Point centre = Point(dft_Filter.rows / 2, dft_Filter.cols / 2);
	double radius;

	// based on the forumla in the IP notes (p. 130 of 2009/10 version)
	// see also HIPR2 on-line
	for (int i = 0; i < dft_Filter.rows; i++)
	{
		for (int j = 0; j < dft_Filter.cols; j++)
		{
			if (filter == "lowpass"  || filter=="lowPassAndNotch")//|| filter == "highpass"
			{
				radius = (double)sqrt(pow((i - centre.x), 2.0) + pow((double)(j - centre.y), 2.0));


				if (radius > D)
				{
					tmp.at<float>(i, j) = 0;
				}
				//tmp.at<float>(i, j) = (float)
				//	(1 / (1 + pow((double)(radius / D), (double)(2 * n))));
			}
		}
	}
	
		for (int i = 0; i < dft_Filter.rows; i++)
		{
			for (int j = 0; j < dft_Filter.cols; j++)
			{
				radius = (double)sqrt(pow((i - centre.x), 2.0) + pow((double)(j - centre.y), 2.0));
				
			
				 if (filter == "highpass" || filter == "highPassAndBandPass")
				 {
					 if (filter == "highPassAndBandPass")
					 {
						 if (radius < D && radius>D3)
						 {
							 tmp.at<float>(i, j) = 0;
						 }
					 }
					 else
					 {
						 if (radius < D)
						 {
							 tmp.at<float>(i, j) = 0;//1 - tmp.at<float>(i, j)
						 }
					 }
				 }
				 if (filter == "bandpass" || filter == "highPassAndBandPass")
				 {
					 if(filter == "highPassAndBandPass")
					 {
						 if (radius < D2)
						 {
							 tmp.at<float>(i, j) = 0;
						 }
						 
					 }
					 else {
						 if (radius < D2)
						 {
							 tmp.at<float>(i, j) = 0;
						 }
						 else if (radius > D3) {
							 tmp.at<float>(i, j) = 0;
						 }
					 }
				 }
				 if (filter == "bandstop"|| filter == "lowPassAndNotch")
				{

					if (radius > D2 && radius < D3)
					{
						tmp.at<float>(i, j) = 0;
					}
				}



			}
		}
	
		
		//cv::imshow("highpass_filter", dft_Filter);

	Mat toMerge[] = { tmp, tmp };
	merge(toMerge, 2, dft_Filter);
}

// A butterworth high pass filter for image sharpening  

void butterworth(Mat& img, int order, int Rad)
{
	Mat tmp = Mat(img.rows, img.cols, CV_32F);


	Point centre = Point(img.rows / 2, img.cols / 2);
	double radius;


	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			radius = (double)sqrt(pow((i - centre.x), 2.0) + pow((double)(j - centre.y), 2.0));
			tmp.at<float>(i, j) = (float)
				(1 / (1 + pow((double)(radius / Rad), (double)(2 * order))));
		}
	}
	//Experimneting with highpass
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			radius = (double)sqrt(pow((i - centre.x), 2.0) + pow((double)(j - centre.y), 2.0));

			if (radius < Rad)
			{
				tmp.at<float>(i, j) = 1 - tmp.at<float>(i, j);
			}
		}
	}

	Mat toMerge[] = { tmp, tmp };
	merge(toMerge, 2, img);
}
/******************************************************************************/

void DFTOperations(Mat & image, Mat image_roi, Mat& previmg, Mat& opWithSpectrum,string filterName, int x1, int y1, int sx1, int sy1, int radiusThreshold , int radiusThreshold2 ,int radiusThreshold3 ) {
	//Mat src;
	//image.copyTo(src);

	//   y,x,y1,x1		
	
	Mat dst;
	//Mat opWithSpectrum = previmg.clone();
	printf("opWithSpectrum type: %d", opWithSpectrum.type());
	//waitKey(150);
	//DFT(image_roi);
	// convert input to grayscale
	Mat padded;		// fourier image objects and arrays
	Mat complexImg, filter, filter1, filterOutput, imgOutput, filterOutput1, imgOutput1;
	Mat planes[2], mag, mag1;

	int Rows, Columns; // fourier image sizes
	Rows = getOptimalDFTSize(image_roi.rows);
	Columns = getOptimalDFTSize(image_roi.cols);

	printf("OPtimal: %d,%d", Rows, Columns);
	int  order = 2;				// butterworth filter parameter

	// setup the DFT images
	copyMakeBorder(image_roi, padded, 0, Rows - image_roi.rows, 0,
		Columns - image_roi.cols, BORDER_CONSTANT, Scalar::all(0));
	planes[0] = Mat_<float>(padded);
	planes[1] = Mat::zeros(padded.size(), CV_32F);
	printf("OPtimal padded: %d,%d", padded.rows, padded.cols);
	merge(planes, 2, complexImg);

	// do the DFT

	dft(complexImg, complexImg);
	printf("TypeMag:%d \n", mag.type());
	// construct the filter (same size as complex image)

	filter = complexImg.clone();
	filter1 = complexImg.clone();

	//create_butterworth_lowpass_filter(filter, radius, order,true); Tried for extra credits

	 butterworth(filter1, order,30);

	create_filter(filter, radiusThreshold, filterName, radiusThreshold2, radiusThreshold3);
	// apply filter
	//Rearranging the quadrants  in fourier domain
	shiftDFT(complexImg);
	// for convolution of dft image and filtered image output
	Mat outputForButterWorth = complexImg.clone();
	mulSpectrums(complexImg, filter, complexImg, 0);
	
	mulSpectrums(outputForButterWorth, filter1, outputForButterWorth, 0);
	shiftDFT(complexImg);
	shiftDFT(outputForButterWorth);
	// create magnitude spectrum for display

	mag = magnitude_display(complexImg);
	mag1 = magnitude_display(outputForButterWorth);
	printf("Mag type: %d",mag.type());
	imshow("originalNameROI", image_roi);
	mag.convertTo(mag, CV_8U, 255.0);
	imshow("originalName", image);

	int x11=0, y11 = 0;
	for (int i = x1; i < sx1; i++) {
		y11 = 0;
		for (int j = y1; j < sy1; j++) {

			//printf("\n %d,%d",x,y);
			opWithSpectrum.at<uchar>(i, j) = mag.at<uchar>(x11, y11);
			y11++;
		}
		x11++;
	}
	printf("x/y:%d,%d,%d,%d", x11, y11,mag.rows,mag.cols);
	imshow(filterName+"SpectrumMagOnImage", opWithSpectrum);
	//opWithSpectrum.convertTo(opWithSpectrum, CV_8U, 255.0);
	imshow("spectrumMagName", mag);
	imshow("spectrumMagNamebUTTERwORTH", mag1);

	// do inverse DFT on filtered image

	printf("TypeMag:%d,%d \n", mag.type(), filter.type());
	// apply inverse DFT
	idft(complexImg, complexImg);
	idft(outputForButterWorth, outputForButterWorth);

	// split into planes and extract plane 0 as output image

	split(complexImg, planes);
	normalize(planes[0], imgOutput, 0, 1, NORM_MINMAX);

	split(outputForButterWorth, planes);
	normalize(planes[0], imgOutput1, 0, 1, NORM_MINMAX);

	split(filter, planes);
	normalize(planes[0], filterOutput, 0, 1, NORM_MINMAX);
	printf("Type:%d,%d", image.type(), filterOutput.type());
	printf("filterop:%d,%d,%d,%c", filterOutput.rows, filterOutput.cols, filterOutput.at<int>(0, 0), image.at<uchar>(0, 0));
	filterOutput.convertTo(filterOutput, CV_8U, 255.0);
	imgOutput.convertTo(imgOutput, CV_8U, 255.0);

	int y, x = 0;
	for (int i = x1; i < sx1; i++) {
		y = 0;
		for (int j = y1; j < sy1; j++) {

			//printf("\n %d,%d",x,y);
			previmg.at<uchar>(i, j) = imgOutput.at<uchar>(x, y);
			y++;
		}
		x++;
	}
	/*Mat finalImg = src.clone();
	for (int i = 0; i < image.rows; i++) {
		for (int j = 0; j < image.cols; j++) {
			finalImg.at<uchar>(i, j) = image.at<uchar>(i, j);
		}
	}*/
	imshow("originalNameWithMag" + filterName, previmg);
	// ***



	imshow(filterName, imgOutput);
	imshow("highpass-bUTTERWORTH", imgOutput1);
	imshow(filterName+" Filter", filterOutput);
	//imwrite("../output/dft.pgm", filter);
	imwrite("../output/" + filterName+"-Mag"+ ".pgm", filterOutput);
	imwrite("../output/"+ filterName+".pgm", imgOutput);
	imwrite("../output/FilterInROI-"+ filterName + ".pgm", previmg);
	imwrite("../output/MagSpectrum-" + filterName + ".pgm", mag);
	imwrite("../output/"+filterName + "SpectrumMagOnImage" +".pgm", opWithSpectrum);

	//imwrite("../output/BWLP.pgm", filterOutput);

	waitKey(0);
}


void DFT(Mat &src) {// get optimal size of DFT
	int  optimalRows = getOptimalDFTSize(src.rows);
	int optimalCols = getOptimalDFTSize(src.cols);
	Mat padded;
	copyMakeBorder(src, padded, 0, optimalRows - src.rows, 0,
		optimalCols - src.cols, BORDER_CONSTANT, Scalar::all(0));

	// use cv.MatVector to distribute space for real part and imaginary part
	Mat plane0;
	padded.convertTo(plane0, CV_32F);
	Mat complexI;
	Mat planes[] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
	merge(planes, 2, complexI);

	// in-place dft transform
	dft(complexI, complexI);
	//imshow("DFT", complexI);
	// compute log(1 + sqrt(Re(DFT(img))**2 + Im(DFT(img))**2))
	split(complexI, planes);
	magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
	Mat mag = planes[0];

	mag += Scalar::all(1);                    // switch to logarithmic scale
	log(mag, mag);
	Mat m1 = Mat::ones(mag.rows, mag.cols, mag.type());

	// crop the spectrum, if it has an odd number of rows or columns
	mag = mag(Rect(0, 0, mag.cols & -2, mag.rows & -2));

	// rearrange the quadrants of Fourier image
	// so that the origin is at the image center
	int cx = mag.cols / 2;
	int cy = mag.rows / 2;
	Mat tmp;

	Mat q0(mag, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
	Mat q1(mag, Rect(cx, 0, cx, cy));  // Top-Right
	Mat q2(mag, Rect(0, cy, cx, cy));  // Bottom-Left
	Mat q3(mag, Rect(cx, cy, cx, cy)); // Bottom-Right

	// exchange 1 and 4 quadrants
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	// exchange 2 and 3 quadrants
	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);

	// The pixel value of cv.CV_32S type image ranges from 0 to 1.
	normalize(mag, mag, 0, 1, NORM_MINMAX);

	imshow("canvasOutput", src);
	imshow("MagnitudeDFTONLY", mag);
}

Mat HistogramPlotting(Mat image1,Mat out1, const int * hsize, const float * hRange) {
	calcHist(&image1, 1, { 0 }, Mat(), out1, 1, hsize, &hRange);
	printf("hist:%d,%d", out1.rows, out1.cols);
	double max_val = 0;
	minMaxLoc(out1, 0, &max_val);
	int bins = 256;
	int const hist_height = 256;

	cv::Mat hist_image = cv::Mat::zeros(hist_height, bins, CV_32F);
	// visualize each bin
	for (int b = 0; b < bins; b++) {
		float const binVal = out1.at<float>(b);
		int   const height = cvRound(binVal * (hist_height / max_val));
		cv::line
		(hist_image
			, cv::Point(b, hist_height - height), cv::Point(b, hist_height)
			, cv::Scalar::all(255)
		);
	}
	return hist_image;
}

/******************************************************************************/

int main (int argc, char** argv)
{
	image src, tgt, tgt2, tgt3, tgt4, tgt5, tgt6, tgt7, bg, fg, bgCombination, fgCombination, bgHist, fgHist, combined;
	ifstream fp(argv[1]);
	char str[MAXLEN], str2[MAXLEN], a[MAXLEN], b[MAXLEN], a1[MAXLEN], b1[MAXLEN], a2[MAXLEN], b2[MAXLEN];
	String inputPath;
	rsize_t strmax = sizeof str;
	char outfile[MAXLEN],outfile1[MAXLEN], outfile2[MAXLEN], outfile3[MAXLEN], outfile4[MAXLEN],outfile5[MAXLEN];
	char *pch, *next_pch;
	int nOP;
	int arr[50];
	if (!fp.is_open()) {
		fprintf(stderr, "Can't open file: %s\n", argv[1]);
		exit(1);
	}

	fp >> nOP;

	for (int OP = 0; OP < nOP; OP++) {
		fp >> str;
	 inputPath= str;
		src.read(str);
		
		fp >> str;
		strcpy_s(outfile, MAXLEN, str);

		/*if (OP == 5) {
			fp >> str;
			strcpy_s(outfile1, MAXLEN, str);
		}*/

		fp >> str;
        if (strncmp(str,"add",3)==0) {
			/* Add Intensity */
			fp >> str;
        	utility::addGrey(src,tgt,atoi(str));
        }

        else if (strncmp(str,"binarize",8)==0) {
			/* Thresholding */
			fp >> str;
			utility::binarize(src,tgt,atoi(str));
		}

		else if (strncmp(str,"scale",5)==0) {
			/* Image scaling */
			fp >> str;
			utility::scale(src,tgt,atof(str));
		}

		/*else if (strncmp(str, "Histograms", 10) == 0) {
			
		
			utility::drawHistogram(src, tgt);
		}*/

		else if (strncmp(str, "HistogramStretching", 19) == 0) {
		
			for (int i = 0; i < 15; i++)
			{
				fp >> str;
				arr[i] = atoi(str);
			}
			fp >> a;
			fp >> b;
			fp >> a1;
			fp >> b1;
			fp >> a2;
			fp >> b2;

			utility::drawHistogram(src, tgt5, tgt6, tgt7, arr);
			strcpy_s(outfile2, MAXLEN, "../output/lena_Before_ROI1.pgm");
			strcpy_s(outfile3, MAXLEN, "../output/lena_Before_ROI2.pgm");
			strcpy_s(outfile4, MAXLEN, "../output/lena_Before_ROI3.pgm");
			tgt5.save(outfile2);
			tgt6.save(outfile3);
			tgt7.save(outfile4);
			utility::histogramStretching(src, tgt, tgt2, tgt3, tgt4, atoi(a), atoi(b), atoi(a1), atoi(b1), atoi(a2), atoi(b2), arr, false);
			//strcpy_s(outfile1, MAXLEN, "../output/StretchedThresholdROI1.pgm");
			strcpy_s(outfile2, MAXLEN, "../output/lena_stretched.pgm");
			strcpy_s(outfile3, MAXLEN, "../output/StretchedThresholdROI2.pgm");
			strcpy_s(outfile4, MAXLEN, "../output/StretchedThresholdROI3.pgm");
			tgt2.save(outfile2);
			tgt3.save(outfile3);
			tgt4.save(outfile4);

		}
		else if (strncmp(str, "OptimalThreshold", 16) == 0) {
			
			for (int i = 0; i < 15; i++)
			{
				fp >> str;
				arr[i] = atoi(str);
			}
			utility::optimalthresholding(src, tgt,bg,fg, arr,false);
			printf("%s", str);
		} 
		else if (strncmp(str, "Combined", 8) == 0) {
			/* Image scaling */
			for (int i = 0; i < 4; i++)
			{
				fp >> str;
				arr[i] = atoi(str);
			}
			utility::Combination(src, tgt, bg, fg, bgCombination, fgCombination, bgHist, fgHist,combined, arr);
			strcpy_s(outfile1, MAXLEN, "../output/lena_fg_AfterOT.pgm");
			strcpy_s(outfile2, MAXLEN, "../output/lena_bg_AfterOT.pgm");
			strcpy_s(outfile3, MAXLEN, "../output/lena_optimal_HistStrectched_fg.pgm");
			strcpy_s(outfile4, MAXLEN, "../output/lena_optimal_HistStrectched_bg.pgm");
			strcpy_s(outfile5, MAXLEN, "../output/lena_combined_fg_bg.pgm");
			bg.save(outfile2);
			fg.save(outfile1);
			fgHist.save(outfile4);
			bgHist.save(outfile3);
			combined.save(outfile5);
		}
		else if (strncmp(str, "SobelD", 6) == 0)
		{
			for (int i = 0; i < 16; i++)
			{
				fp >> str;
				arr[i] = atoi(str);
			}
			fp >> str;
			utility::SobelEdgeDetector(src, tgt,tgt2,tgt3,atoi(str),arr);
			strcpy_s(outfile2, MAXLEN, "../output/sobelAmplitude.pgm");
			strcpy_s(outfile1, MAXLEN, "../output/sobelAgle.pgm");
			tgt2.save(outfile1);
			tgt3.save(outfile2);
		}

	/*http://breckon.eu/toby/teaching/dip/opencv/lecture_demos/
	 Took help from various codes in the above repository for DFT, IDFTand filtering*/
			
		else if (strncmp(str, "DFTFunctions", 12) == 0) {
			/*-------------------------------------------------------------------*/
			char x1[MAXLEN], y1[MAXLEN], sx1[MAXLEN], sy1[MAXLEN], x2[MAXLEN], y2[MAXLEN], sx2[MAXLEN], sy2[MAXLEN], x3[MAXLEN], y3[MAXLEN], sx3[MAXLEN], sy3[MAXLEN], RThres[MAXLEN], RThres1[MAXLEN], RThres2[MAXLEN];
			fp >> x1;
			fp >> y1;
			fp >> sx1;
			fp >> sy1;
			fp >> x2;
			fp >> y2;
			fp >> sx2;
			fp >> sy2;
			fp >> x3;
			fp >> y3;
			fp >> sx3;
			fp >> sy3;
			fp >> RThres; fp >> RThres1; fp >> RThres2;

			
			Mat image = imread(inputPath, IMREAD_GRAYSCALE);\
				DFT(image);
			string DFTfilters[] = { "lowpass","highpass","bandpass","bandstop","lowPassAndNotch","highPassAndBandPass" };
			
		
			for (int i=0; i<6;i++)
			{
				Mat image = imread(inputPath, IMREAD_GRAYSCALE);

				Mat img = image.clone();
				printf("img::%d", img.type());
				Mat img1 = image.clone();
				Mat imgMag = image.clone();
				cv::Rect roi(atoi(x1), atoi(y1), atoi(sx1), atoi(sy1));
				//Create the cv::Mat with the ROI you need, where "image" is the cv::Mat you want to extract the ROI from
				cv::Mat image_roi = Mat(atoi(sx1) - atoi(x1), atoi(sy1) - atoi(y1), CV_8U);
				printf("image_roi::%d", image_roi.type());
				//image_roi=image(roi)  ;
				string filter = DFTfilters[i];
				cv::Rect Region(atoi(x2), atoi(y2), atoi(sx2), atoi(sy2));
				cv::Mat image_roi1 = Mat(atoi(sx2) - atoi(x2), atoi(sy2) - atoi(y2), CV_8U);
				//image_roi1 = img(Region);
				cv::Rect R1(atoi(x3), atoi(y3), atoi(sx3), atoi(sy3));
				cv::Mat image_roi2 = Mat(atoi(sx3) - atoi(x3), atoi(sy3) - atoi(y3), CV_8U);
					//image_roi2 = img1(R1);


					int x = 0, y = 0,x11 = 0, y11 = 0, x22 = 0, y22 = 0;
					for (int i = 0; i < img.rows; i++) {
						y = 0;
						y11 = 0; 
						y22 = 0;
						for (int j = 0; j < img.cols; j++) {

							if (i >= atoi(x1) && j >= atoi(y1) && i < atoi(sx1) && j < atoi(sy1))
							{
								image_roi.at<uchar>(x, y) = img.at<uchar>(i, j);
								y++;
							}
							if (i >= atoi(x2) && j >= atoi(y2) && i < atoi(sx2) && j < atoi(sy2))
							{
								image_roi1.at<uchar>(x11, y11) = img.at<uchar>(i, j);
								y11++;
							}
							if (i >= atoi(x3) && j >= atoi(y3) && i < atoi(sx3) && j < atoi(sy3))
							{
								image_roi2.at<uchar>(x22, y22) = img.at<uchar>(i, j);
								y22++;
							}
						}
						if (i >= atoi(x1) && i < atoi(sx1))
						{
							x++;
						}if (i >= atoi(x2) &&  i < atoi(sx2) )
						{
							x11++;
						}if (i >= atoi(x3)  && i < atoi(sx3))
							{
							x22++;
							}
					}
				// Done All filterings for 3 ROIS
				DFTOperations(image, image_roi,image , imgMag, filter,atoi(x1), atoi(y1), atoi(sx1), atoi(sy1), atoi(RThres), atoi(RThres1), atoi(RThres2));
				DFTOperations(img, image_roi1, image, imgMag,filter, atoi(x2), atoi(y2), atoi(sx2), atoi(sy2), atoi(RThres), atoi(RThres1), atoi(RThres2));
				DFTOperations(img1, image_roi2, image, imgMag, filter, atoi(x3), atoi(y3), atoi(sx3), atoi(sy3), atoi(RThres), atoi(RThres1), atoi(RThres2));
			}

		}
		else if (strncmp(str, "OPENCVFunctions", 15) == 0) {

			Mat image = imread(inputPath, IMREAD_GRAYSCALE);
			Mat image1 = imread(inputPath, IMREAD_GRAYSCALE);

			Mat out, out1;
		int bins = 256;
		int hsize1[] = { bins };
		float xranges[] = { 0, 256 };
		const float* hRange = xranges;
		const int* hsize = hsize1;
		Mat outHist, gray_image, otsu, outHist1, outHist2;
		equalizeHist(image, outHist);
		int const hist_height = 256;
		//cv::Mat hist_image = cv::Mat::zeros(hist_height, bins, CV_32F);
		/*--------------------- Calculate Histogram and plotting----------------------*/
		Mat hist_image=HistogramPlotting(image1, out1,hsize, hRange);
		imshow("histPlot", hist_image);
		namedWindow("Hist");
		imshow("Hist", outHist);
		/*----------------------------------------------------------------------------*/
		namedWindow("Original");
		imshow("Original", image);
		Sobel(image, out, -10, 1, 1); /* -10 describes the ddepth ofthe image*/
		namedWindow("Sobel1");
		imshow("Sobel1", out);
		imwrite("../output/Sobel.pgm", out);
		Sobel(image, out, -5, 1, 0); 
		namedWindow("SobelY");
		imshow("SobelY", out);
		imwrite("../output/SobelY.pgm", out);
		Sobel(image, out, -5, 0, 1);
		namedWindow("SobelX");
		imshow("SobelX", out);
		imwrite("../output/SobelX.pgm", out);
		Sobel(image, out, -10, 1, 1,5);
		namedWindow("Sobel5x5");
		imshow("Sobel5x5", out);
		imwrite("../output/Sobel5x5.pgm", out);
		Sobel(image, out, -5, 1, 0,5);
		namedWindow("Sobel5x5Y");
		imshow("Sobel5x5Y", out);
		imwrite("../output/Sobel5x5Y.pgm", out);
		Sobel(image, out, -5, 0, 1,5); 
		namedWindow("Sobel5x5X");
		imshow("Sobel5x5X", out);
		imwrite("../output/Sobel5x5X.pgm", out);


		Canny(image, out, 100, 200);
		namedWindow("canny");
		imshow("canny", out);
		imwrite("../output/canny.pgm", out);

		waitKey(0);
		imwrite("../output/Sobel.pgm", out);
		imwrite("../output/HistPlot.pgm", outHist);
}
		else if (strncmp(str, "OPENCVHistEqualization_OTSU", 27) == 0) {
		int bins = 256;
		int hsize1[] = { bins };
		float xranges[] = { 0, 256 };
		const float* hRange = xranges;
		const int* hsize = hsize1;
		Mat image = imread(inputPath, IMREAD_GRAYSCALE), otsu,adaptThres, outHist2, outHist,out1;
		/* Equalized Histogram image*/
		equalizeHist(image, outHist);
		imshow("EqualizationHist", outHist);
		Mat hist_image = HistogramPlotting(outHist, out1, hsize, hRange);
		imshow("EqualizationHistPlot", hist_image);

		imwrite("../output/HistogramEq.pgm", hist_image);

		/* OTSU */
		double Threshold = threshold(image, otsu, 0, 255, THRESH_OTSU);
		adaptiveThreshold(image, adaptThres, 255, ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,11,2);
		//threshold(image, otsu, 0, 255, THRESH_OTSU);
		printf(" Thresh,%lf", Threshold);
		imshow("OTSU", otsu);
		imshow("adaptThres", adaptThres);


		vector<Point> pts;
		for (int i = 0; i < otsu.rows; i++) {
			for (int j = 0; j < otsu.cols; j++) {
				if (otsu.at<uchar>(i, j) == 255)
					pts.push_back(Point(i, j));
			}
		}

		Mat mask(otsu.rows, otsu.cols, uchar(0));
		vector<vector<Point>> ptsarray{ pts };
		fillPoly(mask, ptsarray, Scalar(255));

		Mat1b equalized;
		equalizeHistWithMask(image, equalized, mask);
	
		/*-------------------------------------------------------------------*/

		/*-----------------Histogram Equalization for OTSU's Foreground---------------------*/
		Mat otsu1 = otsu.clone();
		for (int i = 0; i < otsu.rows; i++) {
			for (int j = 0; j < otsu.cols; j++) {
				if (otsu.at<uchar>(i, j) > Threshold)
					otsu1.at<uchar>(i, j) = 255;
				else

					otsu1.at<uchar>(i, j) = image.at<uchar>(i, j);
			}
		}

		equalizeHist(otsu1, outHist2);
		printf(" equalizeHist :%d", outHist2.type());
		Mat hist_image1 = HistogramPlotting(outHist2, out1, hsize, hRange);
		imshow("histPlot1", hist_image1);

		imread("../input/13Stretched.pgm", IMREAD_GRAYSCALE);

		Mat hist_image12 = HistogramPlotting(imread("../input/lena_stretched.pgm", IMREAD_GRAYSCALE), out1, hsize, hRange);
		imshow("lena_stretched", hist_image12);
		imwrite("../output/Stretched.pgm", hist_image12);
		imwrite("../output/HistogramForeground.pgm", hist_image1);

		/*	for (int i = 0; i < outHist2.rows; i++) {
			for (int j = 0; j < outHist2.cols; j++) {
				if (outHist2.at<uchar>(i, j) != 255)
					outHist2.at<uchar>(i, j) = image.at<uchar>(i, j);
			}
		}*/

		imshow("HistEqNew", outHist2);
		waitKey(0);
		imwrite("../output/OTSU.pgm", otsu);
		imwrite("../output/HistEqualized.pgm", outHist);
		imwrite("../output/ForegroundEqualized.pgm", outHist2);
}
		else {
			printf("No function: %s\n", str);
			continue;
		}
       
		tgt.save(outfile);
		
	}
	//fclose(fp);
	fp.close();
	return 0;
}

