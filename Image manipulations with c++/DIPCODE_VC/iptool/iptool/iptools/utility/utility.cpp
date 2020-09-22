#include "utility.h"
#include <map>
# include <algorithm>
# include <math.h>
#include <vector>
//#include <bits/stdc++.h> 
using namespace std;
#define MAXRGB 255
#define MINRGB 0
#define MAXLEN 70000
#define PI  3.14159265358979323846 
static int maxvalue = 0;
static int minvalue = 0;

std::string utility::intToString(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

int utility::checkValue(int value)
{
	if (value > MAXRGB)
		return MAXRGB;
	if (value < MINRGB)
		return MINRGB;
	return value;
}

/*-----------------------------------------------------------------------**/
void utility::addGrey(image &src, image &tgt, int value)
{
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++)
		for (int j=0; j<src.getNumberOfColumns(); j++)
		{
			tgt.setPixel(i,j,checkValue(src.getPixel(i,j)+value)); 
		}
}

/*-----------------------------------------------------------------------**/
void utility::binarize(image &src, image &tgt, int threshold)
{
	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	for (int i=0; i<src.getNumberOfRows(); i++)
	{
		for (int j=0; j<src.getNumberOfColumns(); j++)
		{
			if (src.getPixel(i,j) < threshold)
				tgt.setPixel(i,j,MINRGB);
			else
				tgt.setPixel(i,j,MAXRGB);
		}
	}
}

/*------*/

void utility :: histogramStretching(image& src, image & tgt, image& tgt2,image& tgt3,image& tgt4,int a, int b,int a1, int b1,int a2,int b2, int arr[],bool isCombination) {

	int x = arr[0], y = arr[1], x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		Sx = arr[2], Sy = arr[3], Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9];
	tgt2.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt3.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt4.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	/*for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			//ROI1 
			if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
			{
				if (src.getPixel(i, j) > maxval)
				{
					maxval = src.getPixel(i, j);
				}
				if (src.getPixel(i, j) < minval && minval!=0)
				{
					minval = src.getPixel(i, j);
				}
			}
			//ROI2
			else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
				(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3) && !isCombination) {
				if (src.getPixel(i, j) > maxval1)
				{
					maxval1 = src.getPixel(i, j);
				}
				if (src.getPixel(i, j) < minval1 && minval1 != 0)
				{
					minval1 = src.getPixel(i, j);
				}
			}
			//ROI3
			else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
				(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3) && !isCombination) {
				if (src.getPixel(i, j) > maxval2)
				{
					maxval2 = src.getPixel(i, j);
				}
				if (src.getPixel(i, j) < minval2 && minval2 != 0)
				{
					minval2 = src.getPixel(i, j);
				}
			}
		}
	}*/

	
	image hist1, hist2, hist3;

	for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			tgt2.setPixel(i, j, src.getPixel(i, j));
			int currentPixel = src.getPixel(i, j);
			int newPix = 0;
			//ROI1 
			if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
			{
				tgt2.setPixel(i, j, currentPixel);
				if (currentPixel != 0)
					newPix = ((currentPixel - a) * ((255 - 0) / (b - a))) + 0;
				tgt2.setPixel(i, j, newPix);
			}
			//ROI2
			else if (((i<x || i<y || i>Sx || j>Sy) &&
				(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2)&&!isCombination) {
				if (currentPixel != 0)
					newPix = ((currentPixel - a) * ((255 - 0) / (b1 - a1))) + 0;
				tgt2.setPixel(i, j, newPix);
			}
			//ROI3
			else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
				(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3) && !isCombination){
				if (currentPixel != 0)
					newPix = ((currentPixel - a) * ((255 - 0) / (b2 - a2))) + 0;
				tgt2.setPixel(i, j, newPix);
			}
		}
	}
	if(!isCombination)
	drawHistogram(tgt2, tgt,tgt3,tgt4,arr);
	//scale(tgt3,tgt, 0.5);
}

/*------*/
int utility::findMedian(int arr1[], int size) {
	int n = size;
	int median = 0;
	sort(arr1, arr1 + n);
	if (n % 2 == 0) {
		median = arr1[(n - 1) / 2];
	}
	else {
		median = arr1[((n - 1) / 2) + 1];
	}
	return median;
}
/*------*/
int utility:: findMaxMin(image& src,int arr[] ,bool isMax) {
	int x = arr[0], y = arr[1], x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		Sx = arr[2], Sy = arr[3], Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9];
	int maxval = 0;
	int minval = INT16_MAX;
	for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			//ROI1 
			if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
			{
				if (src.getPixel(i, j) > maxval)
				{
					maxval = src.getPixel(i, j);
				}
				if (src.getPixel(i, j) < minval && src.getPixel(i, j) != 0)
				{
					minval = src.getPixel(i, j);
				}
			}
		}
	}
	return isMax ? maxval : minval;
}

/*------*/


void utility::Combination(image & src, image& tgt, image& bg, image& fg, image& bgCombination, image& fgCombination, image& bgHist, image& fgHist,image & combined, int arr[])
{
	int x = arr[0], y = arr[1], Sx = arr[2], Sy = arr[3];
		/*x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		 Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9];*/

	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	fg.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	bg.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	combined.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	
	optimalthresholding(src, tgt, bg, fg,arr,true);
	int max = 0, min = 0;
	max = findMaxMin(fg, arr, true);
	min = findMaxMin(fg, arr, false);
	printf("fg max&min %d,%d:", max, min);
	// just to match the function call !
	image fgdummy,fgdummy2, bgdummy, bgdummy2;
	histogramStretching(fg,fgCombination, fgHist,fgdummy,fgdummy2,min,max,0,0,0,0,arr,true);
	max = findMaxMin(bg, arr, true);
	min = findMaxMin(bg, arr, false);
	histogramStretching(bg, bgCombination, bgHist, bgdummy, bgdummy2,min,max,0,0,0,0,arr,true);

	combined.copyImage(fgHist);

	for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
			{
				if (fgHist.getPixel(i, j) == 0)
					combined.setPixel(i, j, bgHist.getPixel(i, j));
			}
		}
	}
}





void utility::optimalthresholding(image& src, image& tgt, image& combinedBlackImg, image& combinedWhiteImg, int arr[],bool isCombination)
{
	int x = arr[0], y = arr[1], x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		Sx = arr[2], Sy = arr[3], Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9];
	int count = 0, count2 = 0, count3 = 0;
	int ROI1Mean = 0;
	int ROI2Mean = 0;
	int ROI3Mean = 0;
	int ROI1Median = 0;
	int ROI2Median = 0;
	int ROI3Median = 0;

	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	int arr11[MAXLEN], arr2[MAXLEN], arr3[MAXLEN];
		int k = 0,l = 0, m = 0;
	
		for (int i = 0; i < src.getNumberOfRows(); i++)
		{
			for (int j = 0; j < src.getNumberOfColumns(); j++)
			{
				tgt.setPixel(i, j, src.getPixel(i, j));
				if (isCombination)
				{
					combinedBlackImg.setPixel(i, j, src.getPixel(i, j));
					combinedWhiteImg.setPixel(i, j, src.getPixel(i, j));
				}

				if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
				{
					arr11[k++] = tgt.getPixel(i, j);
				}
				else if (((i<x || i<y || i>Sx || j>Sy) &&
					(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2))
				{
					arr2[l++] = tgt.getPixel(i, j);
				}
				else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
					(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3))
				{
					arr3[m++] = tgt.getPixel(i, j);
				}
			}
		}
		int n = k;
		//printf("array size %d,%d,%d:", sizeof(arr11),k,isCombination);
		int n1 = l;
		//printf("array size %d,%d,%d:", sizeof(arr2), l, isCombination);
		int n2 = m;
		//printf("array size %d,%d,%d:", sizeof(arr3), m, isCombination);
		int T1 = 0, T2 = 0, T3 = 0;
	
		
		int ROI1Fgr = 0, e = 0, ROI2Fgr = 0, f = 0, ROI3Fgr = 0, g = 0, ROI1Bgr = 0, h = 0, ROI2Bgr = 0, o = 0, ROI3Bgr = 0, p = 0;
		do {
			if (count == 0)
			{
				ROI1Median = findMedian(arr11,k);
				T1 = ROI1Median;
			}
			else {
				T1 = ROI1Median;
			}

			for (int i = 0; i < src.getNumberOfRows(); i++)
			{
				for (int j = 0; j < src.getNumberOfColumns(); j++)
				{
					if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
					{
						if (ROI1Median > src.getPixel(i, j)) {
							ROI1Bgr += src.getPixel(i, j);
							++e;
						}
						else {
							ROI1Fgr += src.getPixel(i, j);
							++f;
						}
					}
				}
			}
			ROI1Median = ((ROI1Bgr / e + ROI1Fgr / f) / 2);
			ROI1Mean = T1 - ROI1Median;
			count++;
		} while (Thresh < ROI1Mean);

		/*if (!isCombination) {
			do {
				if (count2 == 0)
				{
					ROI2Median = findMedian(arr2,l);
					T2 = ROI2Median;
				}
				else
					T2 = ROI2Median;

				for (int i = 0; i < src.getNumberOfRows(); i++)
				{
					for (int j = 0; j < src.getNumberOfColumns(); j++)
					{
						if (((i<x || i<y || i>Sx || j>Sy) &&
							(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2))
						{
							if (ROI2Median > src.getPixel(i, j)) {
								ROI2Bgr += src.getPixel(i, j);
								++g;
							}
							else {
								ROI2Fgr += src.getPixel(i, j);
								++h;

							}
						}
					}
				}
				ROI2Median = ((ROI2Bgr / g + ROI2Fgr / h) / 2);
				ROI2Mean = T2 - ROI2Median;
				count2++;
			} while (Thresh2 < ROI2Mean);

			do {
				if (count3 == 0)
				{
					ROI3Median = findMedian(arr3,m);
					T3 = ROI3Median;
				}
				else
					T3 = ROI3Median;

				for (int i = 0; i < src.getNumberOfRows(); i++)
				{
					for (int j = 0; j < src.getNumberOfColumns(); j++)
					{
						if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
							(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3))
						{
							if (ROI3Median > src.getPixel(i, j)) {
								ROI3Bgr += src.getPixel(i, j);
								++o;
							}
							else {
								ROI3Fgr += src.getPixel(i, j);
								++p;
							}
						}
					}
				}
				ROI3Median = ((ROI3Bgr / o + ROI3Fgr / p) / 2);
				ROI3Mean = T3 - ROI3Median;
				count3++;
			} while (Thresh3 < ROI3Mean);
		}*/
	
				//printf("%d,%d,%d \n", ROI1Mean, ROI2Mean, ROI3Mean);

		for (int i = 0; i < src.getNumberOfRows(); i++)
		{
			for (int j = 0; j < src.getNumberOfColumns(); j++)
			{
				if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
				{
					if (tgt.getPixel(i, j) > ROI1Median)
					{
						if (isCombination)
						{
							combinedBlackImg.setPixel(i, j, 0);
						}
						tgt.setPixel(i, j, 255);
						
					}
				else
				{
						if (isCombination)
						{
							combinedWhiteImg.setPixel(i, j, 0);
						}
					tgt.setPixel(i, j, 0);
					
				}
				}
				else if (((i<x || i<y || i>Sx || j>Sy) &&
					(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2) &&(!isCombination)) {
					
					if (tgt.getPixel(i, j) > ROI2Median)
						tgt.setPixel(i, j, 255);
					else
						tgt.setPixel(i, j, 0);
				}
				else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
					(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3) && (!isCombination))
				{
					if (tgt.getPixel(i, j) > ROI3Median)
						tgt.setPixel(i, j, 255);
					else
						tgt.setPixel(i, j, 0);
				}
			}
		}


	}
	// Applying the Sobel kernel (3*3) over the image's pixels
	float multiply(image& src, int i, int j, bool isXdirection)
	{
		int sum = 0;
		int SobelKernel[3][3] = { {-1, -2, -1},
								   {0, 0, 0},
								   {1, 2, 1} };

		int myMatrix[3][3] = { {src.getPixel(i - 1,j - 1),src.getPixel(i - 1,j),src.getPixel(i - 1,j + 1)},
									   {src.getPixel(i,j - 1),src.getPixel(i ,j),src.getPixel(i ,j + 1)},
									   {src.getPixel(i + 1,j - 1),src.getPixel(i + 1,j),src.getPixel(i + 1,j + 1)} };
		if (isXdirection) {
			int SobelKernel1[3][3] = { 0 };
			for (i = 0; i < 3; ++i) {
				for (j = 0; j < 3; ++j) {
					SobelKernel1[j][i] = SobelKernel[i][j];
				}
			}
			for (i = 0; i < 3; ++i) {
				for (j = 0; j < 3; ++j) {
					SobelKernel[i][j] = SobelKernel1[i][j];
				}
			}
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
					sum += myMatrix[i][j]
					* SobelKernel[i][j];
			}
		}
		return sum;
	}

	// Applying the Sobel kernel (5*5) over the image's pixels
	float multiply5by5(image& src,int i, int j, bool isXdirection)
	{
		int sum = 0;
		int myMatrix[5][5] = { {src.getPixel(i - 2,j - 2),src.getPixel(i - 2,j - 1),src.getPixel(i - 2,j),src.getPixel(i - 2,j + 1),src.getPixel(i - 2,j + 2)},
								   {src.getPixel(i - 1, j - 2),src.getPixel(i - 1,j - 1),src.getPixel(i - 1,j),src.getPixel(i - 1,j + 1),src.getPixel(i - 1,j + 2)},
								   {src.getPixel(i, j - 2),src.getPixel(i ,j - 1),src.getPixel(i ,j),src.getPixel(i ,j + 1),src.getPixel(i ,j + 2)},
								   {src.getPixel(i + 1, j - 2),src.getPixel(i + 1 ,j - 1),src.getPixel(i + 1 ,j),src.getPixel(i + 1 ,j + 1),src.getPixel(i + 1 ,j + 2)},
								   {src.getPixel(i + 2, j - 2),src.getPixel(i + 2 ,j - 1),src.getPixel(i + 2 ,j),src.getPixel(i + 2 ,j + 1),src.getPixel(i + 2 ,j + 2)} };

		int SobelKernel[5][5] = { {-5,-8,-10,-8,-5},{-4,-10,-20,-10,-4},{0,0,0,0,0},{4,10,20,10,4},{5,8,10,8,5} };
		if (isXdirection) {
			int SobelKernel1[5][5] = { 0 };
			for (i = 0; i < 5; ++i) {
				for (j = 0; j < 5; ++j) {
					SobelKernel1[j][i] = SobelKernel[i][j];
				}
			}
			for (i = 0; i < 5; ++i) {
				for (j = 0; j < 5; ++j) {
					SobelKernel[i][j] = SobelKernel1[i][j];
				}
			}
		}
		
		for (int i = 0; i < 5; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				sum += myMatrix[i][j]
					* SobelKernel[i][j];
			}
		}
		return sum;
	}
	/*------*/
	void SobelUtility(image& src, image& tgt, image& tgt2, image& tgt3, int i, int j, int angle, int Threshold,int order,int max) {
		int startValue, endValue;
		if (order == 5) {
			startValue = 1;
			endValue = 3;
		}
		else
		{
			startValue = 0;
			endValue = 1;
		}
		if (i > startValue && j > startValue && i < src.getNumberOfRows() - endValue && j < src.getNumberOfColumns() - endValue) {
			float sumX = 0;
			float sumY = 0;
			if (order == 5) {
				sumX = multiply5by5(src, i, j, true);
				sumY = multiply5by5(src, i, j, false);
			}
			else {
				sumX = multiply(src, i, j, true);
				sumY = multiply(src, i, j, false);
			}
			float lessThanSpecified = 0;
			float MoreThanSpecified = 0;
			if ((angle - 10) < 90) {
				lessThanSpecified = (angle - 10) * ((float)PI / 180);
			}
			else {
				lessThanSpecified = ((90 - angle) - 10) * ((float)PI / 180);
			}
			if ((angle + 10) < 90) {
				MoreThanSpecified = (angle + 10) * ((float)PI / 180);
			}
			else {
				MoreThanSpecified = ((90 - angle) + 10) * ((float)PI / 180);
			}


			float magnitude = sqrt(pow(sumX, 2) + pow(sumY, 2));
			magnitude = ((magnitude / max) * 255);

			tgt3.setPixel(i, j, (int)magnitude);

			if (magnitude > Threshold)
			{
				if (sumX == 0)
				{
					sumX = 0.0000001;
				}
				if (lessThanSpecified < atan(sumY / sumX) && MoreThanSpecified > atan(sumY / sumX))
				{
					tgt2.setPixel(i, j, 255);
				}
				else{
					tgt2.setPixel(i, j, 0);
			}
				tgt.setPixel(i, j,255);
			}
			else {
				tgt2.setPixel(i, j, 0);
				tgt.setPixel(i, j, 0);
			}

		}
	}

	float getMaxmagnitude(image& src,int angle, int order) {
		float max = 0;
		for (int i = 0; i < src.getNumberOfRows(); i++)
		{
			for (int j = 0; j < src.getNumberOfColumns(); j++)
			{
				int startValue, endValue;
				if (order == 5) {
					startValue = 1;
					endValue = 3;
				}
				else
				{
					startValue = 0;
					endValue = 1;
				}
				if (i > startValue && j > startValue && i < src.getNumberOfRows() - endValue && j < src.getNumberOfColumns() - endValue) {
					float sumX = 0;
					float sumY = 0;
					if (order == 5) {
						sumX = multiply5by5(src, i, j, true);
						sumY = multiply5by5(src, i, j, false);
					}
					else {
						sumX = multiply(src, i, j, true);
						sumY = multiply(src, i, j, false);
					}
					float lessThanSpecified = 0;
					float MoreThanSpecified = 0;
					if ((angle - 10) < 90) {
						lessThanSpecified = (angle - 10) * ((float)PI / 180);
					}
					else {
						lessThanSpecified = ((90 - angle) - 10) * ((float)PI / 180);
					}
					if ((angle + 10) < 90) {
						MoreThanSpecified = (angle + 10) * ((float)PI / 180);
					}
					else {
						MoreThanSpecified = ((90 - angle) + 10) * ((float)PI / 180);
					}


					float magnitude = sqrt(pow(sumX, 2) + pow(sumY, 2));
					if (max < magnitude)
					{
						max = magnitude;
					}

				}
			}	
		}
		return max;
	}
	

void utility::SobelEdgeDetector(image& src, image& tgt, image& tgt2, image& tgt3, int order,int arr[]) {
	int x = arr[0], y = arr[1], x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		Sx = arr[2], Sy = arr[3], Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9],Angle = arr[15];
	int Threshold = 105;
	int angle = 45;

	tgt.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt2.resize(src.getNumberOfRows(), src.getNumberOfColumns());
	tgt3.resize(src.getNumberOfRows(), src.getNumberOfColumns());

	float max = getMaxmagnitude(src,angle,order);
	printf("maximum: %f ",max);
	if(order==5){
	for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			tgt.setPixel(i, j, src.getPixel(i, j));
			tgt2.setPixel(i, j, src.getPixel(i, j));
			tgt3.setPixel(i, j, src.getPixel(i, j));
			if (i > 1 && j > 1 && i < src.getNumberOfRows() - 3 && j < src.getNumberOfColumns() - 3) {
				float sumX = 0;
				float sumY = 0;
				if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
				{
					SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh, 5,max);
				}
				else if (((i<x || i<y || i>Sx || j>Sy) &&
					(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2))
				{
					SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh2, 5,max);
				}
				else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
					(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3))
				{
					SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh3, 5,max);
				}
			}
		}
	}
	}	
	else {
		for (int i = 0; i < src.getNumberOfRows(); i++)
		{
			for (int j = 0; j < src.getNumberOfColumns(); j++)
			{
				tgt.setPixel(i, j, src.getPixel(i, j));
				tgt2.setPixel(i, j, src.getPixel(i, j));
				tgt3.setPixel(i, j, src.getPixel(i, j));
				if (i > 0 && j > 0 && i < src.getNumberOfRows() - 1 && j < src.getNumberOfColumns() - 1) {
					float sumX = 0;
					float sumY = 0;
					if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
					{
						SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh, 3, max);
					}
					else if (((i<x || i<y || i>Sx || j>Sy) &&
						(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2))
					{
						SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh2, 3, max);
					}
					else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
						(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3))
					{
						SobelUtility(src, tgt, tgt2, tgt3, i, j, angle, Thresh3, 3, max);
					}

				}
			}
		}
	}
}

/*------*/
void utility::drawHistogram(image& src, image& tgt, image& tgt2, image& tgt3, int arr[])
{
	int x = arr[0], y = arr[1], x2 = arr[5], y2 = arr[6], x3 = arr[10], y3 = arr[11],
		Sx = arr[2], Sy = arr[3], Sx2 = arr[7], Sy2 = arr[8], Sx3 = arr[12], Sy3 = arr[13],
		Thresh = arr[4], Thresh3 = arr[14], Thresh2 = arr[9];
	int count = 1;
	int max = 0,max1=0,max2=0;
	std::map<int, int> mapOfOccur, mapOfOccur1, mapOfOccur2;
	for (int i = 0; i <256; i++)
	{
		mapOfOccur[i] = 0;
		mapOfOccur1[i] = 0;
		mapOfOccur2[i] = 0;
	}


	for (int i = 0; i < src.getNumberOfRows(); i++)
	{
		for (int j = 0; j < src.getNumberOfColumns(); j++)
		{
			//ROI1
			if ((i >= x && i <= Sx) && (j >= y && j <= Sy))
			{
				mapOfOccur[src.getPixel(i, j)] = ++mapOfOccur[src.getPixel(i, j)];
				if (max < mapOfOccur[src.getPixel(i, j)])
				{
					max = mapOfOccur[src.getPixel(i, j)];
				}
			}
			//ROI2
			else if (((i<x || i<y || i>Sx || j>Sy) &&
				(x2 > Sx || y2 > Sy || x2 < x || y2 < y)) && (i >= x2 && i <= Sx2) && (j >= y2 && j <= Sy2) ) {
				mapOfOccur1[src.getPixel(i, j)] = ++mapOfOccur1[src.getPixel(i, j)];
				if (max1 < mapOfOccur1[src.getPixel(i, j)])
				{
					max1 = mapOfOccur1[src.getPixel(i, j)];
				}
			}
			//ROI3
			else if (((i<x2 || i<y2 || i > Sx2 || j > Sy2) &&
				(x3 > Sx2 || y3 > Sy2 || x3 < x2 || y3 < y2)) && (i >= x3 && i <= Sx3) && (j >= y3 && j <= Sy3)) {
				mapOfOccur2[src.getPixel(i, j)] = ++mapOfOccur2[src.getPixel(i, j)];
				if (max2 < mapOfOccur2[src.getPixel(i, j)])
				{
					max2 = mapOfOccur2[src.getPixel(i, j)];
				}
			}
		
		}
	}
	printf("%d,%d,%d \n maxed :", max,max1,max2);
	tgt.resize(max + 1, 256);
	tgt2.resize(max1 + 1, 256);
	tgt3.resize(max2+ 1, 256);
	for (int i = 255; i >= 0; i--)
	{
		if (mapOfOccur[i] > 0)
		{
			for (int j= max; j>= max-mapOfOccur[i];j--)
			tgt.setPixel(j, i, 255);
		}
		if (mapOfOccur1[i] > 0)
		{
			for (int j = max1; j >= max1- mapOfOccur1[i]; j--)
				tgt2.setPixel(j, i, 255);
		}
		if (mapOfOccur2[i] > 0)
		{
			for (int j = max2; j >= max2- mapOfOccur2[i]; j--)
				tgt3.setPixel(j, i, 255);
		}
	}


	/*for (int i = 0, j = 255; i <= 255, j >= 0; i++, j--)
	{

		if (mapOfOccur[i] > 0 && count1 == 0)
		{
			min = i;
			count1 = 1;
		}
		if (mapOfOccur[j] > 0 && count2 == 0)
		{
			max1= j;
			count2 = 1;
		}
		if (count1 != 0 && count2 != 0)
		{
			break;
		}
	}

	printf("%d,%d \n in histogram fn :",+maxvalue,minvalue);*/

}

/*-----------------------------------------------------------------------**/
void utility::scale(image &src, image &tgt, float ratio)
{
	int rows = (int)((float)src.getNumberOfRows() * ratio);
	int cols  = (int)((float)src.getNumberOfColumns() * ratio);
	tgt.resize(rows, cols);
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{	
			/* Map the pixel of new image back to original image */
			int i2 = (int)floor((float)i/ratio);
			int j2 = (int)floor((float)j/ratio);
			if (ratio == 2) {
				/* Directly copy the value */
				tgt.setPixel(i,j,checkValue(src.getPixel(i2,j2)));
			}

			if (ratio == 0.5) {
				/* Average the values of four pixels */
				int value = src.getPixel(i2,j2) + src.getPixel(i2,j2+1) + src.getPixel(i2+1,j2) + src.getPixel(i2+1,j2+1);
				tgt.setPixel(i,j,checkValue(value/4));
			}
		}
	}
}
