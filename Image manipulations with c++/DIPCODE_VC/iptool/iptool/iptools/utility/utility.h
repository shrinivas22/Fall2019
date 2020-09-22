#ifndef UTILITY_H
#define UTILITY_H

#include "../image/image.h"
#include <sstream>
#include <math.h>

class utility
{
	public:
		utility();
		virtual ~utility();
		static std::string intToString(int number);
		static int checkValue(int value);
		static void addGrey(image &src, image &tgt, int value);
		static void histogramStretching(image& src, image& tgt, image& tgt2, image& tgt3,image& tgt4, int a, int b, int a1, int b1, int a2, int b2, int arr[], bool isCombination);
		static void binarize(image &src, image &tgt, int threshold);
		static void drawHistogram(image& src, image& tgt, image& tgt2, image& tgt3, int arr[]);
		static void optimalthresholding(image& src, image& tgt, image& combinedBlackImg, image& combinedWhiteImg, int arr[],bool isCombination);
		static void Combination(image& src, image& tgt, image& bg, image& fg, image& bgCombination, image& fgCombination, image& bgHist, image& fgHist, image& combined, int arr[]);
		static int findMaxMin(image& src, int arr[],bool isMax);
		static void SobelEdgeDetector(image& src, image& tgt, image& tgt2, image& tgt3, int order,int arr[]);
		static void sobelUtility(image& src, image& tgt, image& tgt2, image& tgt3, int order,int max);
		static int findMedian(int arr1[], int size);
		static void scale(image &src, image &tgt, float ratio);
};

#endif

