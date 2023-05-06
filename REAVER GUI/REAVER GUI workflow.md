# Recommended workflow for segmentation with REAVER GUI

## Disclaimer:
* Man-made segmentation is always the best segmentation- 
if you can't see the vessel the computer can't either

Introduction to the original SW: *[REAVER](https://github.com/uva-peirce-cottler-lab/public_REAVER)*

## What are the changes from the original SW?
* Works on any type of tiff images (8, 10 or 16 bit)
* Added an extra hyperparameter- "background subtraction factor"- 
which multiplies the background image before its subtracted from the 8-pixel-averaged image
* Replaced the "grey neighbors" image display with the "background subtracted" image.
* Added an extra morphological "closing" after vessel filling (before the compensatory dilation the BW image)
* The curosor size adapts to zooming

<img src="../resources/REAVER GUI start screen.png">
## Recommended workflow
1. Load the tiff files (File > Load directory)
2. Open a representative image (by pressing the image name in the files menu on the left)
3. Choose the FITC (green) channel as input chnnale (right menu)
4. Open the gearbox tool and specify the default hyperparameters:
- Averaging filter size = ~1.5 times the diameter of the largest vessel (in px)
- Minimum connected components area = slightly less than the expected smallest vessel cross-section (in px)
for example, if the smallest vessel radius is 6 pixels choose (pi*6^2) or approx. 100 (px)
- Wire dilation = 0
- Vessel thickness threshold = 0
- background subtraction factor = 0.5
5. Choose an initial grey threshold = 0.06 (on the right menu)
6. Segment the image (click the left-side button)

* This will create an initial segmentation, which can be improved iteratively (a new segmentation should be done for each iteration):
1. Look at the grey image and the background subtraction image (left radio-button menu)
- The background subtraction image should be better then the gray image by:
  - having an enhanced SNR
  - having sharper vessel edges
  - having more uniform vessel intenity across different regions of the frame

<img src="../resources/REAVER GUI bkg sub screen.png">
ways to improve the background subtraction image:
- **If the SNR is lower than in the grey image:** increase the averaging filter size / decrease subtraction factor
- **If only the vessel edges are seen/ or only some vessels are seen:** increase the averaging filter size / decrease subtraction factor

Once you are happy with the background subtraction image (all vessels are clearly brighter than the background):
2. Look at the First binary image
- The vessels should all be segmented in white, some small amount of scattered noise should also be segmented in white

<img src="../resources/REAVER GUI 1st bin screen.png">
ways to improve the First binary image:
- **If there is a lot of scattered noise in the segmentation:** increase the grey threshold
- **If some of the vessels are not segmeted:** decrease the grey threshold

Now that you are happy with the first binary segmentation, the 2nd binary image should look good as well.

<img src="../resources/REAVER GUI seg screen.png">
From here on its all about manual segmentation fixes:
