# Image-Treatment-Data-Extraction
Scripts in development implementing some (semi-)automatic computer vision approaches to extracting data from scientific images.

The scripts are developed based on a large dataset of high-resolution images of butterflies downloaded from the London National History Museum (directly from the website) and Museum d'Histoire Naturelle de Toulouse.

The shell / ImageMagick script "Script_ImageMagick.sh" is aimed at splitting images containing several objects of interest into images containing single objects of interest, and convert them to grayscale images.
The R script "Streamlind_papillons.R" contains custom functions to segment grayscale images into foreground and background, remove small unwanted objects, and turn them into binary (black / white) images.
The R script "Papillons_morphometry.R" contains R script based mainly on the Momocs package, and deals with extracting contours from binary images, and (after checking and removing "bad" contours) automatically place landmarks and analyze contours.
