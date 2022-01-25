###Image segmentation part####
###Using package EBImage : install via Bioconductor###
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("EBImage")
###Things to load###
library(Momocs)
library(EBImage)
source("/home/samuel/Documents/Rfunctions1.txt") #Text file containing code of function from Claude (2008). Morphometrics with R.
##########################################################################################

##########################
### Floodfill approach ###
##########################

###Works with background of heterogeneous colors###

setwd("Bureau/collab_Jules_papillons/test/gray") #Working directory contains grayscale images only.

#List of files
files <- list.files()

#Create lists to fill with output from loop
Gray_images <- vector("list", length=length(files))
BW_images <- vector("list", length=length(files))
BW_final <- vector("list", length=length(files))
shapes.1000 <- vector("list", length=length(files))

start <- Sys.time() #Starting time of image treatment

for (i in 1:length(files)) { #Start of loop image treatment

   message(i, " ", files[i]) #To check advancement of loop

   x2 <- readImage(files[i]) #Input grayscale image, with any type of background

#Floodfill background, starting from top-left-most pixel (assuming it is part of background). Some aliasing will happen for clear backgrounds
   e <- floodFill(x2, pt=c(1,1), col="black", tol=0.1) #Default tolerance = 0 not good for heterogeneous background.

#Convert grayscale with black background to B/W image
   f <- bwlabel(e)

#Produce vector with all 'objects' in image, with corresponding number of px.
   tf <- table(f)[-1] #First object is background, remove from vector

#Index of all objects to remove, i.e. all, except background and object with max number of px (hopefully the shape of interest).
   obj_rm <- c(1:max(f))[-which.max(tf)] 

#Remove parasite objects
   frm <- rmObjects(f, obj_rm)

   Gray_images[[i]] <- e 
   BW_images[[i]] <- f
   BW_final[[i]] <- frm #Keeping these lists in memory is a bad idea... ~1GB for 5 pictures...

} #End of image treatment loop.

end <- Sys.time() #Ending time of image treatment
runtime <- end-start #Total time

##########################################
### Finding points of contour of shape ###
##########################################

###With Momocs###
for (i in 1:length(BW_final)) {
#Negative (for Momocs, object is black on white background)
   nh <- max(BW_final[[i]])-BW_final[[i]]

#Extract matrix from 'Image' object to switch from EBImage to Momocs
   img <- t(imageData(nh))

#Find first black pixel in matrix = start point within shape
   x <- which(img==0, arr.ind=T)[1,]

#Extract contour, using Momocs function 'import_Conte'
   contour1<-import_Conte(img=img, x=x)

#Subsample points
   shapes.1000[[i]] <- coo_sample(contour1, 200)
}


##########################
### Synthetic function ###
##########################

# x : vector with file names, eg list.files()[1:10]
# NLM : Number of landmarks to represent the contour
# tol : tolerance level for floodfilling step

auto.LM.cont <- function(x=list.files(), NLM=200, tol=0.1) {
   shapes <- vector("list", length=length(x))
   start <- Sys.time() #Starting timer
   for (i in 1:length(x)) { #Start of loop image treatment
      message(i, " ", x[i]) #To check advancement of loop
      x2 <- readImage(x[i]) #Input grayscale image, with any type of background
#Floodfill background, starting from top-left-most pixel.
      e <- floodFill(x2, pt=c(1,1), col="black", tol=tol)
#Convert grayscale with black background to B/W image
      f <- bwlabel(e)
#Produce vector with all 'objects' in image, with corresponding number of px.
      tf <- table(f)[-1] #First object is background, remove from vector
#Index of all objects to remove, i.e. all, except background and object with max number of px.
      obj_rm <- c(1:max(f))[-which.max(tf)] 
#Remove parasite objects
      frm <- rmObjects(f, obj_rm)
#Negative (for Momocs, object is black on white background)
      nh <- max(frm)-frm
#Extract matrix from 'Image' object to switch from EBImage to Momocs
      img <- t(imageData(nh))
#Find first black pixel in matrix = start point within shape
      first.px <- which(img==0, arr.ind=T)[1,]
#Extract contour, using Momocs function 'import_Conte'
      contour <- import_Conte(img=img, x=first.px)
#Subsample points
      shapes[[i]] <- coo_sample(contour, NLM)
   } #End of 'for' loop

   end <- Sys.time() #Ending timer
   runtime <- end-start
   list(runtime=runtime, shapes=shapes)

} #End of function

