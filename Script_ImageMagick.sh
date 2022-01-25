#Image Magick treatment

#Fred's ImageMagick Scripts : 'corners'
#Fred's ImageMagick Scripts : 'multicrop2'
#Fred's ImageMagick Scripts : 'multicrop'

#If some names have spaces use :

find -name "* *" -type f | rename 's/ /_/g'

#Make directory to hold splitted images

mkdir splitted

#Loop to split images with several objects (e.g. 2 views of same butterfly + scale bar)
#With timer. 
#Options : -d = any object smaller is discarded, 
#-u = rotation (3, no rotation), 
#-e = number of px added to bounding box of the objects.

time for f in *.jpg
   do
   echo $f
   multicrop2 -d 2000 -u 3 -e 10 $f split_$f.jpg
done

mv split_* splitted


time for f in *.jpg
   do
   echo $f
   multicrop2 -f 20 -d 200 -u 3 -e 5 $f split2_$f
done

mv split2_* splitted2

#Different version, first one stopped at image n°4
mkdir splitted

time for f in *.jpg
   do
   echo $f
   multicrop -d 200 -u 3 -e 5 $f split2_$f
done

mv split2_* splitted2

cd splitted2
find *.jpg -size -100k -delete

time for f in *.jpg
   do
   echo $f
   multicrop -c West -f 15 -u 3 -e 5 $f split2_$f
done


mv split2_* splitted2

cd splitted2
find *.jpg -size -100k -delete

#Renaming
for file in * ; do
    mv -v "$file" "${file#*_}"
done

for f in preview.[1-9]* #Specific 0000 padding for image downloaded from London NHM
do
   mv $f `printf preview.%04d ${f#preview.}`
done

for file in prev*
do read line
echo mv -v "${file}" "${line}"
done < ~/Bureau/collab_jules/data/data_londres/Birdwing_dataset/IDs_img3.txt

for file in *
do
 mv -v "${file}" "${file}.jpg"
done
#Change to grayscale
cp -a splitted gray
cd gray

time mogrify -type Grayscale *.jpg 



###Script for Morpho dataset from V. Debat

cd mnt/832c137a-8955-4ded-b866-7120208caee6/Heavy\ stuff/Images/
mkdir all_images_morpho
find Filemail.com\ files\ 2021-1-26\ ilbdxzbkmbutvky -type f -print0 | xargs -0 mv -t all_images_morpho/
find Filemail.com\ files\ 2021-1-26\ urwdghorqppshkl -type f -print0 | xargs -0 mv -t all_images_morpho/
cd all_images_morpho/

mogrify -format jpg *.tif
rm *-1.jpg
mogrify -gravity East -crop 86x100%+0+0 *.jpg #86% is not enough, a few butterflies are cutoff
time mogrify -type Grayscale *.jpg 

mkdir jpgs
mv *.jpg jpgs/


