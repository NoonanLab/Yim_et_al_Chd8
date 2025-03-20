#!/bin/sh
### Minimum requirements to run script: ImageMagick(6.8)
# Input: cropped rectangular strips of full IHC images, spanning the apicobasal axis of the developing cortex

mkdir sliced
for f in $PWD/*.tif; do
    name=`echo "$f" | sed s/\.tif$//`
    filename=$(basename $f .tif)
convert ${f} -resize "1000x1000>" -crop 100%x10% +repage "$PWD/sliced/${filename}.jpg";
done
cp count.sh sliced/
cd sliced
sh count.sh
exit
