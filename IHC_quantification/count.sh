#!/bin/sh
### Minimum requirements to run script: ImageMagick (6.8) and Perl (5) installation
# Input: Color image files (jpg format), located in the same directory as the count.sh and process.pl scripts; this script can handle single-channel images and merged channel images as inputs

echo "" >> data.csv
date >> data.csv 
for file in *.jpg; do
name=`echo "$file" | sed s/\.jpg$//`
convert "$name.jpg" -gaussian-blur 3 -separate -contrast-stretch 2% -combine -set option:convolve:scale '50%!' -bias 50% \( -clone 0 -morphology Convolve "1,2,1,0,0,0,-1,-2,-1" \) \( -clone 0 -morphology Convolve "1,0,-1,2,0,-2,1,0,-1" \) -delete 0 -solarize 50% -level 50,0% +level 0,70% -gamma 0.5 -compose plus -composite -gamma 2 -auto-level txt:- | tail -n +2 | tr -cs '0-9\n'  ' ' | while read x y r g b junk; do 
	if [ $r -gt 0 ] || [ $g -gt 0 ] || [ $b -gt 0 ] 
	then
		echo "$r,$g,$b" >> $name-pixels.csv
	fi
done
perl process.pl "$name-pixels.csv" "data.csv" &
#rm "$name-pixels.csv"
done
