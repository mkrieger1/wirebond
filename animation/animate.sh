#!/bin/sh

# place postscript files named "bondspstest_xxx.ps" (xxx: sequentially
# numbered) in this directory and run this script to create an animated GIF


for file in bondspstest*.ps
do
    echo converting $file
    gs -sDEVICE=png16m -sOutputFile=$file.png -r150 -dGraphicsAlphaBits=4 \
       -dBATCH -dNOPAUSE $file $2>/dev/null
done

echo creating animated GIF...
convert bondspstest*.png -trim +repage animation.gif

echo cleaning up...
rm bondspstest*.png
rm bondspstest*.ps

