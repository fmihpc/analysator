#!/bin/bash

/proj/ykempf/visit3.3.3-build/bin/visit -nowin -cli -np 64 -l srun -la "--mpi=pmix_v3" -s visit_dayside_cutoff.py

/proj/ykempf/visit3.3.3-build/bin/visit -nowin -cli -np 64 -l srun -la "--mpi=pmix_v3" -s visit_nightside_cutoff.py

module load ImageMagick

#fc-cache -vf

export font="DejaVu-Sans"
export ptsize=65

arrow_head="path 'M 0,0  l -30,-30, +80,+30 -80,+30 +30,-30 z'"

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'(a) <i>R</i><span size="35000" rise="-6000">cutoff</span>=3<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +115+0 \
    -  FHA_dayside_fluxropes_cutoff_1600_3.0.png   tmp0.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp0.png   tmp1.png

convert tmp1.png -gravity east -extent 3280x1800 tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Z</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png  tmp3.png

convert -composite tmp3.png  xc: \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 160,1360 300,1360' -strokewidth 5 -draw "translate 300,1360 rotate 0 $arrow_head " \
FHA_dayside_fluxropes_cutoff_1600_3.0_labelled.png



magick -background white -fill black -font $font -pointsize $ptsize pango:'(b) <i>R</i><span size="35000" rise="-6000">cutoff</span>=5<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +115+0 \
    -  FHA_dayside_fluxropes_cutoff_1600_5.0.png   tmp0.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp0.png   tmp1.png

convert tmp1.png -gravity east -extent 3280x1800 tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Z</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png   tmp3.png

convert -composite tmp3.png  xc: \
-stroke black -strokewidth 10 -draw 'line 1100,170 1200,270' -strokewidth 5 -draw "translate 1200,270 rotate 45 $arrow_head " \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 160,1360 300,1360' -strokewidth 5 -draw "translate 300,1360 rotate 0 $arrow_head " \
  FHA_dayside_fluxropes_cutoff_1600_5.0_labelled.png



magick -background white -fill black -font $font -pointsize $ptsize pango:'(c) <i>R</i><span size="35000" rise="-6000">cutoff</span>=7<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +115+0 \
    -  FHA_dayside_fluxropes_cutoff_1600_7.0.png   tmp0.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp0.png   tmp1.png

convert tmp1.png -gravity east -extent 3280x1800 tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Z</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png   tmp3.png

convert -composite tmp3.png  xc: \
-stroke black -strokewidth 10 -draw 'line 1100,170 1200,270' -strokewidth 5 -draw "translate 1200,270 rotate 45 $arrow_head " \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 160,1360 300,1360' -strokewidth 5 -draw "translate 300,1360 rotate 0 $arrow_head " \
  FHA_dayside_fluxropes_cutoff_1600_7.0_labelled.png




magick -background white -fill black -font $font -pointsize $ptsize pango:'(a) <i>R</i><span size="35000" rise="-6000">cutoff</span>=3<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +142+0 \
    -  FHA_nightside_fluxropes_cutoff_1600_3.0.png   tmp0.png

convert tmp0.png -gravity northeast -extent 3280x1400 tmp1.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>X</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp1.png   tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png   tmp3.png

convert -composite tmp3.png  xc: \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 2630,1330 2630,1180' -strokewidth 5 -draw "translate 2630,1180 rotate -90 $arrow_head " \
-stroke SpringGreen2 -fill SpringGreen2 -strokewidth 10 -draw 'line 1845,1350 1770,1200' -strokewidth 5 -draw "translate 1770,1200 rotate -120 $arrow_head " \
FHA_nightside_fluxropes_cutoff_1600_3.0_labelled.png



magick -background white -fill black -font $font -pointsize $ptsize pango:'(b) <i>R</i><span size="35000" rise="-6000">cutoff</span>=5<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +142+0 \
    -  FHA_nightside_fluxropes_cutoff_1600_5.0.png   tmp0.png

convert tmp0.png -gravity northeast -extent 3280x1400 tmp1.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>X</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp1.png   tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png   tmp3.png

convert -composite tmp3.png  xc: \
-stroke black -strokewidth 10 -draw 'line 2850,750 3000,750' -strokewidth 5 -draw "translate 3000,750 rotate 0 $arrow_head " \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 2630,1330 2630,1180' -strokewidth 5 -draw "translate 2630,1180 rotate -90 $arrow_head " \
-stroke SpringGreen2 -fill SpringGreen2 -strokewidth 10 -draw 'line 1845,1350 1770,1200' -strokewidth 5 -draw "translate 1770,1200 rotate -120 $arrow_head " \
FHA_nightside_fluxropes_cutoff_1600_5.0_labelled.png



magick -background white -fill black -font $font -pointsize $ptsize pango:'(c) <i>R</i><span size="35000" rise="-6000">cutoff</span>=7<i>R</i><span size="35000" rise="-6000">c</span>' -compose Copy -bordercolor Black -border 5 miff:- |\
  magick composite -gravity northwest -geometry +142+0 \
    -  FHA_nightside_fluxropes_cutoff_1600_7.0.png   tmp0.png

convert tmp0.png -gravity northeast -extent 3280x1400 tmp1.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>X</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -compose Copy miff:- |\
  magick composite -gravity south -geometry +0+0 \
    -  tmp1.png   tmp2.png

magick -background white -fill black -font ${font} -pointsize $ptsize pango:'<i>Y</i> (<i>R</i><span size="35000" rise="-6000">E</span>)' -rotate -90 -compose Copy miff:- |\
  magick composite -gravity west -geometry +0+0 \
    -  tmp2.png   tmp3.png

convert -composite tmp3.png  xc: \
-stroke black -strokewidth 10 -draw 'line 2850,750 3000,750' -strokewidth 5 -draw "translate 3000,750 rotate 0 $arrow_head " \
-stroke grey40 -fill grey40 -strokewidth 10 -draw 'line 2630,1330 2630,1180' -strokewidth 5 -draw "translate 2630,1180 rotate -90 $arrow_head " \
-stroke SpringGreen2 -fill SpringGreen2 -strokewidth 10 -draw 'line 1845,1350 1770,1200' -strokewidth 5 -draw "translate 1770,1200 rotate -120 $arrow_head " \
FHA_nightside_fluxropes_cutoff_1600_7.0_labelled.png


convert -append FHA_dayside_fluxropes_cutoff_1600_?.0_labelled.png FHA_dayside_fluxropes_cutoff_1600_labelled.png
convert -append FHA_nightside_fluxropes_cutoff_1600_?.0_labelled.png FHA_nightside_fluxropes_cutoff_1600_labelled.png

