#!/bin/bash
MYDATE=`date +%d%m%y`
LOGFILE=`echo "$3.log.$MYDATE"`

cp $2  data.zip

replot.r  $1 -I data.zip   -O  $3  -SE  $4  -MW  $5  -H $6  -P 12  >$LOGFILE 2>&1

if (grep -q "Error" $LOGFILE)||(grep -q "this is not a BAM file" $LOGFILE)
then
  $NGSPLOT/galaxy/ngsplot/redirect.pl $LOGFILE
fi
old_replot=`echo "$3.pdf"`
#new_replot=`echo "$3.png"`

#convert  $old_replot  $new_replot
 
#mv $new_replot  $3
#rm $old_replot
mv  $old_replot  $3


