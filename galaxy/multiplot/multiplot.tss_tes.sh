#!/bin/bash
MYDATE=`date +%d%m%y`
LOGFILE=`echo "$4.log.$MYDATE"`

ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6    -D $7   -L $8  -S  $9  -GO  ${10}  -CS  ${11}  -FL  ${12}  -MQ   ${13} -SE   ${14}  -RB   ${15}  -FC  ${16}  -MW    ${17}  -H   ${18}    >$LOGFILE 2>&1
if (grep -q "Error" $LOGFILE)||(grep -q "this is not a BAM file" $LOGFILE)
then
  $NGSPLOT/galaxy/ngsplot/redirect.pl $LOGFILE
fi

zip  -q  data.zip  data/*
mv  data.zip  $4

#new_avg=`echo "$5.png"`
#new_heatmap=`echo "$6.png"`
#convert  $5  $new_avg
#convert  $6  $new_heatmap 

#mv $new_avg  $5  
#mv $new_heatmap  $6

