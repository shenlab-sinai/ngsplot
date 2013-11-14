#!/bin/bash
MYDATE=`date +%d%m%y`
LOGFILE=`echo "$4.log.$MYDATE"`


ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6   -T "$7"  -D $8   -L $9  -S  ${10}  -GO  ${11}  -CS  ${12}  -FL  ${13}  -MQ   ${14} -SE   ${15}  -RB   ${16}  -FC  ${17}  -MW    ${18}  -H   ${19}    >$LOGFILE 2>&1
if (grep -q "Error" $LOGFILE)||(grep -q "this is not a BAM file" $LOGFILE)
then
  $NGSPLOT/galaxy/ngsplot/redirect.pl $LOGFILE
fi

zip  -q  data.zip  data/*
mv  data.zip  $4

#new_avg=`echo "$5.pdf"`
#new_heatmap=`echo "$6.pdf"`
#convert  $5  $new_avg
#convert  $6  $new_heatmap 

#mv $new_avg  $5  
#mv $new_heatmap  $6

