#!/bin/bash
MYDATE=`date +%d%m%y`
LOGFILE=`echo "$4.log.$MYDATE"`


Flanking_region_option=${18}
Interval_size=${20}

if [ "$Flanking_region_option" == "flanking_region_size" ] 
then
if [ "$Interval_size" == "automatic" ] 
then
ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6    -D $7     -S  $8  -GO  $9  -CS  ${10}  -FL  ${11}  -MQ   ${12} -SE   ${13}  -RB   ${14}  -FC  ${15}  -MW    ${16}  -H   ${17}   -F ${19} -L ${21}  >$LOGFILE 2>&1
else  ngs.plot.r  -Galaxy 1 -P 0  -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6    -D $7     -S  $8  -GO  $9  -CS  ${10}  -FL  ${11}  -MQ   ${12} -SE   ${13}  -RB   ${14}  -FC  ${15}  -MW    ${16}  -H   ${17}   -F ${19}   -L ${21}  -I ${20}  >$LOGFILE 2>&1
fi
fi

if [ "$Flanking_region_option" == "flanking_floating_size" ]
then
if [ "$Interval_size" == "automatic" ] 
then
ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6    -D $7     -S  $8  -GO  $9  -CS  ${10}  -FL  ${11}  -MQ   ${12} -SE   ${13}  -RB   ${14}  -FC  ${15}  -MW    ${16}  -H   ${17}  -F ${19}   -N ${21}   >$LOGFILE 2>&1
else  ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R $2 -C $3 -O $4 -O2 $5  -O3 $6    -D $7     -S  $8  -GO  $9  -CS  ${10}  -FL  ${11}  -MQ   ${12} -SE   ${13}  -RB   ${14}  -FC  ${15}  -MW    ${16}  -H   ${17}   -N ${21} -F ${19} -I ${20}  >$LOGFILE 2>&1
fi
fi



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

