#!/bin/bash
MYDATE=`date +%d%m%y`
LOGFILE=`echo "$4.log.$MYDATE"`



Flanking_region_option=${19}
Interval_size=${21}

if [ "$Flanking_region_option" == "flanking_region_size" ] 
then
if [ "$Interval_size" == "automatic" ] 
then
ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R bed -E $2 -C $3 -O $4 -O2 $5  -O3 $6   -T "$7"  -D $8     -S  $9  -GO  ${10}  -CS  ${11}  -FL  ${12}  -MQ   ${13} -SE   ${14}  -RB   ${15}  -FC  ${16}  -MW    ${17}  -H   ${18}  -L ${20}   >$LOGFILE 2>&1
else  ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R bed -E $2 -C $3 -O $4 -O2 $5  -O3 $6   -T "$7"  -D $8     -S  $9  -GO  ${10}  -CS  ${11}  -FL  ${12}  -MQ   ${13} -SE   ${14}  -RB   ${15}  -FC  ${16}  -MW    ${17}  -H   ${18}  -L ${20}  -I ${21}  >$LOGFILE 2>&1
fi
fi

if [ "$Flanking_region_option" == "flanking_floating_size" ]
then
if [ "$Interval_size" == "automatic" ] 
then
ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R bed -E $2 -C $3 -O $4 -O2 $5  -O3 $6   -T "$7"  -D $8     -S  $9  -GO  ${10}  -CS  ${11}  -FL  ${12}  -MQ   ${13} -SE   ${14}  -RB   ${15}  -FC  ${16}  -MW    ${17}  -H   ${18}  -N ${20}   >$LOGFILE 2>&1
else  ngs.plot.r  -Galaxy 1 -P 0 -G $1 -R bed -E $2 -C $3 -O $4 -O2 $5  -O3 $6   -T "$7"  -D $8     -S  $9  -GO  ${10}  -CS  ${11}  -FL  ${12}  -MQ   ${13} -SE   ${14}  -RB   ${15}  -FC  ${16}  -MW    ${17}  -H   ${18}  -N ${20}  -I ${21}  >$LOGFILE 2>&1
fi
fi



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

