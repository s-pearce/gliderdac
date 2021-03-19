#!/bin/bash

# auto_ftp_tx.sh  Auto GliderDAC ftp transfer
#   given a gliderid
#

app=$(basename $0)
usage="Usage: $app GLIDER_REF_DESIGNATOR DEPLOYNUM"
if [ $# -ne 2 ]; then
    echo "Requires 2 inputs" >&2
    echo "$usage" >&2
    exit 1
fi

# make sure GLIDERID is ref designator ce05moas-gl### with regex
if [[ $1 =~ ^ce05moas-gl[0-9]{3}$ ]]
then
    GLIDER_RD=$1
else 
    echo "GLIDERID must be of the form ce05moas-gl###" >&2
    echo "$usage" >&2
    exit 1
fi

# make sure DEPLOYNUM is (R|D)##### with regex
if [[ $2 =~ ^(R|D)([0-9]{5})$ ]]
then
    DEPLOYMENT=$2
else
    echo "DEPLOYNUM must be of the form (R|D)#####" >&2
    echo "$usage" >&2
    exit 1
fi

PROCDIR=/home/ooiuser/data/proc
NC_DIR=$PROCDIR/$GLIDER_RD/$DEPLOYMENT/nc
CONFIG_DIR=$PROCDIR/$GLIDER_RD/$DEPLOYMENT/config
LOG_DIR=$PROCDIR/$GLIDER_RD/$DEPLOYMENT/logs
CODE_DIR=/home/ooiuser/code/cgsn-gliders/
FTP_SCRIPT=$CODE_DIR/autodacscripts/dac_ftp_uploads.py

if [ ! -d "$NC_DIR" ]
then
    echo "Directory $NC_DIR does not exist, check your inputs" >&2 
    echo "or add new glider directory" >&2
    exit 1
fi

# Log the uploading results to an FTP log
LOG=$LOG_DIR/${GLIDER_RD}-${DEPLOYMENT}_rtime_ftp.log
# dump date, time, and script call into the log
echo `date +"%F %T"` >> $LOG
echo "auto_ftp_tx.sh $GLIDER_RD $DEPLOYMENT" >> $LOG

#source activate gdac
#try this
gdac_python=/home/cgsnmo/anaconda3/envs/gdac/bin/python
$gdac_python $FTP_SCRIPT $CONFIG_DIR $NC_DIR 2>&1 | tee -a $LOG

# add a blank line between runs
echo "" >> $LOG