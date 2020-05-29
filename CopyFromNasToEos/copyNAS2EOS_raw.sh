#!/bin/bash

# copyNAS2EOS_dat.sh
#
# Created by: Kalpanie Madara Liyanage
#
# Usage: Copy QC8 uFEDkit run files from NAS to EOS

RUN=$(printf "%06d" $1)
DIRECTORY=/eos/cms/store/group/dpg_gem/comm_gem/QC8_Commissioning/run$RUN
src=/data/bigdisk/GEM-Data-Taking/GE11_QC8/Cosmics/run$RUN/run$RUN
COUNT=0
COPY=0
size=0

if [ ! -d "$DIRECTORY" ]; then
echo "Creating Directory: "$DIRECTORY
mkdir $DIRECTORY
else
echo "Existing Directory: "$DIRECTORY
fi

#Get the list of files to send
FILES=$(ssh gemuser@gem904qc8daq ls ${src}"_ls"*".raw")
max=-1
#Tell the user which group of files will be copied
echo "The list of files to be copied is:"
for f in $FILES\;
do
echo $f

i=$(echo "$f" | cut -d'_' -f 3)
n=${i:3}
n=$(expr $n + 0)
[[ $n -gt $max ]] && max=$n
done
echo "Number of lumi sections: $max"

MODTIME=$(ssh gemuser@gem904qc8daq date -r ${src}"_ls0001_index000000.raw" +"%Y-%m-%d_%H-%M-%S")
#echo $MODTIME

#Copy the files
echo "Copying files"

for (( l=1; l<=$max; l++ ))
do
((COUNT++))
echo "Processing lumi section: " $l
j=$(printf "%04d" $l)
name=run${RUN}"_Cosmic_CERNQC8_"${MODTIME}"_ls"${j}"_allindex.raw"
lines=$(ssh gemuser@gem904qc8daq cat ${src}_ls${j}_index*.raw | wc -l)
echo "Number of lines in flies of lumi${j} : "$lines
ssh gemuser@gem904qc8daq cat ${src}"_ls"${j}"_index"*".raw" >> $name
echo 'Processed lumi section' $l
echo $name
all_lines=$(cat $name | wc -l)
echo "Number of lines in $DIRECTORY/$name : "$all_lines

if [ $lines -eq $all_lines ]; then
echo "All files have been merged!"
scp $name $DIRECTORY/$name
else
echo "Some files have not been merged!"
fi
rm $name
done

#Confirm for the user that the files were copied
echo "Showing the contents of " $DIRECTORY
ls -1 | wc -l $DIRECTORY"/"*"_allindex.raw"
COPIED=$(ls $DIRECTORY"/"*"_allindex.raw")

for c in $COPIED;
do
((COPY++))
done

if [ $COUNT -eq $COPY ]; then
echo "All files have been copied!"
else
echo "Some files are missing!"
fi
