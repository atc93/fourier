#!/bin/bash

# Argument 1: full paht directory with jobs output

cd $1

echo "Working in $1"
 
		 nFile=`ls | wc -l`
 firstFile=`ls  | sed '1!d'`
	lastFile=`ls  | sed ${nFile}!'d'`

echo "Number of file to process: $nFile"
echo "First file to process: $firstFile"
echo "Last file to process:  $lastFile"
echo "Is that righ? Press any key if yes"
read 
echo ""
echo " -- PROCESSING -- "

/nfs/gm2/data1/achapela/CodeForBMAD/bin/textFileToRootFile $firstFile $lastFile $1
