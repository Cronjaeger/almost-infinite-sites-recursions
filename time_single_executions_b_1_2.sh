#!/bin/bash
PATH=./testdata/fixed-number-of-mutations/

for B in 1 2; do
  OUTPUT_FILE=timing-output-with-$B-extra-mutations.txt
  /bin/touch $OUTPUT_FILE
  for FILE in $PATH*.csv; do
    #/usr/bin/python benchmarking_singleExecution.py $FILE $B >> OUTPUT_FILE
    echo '================================================================================' | /usr/bin/tee -a $OUTPUT_FILE
    /usr/bin/python benchmarking_singleExecution.py $FILE $B | /usr/bin/tee -a $OUTPUT_FILE
  done
done
