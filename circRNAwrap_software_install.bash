#!/bin/bash

# test whether install the software


command='awk'
m=$(command -v $command)
if [ -z "$m" ]; then
    echo "no, not install   $command"
else  
    echo "yes, install well $command in" command -v $command	
fi

readonly m
unset m
echo $m

for i in 1 2 3 
  do 
    echo $i
done

for i in 1 2 3;do  echo $i;done