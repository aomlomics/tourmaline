#!/bin/bash
while read line
do
    if [[ $line =~ ^\>.* ]]
    then
        echo $line | sed 's/>//' | tr -d '\n'
        echo -e -n '\t'
    else
        echo $line | sed 's/[^-]//g' | awk '{{ print length }}'
    fi
done < "${1:-/dev/stdin}"
