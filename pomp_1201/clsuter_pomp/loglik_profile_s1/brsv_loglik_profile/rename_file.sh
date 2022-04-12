#!/bin/bash

## rename file
for old in hhs*_arsv*
do 
	new=$(echo $old | sed 's/arsv/brsv/') 
	mv -v "$old" "$new" 

done