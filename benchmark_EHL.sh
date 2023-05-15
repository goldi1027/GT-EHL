#!/bin/bash

for file in dataset/merged-mesh/$1/* ;do
  echo Processing query for $file ..........
  array=(${file//// });
  directory_name=${array[2]};
  f="$(basename -- $file)";
  array2=(${f//-/ });
  map_name=${array2[0]};
  echo $map_name

./bin/testEHL $directory_name $map_name 1
./bin/testEHL $directory_name $map_name 4
./bin/testEHL $directory_name $map_name 16
./bin/testEHL $directory_name $map_name 64

done