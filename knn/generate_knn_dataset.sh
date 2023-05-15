#!/bin/bash
for file in dataset/merged-mesh/$1/* ;do
  echo Generating knn dataset for $file ..........
  array=(${file//// });
  directory_name=${array[2]};
  f="$(basename -- $file)";
  array2=(${f//-/ });
  map_name=${array2[0]};

#  ./bin/knn_dataset_generator $directory_name $map_name 0.001
#  ./bin/knn_dataset_generator $directory_name $map_name 0.01
#  ./bin/knn_dataset_generator $directory_name $map_name 0.1

  ./bin/knn_dataset_generator $directory_name $map_name 0.001 10 1
  ./bin/knn_dataset_generator $directory_name $map_name 0.001 10 4
  ./bin/knn_dataset_generator $directory_name $map_name 0.001 10 16
  ./bin/knn_dataset_generator $directory_name $map_name 0.01 10 1
  ./bin/knn_dataset_generator $directory_name $map_name 0.01 10 4
  ./bin/knn_dataset_generator $directory_name $map_name 0.01 10 16
  ./bin/knn_dataset_generator $directory_name $map_name 0.1 10 1
  ./bin/knn_dataset_generator $directory_name $map_name 0.1 10 4
  ./bin/knn_dataset_generator $directory_name $map_name 0.1 10 16

done