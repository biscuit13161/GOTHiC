#!/bin/bash

set -x

for i in *.cpp
do
  cat ../copyright.txt $i > $i.new
  mv $i.new $i
done
