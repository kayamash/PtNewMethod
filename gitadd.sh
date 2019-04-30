#!/bin/bash
branch="master"
file1="run.cpp"
file2="PtNewMethod.cpp"
file3="PtNewMethod.chh"
file4="kayamashForLUT.cpp"
file5="kayamashForLUT.chh"
file6="gitadd.sh"
add="git add "
message="Add New Method"
push="git push origin "

eval $add$file1
eval $add$file2
eval $add$file3
eval $add$file4
eval $add$file5
eval $add$file6
#git add -A
git commit -m "${message}"
eval $push$branch




