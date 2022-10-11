#!/bin/bash
# get the class name lists
ls -l | awk '{print $NF}' | awk '{sub(/.md/,"")sub(/class/,"")sub(/__/,"_")sub(/__/,"_");print "'\''"$1"'\','"}'
# get the class file list
ls -l | awk '{print "'\''"$NF"'\','"}'