#!/bin/bash
#
# --- configure a bggen run

 PROG=`basename $0`
# echo Start $PROG

 ffp=particles.ffr  
 ffw=run_mcwrapper.ffr

# ==========   Help routine start
 help_dis () {
  cat <<help_doc
 $PROG - a script to configure the bggen input files.
  Call parameters:
 $PROG main_control_file_name.ffr
 it links this file and 2 static data files to FFREAD-readable names

 Example: $PROG run_pyth.ffr
help_doc
 }
# ==========   Help routine end

  if [ "$#" -lt 1 ]; then
    echo "***" Error - too few parameters. Usage:
    help_dis
    exit 2
  fi
 
  if [ "$1" = "?" -o "$1" = "-?" ]; then
    help_dis
    exit 1
  fi

  ffm=$1
  
  ok=1
  if [ ! -f $ffm ]; then
    echo "***" No such file: $ffm
    ok=0
  fi
  if [ ! -f $ffp ]; then
    echo "***" No such file: $ffp
    ok=0
  fi
  if [ ! -f $ffw ]; then
    echo "***" No such file: $ffw
    ok=0
  fi
  if [ $ok -eq 0 ]; then
    help_dis
    exit 2
  fi

  ok=1
  grep '^READ\s*16' $ffm > /dev/null
  if [ $? -ne 0 ]; then
    echo "***" No command READ 16 at the beginning of $ffm - should read the particle list
    ok=0
  fi 
  grep '^READ\s*17' $ffm > /dev/null
  if [ $? -ne 0 ]; then
    echo "***" No command READ 17 at the end of $ffm - should read the MCwrapper-changeable variables 
    ok=0
  fi 
  if [ $ok -eq 0 ]; then
    help_dis
    exit 2
  fi

  ln -s -f $ffm fort.15
  ln -s -f $ffp fort.16
  ln -s -f $ffw fort.17

