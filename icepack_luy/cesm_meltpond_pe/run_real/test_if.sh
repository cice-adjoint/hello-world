#!/bin/sh -f

iflag=`grep "IFLAG" icepack.runlog.191215-114747 | cut -f 2 -d"=" | sed 's/ //g'`
    echo "iflag= $iflag"
    if [[ $iflag -eq 0 ]]
    then
      echo 'iflag=0'
      echo 1 >> iflag.txt
     sed -n "` grep -n "FINAL POINT X=" icepack.runlog.191215-114747 | awk -F ":" '{print $1}' `,+1p" icepack.runlog.191215-114747 | tail -n 1 >> X.txt
    else
      echo 'iflag=-1'
      echo -1 >> iflag.txt
      echo NAN >> X.txt
    fi

