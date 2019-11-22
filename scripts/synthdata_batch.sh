stamp=100-triplets
workdir=/Users/rfabbri/cprg/vxlprg/lemsvpe/minus/scripts/results-synth/work-$stamp
mkdir $workdir 

# echoerr() { echo "LOG $@" 1>&2; }
echoerr() { printf "LOG %s\n" "$*" >&2; }

echoerr started sequential minus tester
set -x

while IFS= read -r sample_id || [ -n "$sample_id" ]
do
  sample_id_tr=`echo $sample_id |tr \  -`
  #mstdout=$workdir/${sample_id_tr}-$stamp-minus-synth-batch-stdout.$$;
  mstdout=/dev/null
  mstderr=$workdir/${sample_id_tr}-$stamp-minus-synth-batch-stderr.$$;
  minus-chicago-synth $sample_id 1>$mstdout 2>$mstderr
  solver_fail="$?"
  failfile=$workdir/${sample_id_tr}-$stamp-fail.$$;
  echo $solver_fail > $failfile
  tfile=$workdir/${sample_id_tr}-$stamp-time.$$;
  grep 'Time of solver' $mstderr |grep -o '[0-9][0-9]*ms'|grep -o '[0-9][0-9]*' > $tfile
done < 100-configurations-synthdata #tiny-configurations-synthdata 
