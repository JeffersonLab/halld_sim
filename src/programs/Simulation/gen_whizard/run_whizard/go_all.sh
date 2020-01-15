wf=gluex_primex_qed_bkgs
out_dir=/volatile/halld/home/ijaegle/qed_whizard
swif create $wf 
for fs in ae eee aeee; do
    ./run_whizard 1e7 whizard_$fs.cfg 61321 $wf $out_dir 
done
swif run $wf
