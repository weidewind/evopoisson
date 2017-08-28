foldername="fakes_h1_krya_skipstoppers_all"
for i in {0..9}
do
                for number in {1..10}
                do
                        it=$((number+10*i))
                                out="output/${foldername}/${it}_fake/nsyn/maxpath_not_subtracted/skip_stoppers"
                                perl grep_newsim_stats.pl -i $out &

                done
                wait
done

perl grep_fake_fdr.pl --skip_stoppers -i output/${foldername}
perl  metagrep_newsim_stats.pl --skip_stoppers -i output/${foldername}

