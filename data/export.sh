array=()
for out in $(find . -name 'observed_shifts_*.txt')
do
    test=$(echo $out| cut  -d'_' -f 3 | cut -d'.' -f 1)
    printf "%s\n" $test
    lmd="larmord_$test.txt"
    printf "%s\n" $lmd
done
