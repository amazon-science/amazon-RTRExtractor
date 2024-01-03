END=466

for ((i=1; i<=END; i++)); do
    python3 sanitize.py ../../Data/bug_bash/95_metis3agreement/ $i.txt -m
done
