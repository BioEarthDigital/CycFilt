target/release/cyc_filt -i example/sample.fastq.gz -o example/sample.clean.fastq.gz -q 15 -l 3000 -c 2 -a example/adapter.fa -s
python3 scripts/CycFqStatPlot.py --csv example/sample.clean.fastq.stat.csv.gz --min-quality 7 --min-length 3000 --out-dir example
