# tandem_repeats

Script to find tandem repeats in sequences/contigs/scaffolds. 

Requirements:

python3 (tested with 3.8.2)\
biopython (tested with 1.79)\
pandas (tested with 1.4.1)\
smart_open (tested with 5.2.1)\
MUMmer (tested with 3.23)

Example usage with 10 processes:

`python3 tandem_repeats.py --force -i input.fa --min_region_len 50 --min_repeat_count 3 -p 10 -o input.fa.tandem-repeats`
