## Simulation

usage: src/simulation.py -m n_mutations -g n_genomes -n n_cells -c coverage [-t threshold] -o output_prefix

optional arguments:
  -h, --help  		show this help message and exit
  
  -m, --mutations     	number of mutations 
  
  -n, --cells        	number of cells
  
  -g, --genomes		number of genomes / clones
  
  -c, --coverage 	expected sequencing coverage
  
  -t, --threshold	minimum variant allele frequency
  
  -o, --output        	output prefix

synthetic datasets can be found in directory /ground_truth
