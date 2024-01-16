
#!/bin/bash

bsub -R'select[mem>150000] rusage[mem=150000]' -M150000 -J "mapRef[1-8]" -q long -o log/mapRef.%J.%I.out -e log/mapRef.%J.%I.err ./testReductions.sh testReductions.txt