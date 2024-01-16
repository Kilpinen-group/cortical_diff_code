
#!/bin/bash


## Generate ped file with the PLINK plugin from GenomeStudio Software (ReportWizard)

## Use Plink1.9 to generate the vcf (not matched ref alleles)
/soft/plink-1.9/plink --map KIL_JUN_2021.map --ped KIL_JUN_2021_Plus.ped --recode vcf

## Curate vcf and match ref alleles, remove homozygous positions
/soft/Rscript vcfCuration.R

## Vcf needs to be recoded to the "chrX" format from NCBI

awk '{ 
        if($0 !~ /^#/) 
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0 
      }' curatedGeno.vcf > curatedGeno.recoded.vcf




