use strict;
use warnings;



my @samples = ('Basal-M_BY-10-2225L',
'Basal-M_BY-11-T1',
 'Basal-M_BY-12-2153L',
 'Basal-M_BY-13-2224L',
 'ClaudinLow-M_BY-7-2247R',
 'ClaudinLow-M_BY-8-T11',
 'ClaudinLow-M_BY-9-2151R',
 'Her2-M_BY-6-2304L2',
 'Luminal-M_BY-1-2250L',
 'Luminal-M_BY-2-2243L',
 'Luminal-M_BY-3-2208L',
 'p53null-M_BY-4-2151L');
 
 
 my @alleles = ('H-2_Ld', 'H-2_Kd','H-2_Dd');
 my @lengths = ('8','9','10');
 
 
 
 
foreach my $sample(@samples)
{
	my $command1 = join (
                        q/ /,
                        'bsub -q techd',
                        '\'gmt epitope-prediction get-wildtype --input-tsv-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/snvs.'
                        .$sample.
                        '.cluster-filter.tier1.bed.transcript-annotation --output-tsv-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.
                        $sample.'/'.$sample.'_WT.tsv --anno-db=NCBI-mouse.combined-annotation --version=58_37k_v2\'');

#print $command1 . "\n";

my $command2 = join (
                        q/ /,
                        'bsub -q techd',
                        '\'gmt epitope-prediction generate-variant-seq21mer --input-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.
                        $sample.'/'.$sample.'_WT.tsv --output-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.
                        $sample.'/'.$sample.'_21mer\'');

#print $command2 . "\n";


my $command3 = join (
                        q/ /,
                        'perl /gscmnt/sata141/techd/jhundal/Schreiber_mouse/TESTING_PIPELINE/remove_stars.pl /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.
                        $sample.'/'.$sample.'_21mer /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.
                        $sample.'/'.$sample.'_21mer.fa');

#print $command3 . "\n";


	foreach my $allele(@alleles)
#H-2_Ld,H-2_Kd,H-2_Dd


			{
				
				foreach my $length(@lengths)
#8,9,10
					{
						
						my $command4 = join (
                        q/ /,
                        'bsub -q techd',
                        'gmt epitope-prediction run-netmhc --allele='.$allele,
                        '--epitope-length='.$length,
                        '--fasta-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_21mer.fa',
                        '--output-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_netmhc.'.$allele.'.'.$length);
						
						#print $command4."\n";
						
						
					#	gmt epitope-prediction parse-netmhc-mutant --netmhc-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Kd/Basal-M_BY-11-T1.netmhc.Kd.10 --parsed-file=test --output-type=all
						
						my $command5 = join (
                        q/ /,
                        
                        'gmt epitope-prediction parse-netmhc-mutant',

                        '--netmhc-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_netmhc.'.$allele.'.'.$length,
                        '--parsed-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_netmhc.'.$allele.'.'.$length.'.parsed.all',
                        '--output-type=all');

						
					#	print $command5."\n";
						
						my $command6 = join (
                        q/ /,
                        
                        'gmt epitope-prediction parse-netmhc',

                        '--netmhc-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_netmhc.'.$allele.'.'.$length,
                        '--parsed-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/'.$sample.'/'.$sample.'_netmhc.'.$allele.'.'.$length.'.parsed.all.nofilter',
                        '--output-type=all');

						
						print $command6."\n";
						
						
						
						

					}
	
			}


#system( "$command" );

} 
 





__END__





gmt epitope-prediction get-wildtype --input-tsv-file=/gscmnt/sata809/info/medseq/TP53-Knockout/analysis/February2012/snvs/snv_files/snvs.Basal-M_BY-11-T1.cluster-filter.tier1.bed.transcript-annotation --output-tsv-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_WT.tsv --anno-db=NCBI-mouse.combined-annotation --version=58_37k_v2


EpitopePrediction$ gmt epitope-prediction generate-variant-seq21mer --input-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_WT.tsv --output-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer.fa

EpitopePrediction$ perl /gscmnt/sata141/techd/jhundal/Schreiber_mouse/TESTING_PIPELINE/remove_stars.pl /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer.fa /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer_ns.fa


Tools$ bsub -q techd -u jhundal@genome.wustl.edu 'gmt epitope-prediction run-netmhc --allele=H-2_Ld --epitope-length=8 --fasta-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer_ns.fa --output-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1.netmhc.Ld.8'
Job <7209323> is submitted to queue <techd>.
head: cannot open `o' for reading: No such file or directory
Basal-M_BY-11-T1$ more notes 
gmt epitope-prediction get-wildtype --input-tsv-file=/gscmnt/sata809/info/medseq/TP53-Knockout/analysis/February2012/snvs/snv_files/snvs.Basal-M_BY-11-T1.cluster-filter.tier1.b
ed.transcript-annotation --output-tsv-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_WT.tsv --anno-db=NCBI-mouse.combined-annotation --
version=58_37k_v2


EpitopePrediction$ gmt epitope-prediction generate-variant-seq21mer --input-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_WT.tsv --out
put-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer.fa

EpitopePrediction$ perl /gscmnt/sata141/techd/jhundal/Schreiber_mouse/TESTING_PIPELINE/remove_stars.pl /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M
_BY-11-T1_21mer.fa /gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer_ns.fa


Tools$ bsub -q techd -u jhundal@genome.wustl.edu 'gmt epitope-prediction run-netmhc --allele=H-2_Ld --epitope-length=8 --fasta-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_
NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1_21mer_ns.fa --output-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1.netmhc.Ld.8'
Job <7209323> is submitted to queue <techd>.

Tools$ gmt epitope-prediction parse-netmhc --netmhc-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1.netmhc.Ld.8 --output-type=top --pars
ed-file=/gscmnt/sata141/techd/jhundal/KOMEN_BRC/p53_NULL/Basal-M_BY-11-T1/Basal-M_BY-11-T1.netmhc.Ld.8.parsed

