
package Genome::Model::Tools::MethylArray::ProbeDistribution;

use strict;
use Data::Dumper;
use Genome;

class Genome::Model::Tools::MethylArray::ProbeDistribution {
    is => ['Genome::Model::Tools::MethylArray::Base'],
    has_input => [
        probe_list_file => {
        	is => 'Text',
            doc => 'Input list of probes of interest',
        },
        annotation_tsv_file => {
            is => 'Text',
            doc => 'Input annotation TSV file from illumina',
        },
        output_file => {
            is => 'Text',
            is_output=> 1,
            doc => 'Raw Output file from Dumper containing counts of features',
        },
    ],
};

sub execute {
    my $self = shift;
    my $probe_list_file = $self->probe_list_file;
    my $illumina_annotation_file = $self->annotation_tsv_file;

	my %probe_hash;
    		
	my $probe_list_fh = Genome::Sys->open_file_for_reading($probe_list_file);		
		
	while (my $probe_row = $probe_list_fh->getline) 
		{
			chomp $probe_row;
			$probe_hash{$probe_row} = ();
		}

    my $illumina_fh = Genome::Sys->open_file_for_reading($illumina_annotation_file);
    my %feature_count;
    while (my $anno_line = $illumina_fh->getline) 
    {

		my @annotation = split('\t',$anno_line);
		my $probe_name = $annotation[0];
		my @ucsc_refgene_groups = split(';',$annotation[23]);
		my $gene_group = $ucsc_refgene_groups[0];
	
    	if ( exists($probe_hash{$probe_name}) ) {
    		$feature_count{$gene_group}++; 
       		 #print $probe_name."\t".$gene_group."\t".$gene_group."\n";
    	}
	}
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_file);
    while (my ($key, $value) = each(%feature_count))
    {
    	if ($key eq "")
     	{
     		print $output_fh "Intergenic"."\t ".$value."\n";
    	 }
     	else
     	{
    		 print $output_fh $key."\t ".$value."\n";
     	}
    }
	#print $output_fh Data::Dumper::Dumper(%feature_count) ."\n";

	return 1;   
}

1;



__END__

1	IlmnID
2	Name
3	AddressA_ID
4	AlleleA_ProbeSeq
5	AddressB_ID
6	AlleleB_ProbeSeq
7	Infinium_Design_Type
8	Next_Base
9	Color_Channel
10	Forward_Sequence
11	Genome_Build
12	CHR
13	MAPINFO
14	SourceSeq
15	Chromosome_36
16	Coordinate_36
17	Strand
18	Probe_SNPs
19	Probe_SNPs_10
20	Random_Loci
21	Methyl27_Loci
22	UCSC_RefGene_Name
23	UCSC_RefGene_Accession
24	UCSC_RefGene_Group
25	UCSC_CpG_Islands_Name
26	Relation_to_UCSC_CpG_Island
27	Phantom
28	DMR
29	Enhancer
30	HMM_Island
31	Regulatory_Feature_Name
32	Regulatory_Feature_Group
33	DHS
