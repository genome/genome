package Genome::Model::Tools::EpitopePrediction::ParseNetmhcpanKey;

use strict;
use warnings;
use Data::Dumper;
use Genome;

class Genome::Model::Tools::EpitopePrediction::ParseNetmhcpanKey {
    is => ['Genome::Model::Tools::EpitopePrediction::Base'],
    has_input => [
        netmhc_file => {
        	is => 'Text',
            doc => 'Raw output file from Netmhc',
        },
        parsed_file => {
        	is => 'Text',
            doc => 'File to write the parsed output',
        },
        output_type => {
            is => 'Text',
            doc => 'Type of epitopes to report in the final output - select \'top\' to report the top epitopes in terms of fold changes,  \'all\' to report all predictions ',
            valid_values => ['top','all'],
        },
        key_file => {
            is => 'Text',
            doc => 'Key file for lookup of FASTA IDs'
        },
        
    ],
};

sub help_brief {
    "Runs NetMHC for MHC Class I epitope prediction

",
}



sub execute {
    my $self = shift;
	my %netmhc_results;
	my %epitope_seq;
	my %position_score;
	
	my $type = $self->output_type;
	my $input_fh = Genome::Sys->open_file_for_reading($self->netmhc_file);
	my $output_fh = Genome::Sys->open_file_for_writing($self->parsed_file);
	my $key_fh	= Genome::Sys->open_file_for_reading($self->key_file);
	
	
	### HASH FOR KEY FILE LOOKUP#####
	
	my %list_hash;
	while (my $keyline = $key_fh->getline) 
	{
		chomp $keyline;
		
		# >WT_8	>WT.GSTP1.R187W	   
	#	$list_hash{$keyline}  = ();	
		my ($pan_name, $original_name)    = split ( /\t/, $keyline );	
    		$list_hash{$pan_name} =();
		$list_hash{$pan_name}{'name'}   = $original_name;
		
	}
		#print Dumper %list_hash;
	#########
	
	
	while (my $line = $input_fh->getline) {
	chomp $line;
    if ($line =~ /^\s+\d+/) 
    {
	#	print $line."\n";
		my @result_arr = split (/\s+/,$line);
    	#print Dumper(@result_arr);
		# <space>  6  HLA-A*03:01  RLSAWPKLK            MT_8         0.656        41.45     0.12 <= SB
    	    	
    	my $position = $result_arr[1];
    	my $score = $result_arr[6];
    	my $epitope = $result_arr [3];
    	my $protein_pan_name = $result_arr[4];
		
		if ( exists( $list_hash{$protein_pan_name} ) )
		{
			#print $_."\n";
			my $protein = $list_hash{$protein_pan_name}{'name'};
#>MT.FBN3.P1547L
	    		my @protein_arr = split (/\./,$protein);
    	
    	#print Dumper(@protein_arr);
    	
    	my $protein_type = $protein_arr[0];
    	my $protein_name = $protein_arr[1];
    	my $variant_aa =  $protein_arr[2];
       
        $netmhc_results{$protein_type}{$protein_name}{$variant_aa}{$position} = $score;
        $epitope_seq{$protein_type}{$protein_name}{$variant_aa}{$position} = $epitope;
        
   	 }
	}

	}

	print $output_fh join("\t","Gene Name","Point Mutation","Sub-peptide Position","MT score", "WT score","MT epitope seq","WT epitope seq","Fold change")."\n";
	my $rnetmhc_results = \%netmhc_results;
	my $epitope_seq = \%epitope_seq;
	my @score_arr;
	for my $k1 ( sort keys %$rnetmhc_results ) {
	 my @positions;
	   if ($k1 eq 'MT')
	   {
      #  print "$k1\t";

        for my $k2 ( sort keys %{$rnetmhc_results->{ $k1 }} ) {
            #print "$k2\t";
     
            for my $k3 ( sort keys %{$rnetmhc_results->{ $k1 }->{ $k2 }} ) {
                #print "\t$k3";
				%position_score = %{$netmhc_results{$k1}{$k2}{$k3}};
				@positions = sort {$position_score{$a} <=> $position_score{$b}} keys %position_score;
				my $total_positions = scalar @positions; 
				
				if ($type eq 'all')
				{
					
				
					for (my $i = 0; $i < $total_positions; $i++){
					
				
						print $output_fh join("\t",$k2,$k3,$positions[$i],$position_score{$positions[$i]})."\t";
						print $output_fh $rnetmhc_results->{ 'WT'}->{ $k2 }->{ $k3 }->{$positions[$i]}."\t";
				
						print $output_fh $epitope_seq->{'MT'}->{$k2}->{$k3}->{$positions[$i]}."\t";
						print $output_fh $epitope_seq->{'WT'}->{$k2}->{$k3}->{$positions[$i]}."\t";
						my $fold_change = $rnetmhc_results->{ 'WT'}->{ $k2 }->{ $k3 }->{$positions[$i]}/$position_score{$positions[$i]};
						my $rounded_FC = sprintf("%.3f", $fold_change);
						print $output_fh $rounded_FC."\n";	
					}
				}
				if ($type eq 'top')
				{
					
					print $output_fh join("\t",$k2,$k3,$positions[0],$position_score{$positions[0]})."\t";
					print $output_fh $rnetmhc_results->{ 'WT'}->{ $k2 }->{ $k3 }->{$positions[0]}."\t";
				
					print $output_fh $epitope_seq->{'MT'}->{$k2}->{$k3}->{$positions[0]}."\t";
					print $output_fh $epitope_seq->{'WT'}->{$k2}->{$k3}->{$positions[0]}."\t";
					my $fold_change = $rnetmhc_results->{ 'WT'}->{ $k2 }->{ $k3 }->{$positions[0]}/$position_score{$positions[0]};
					my $rounded_FC = sprintf("%.3f", $fold_change);
					print $output_fh $rounded_FC."\n";	
				}
		
			}
				
           }
          }
        }

	
    
    return 1;
}



1;

__END__
### PREVIOUS FILE###


Pos	Peptide	ID	HLA-A03:01	Ave
1	EHLVITLGS	WT_NRIP3_p_R179	0.0182	0.0182
2	HLVITLGSL	WT_NRIP3_p_R179	0.0680	0.0680
3	LVITLGSLR	WT_NRIP3_p_R179	0.2807	0.2807
4	VITLGSLRL	WT_NRIP3_p_R179	0.0810	0.0810
5	ITLGSLRLD	WT_NRIP3_p_R179	0.0475	0.0475
6	TLGSLRLDC	WT_NRIP3_p_R179	0.0391	0.0391
7	LGSLRLDCP	WT_NRIP3_p_R179	0.0074	0.0074
8	GSLRLDCPA	WT_NRIP3_p_R179	0.0402	0.0402
9	SLRLDCPAA	WT_NRIP3_p_R179	0.0794	0.0794

####### KEY PAN FILE###

# NetMHCpan version 2.4

# Input is in FSA format

# Peptide length 9

HLA-A03:01 : Estimated prediction accuracy  0.853 (using nearest neighbor HLA-A03:01)

# Threshold for Strong binding peptides  50.000
# Threshold for Weak binding peptides 500.000
-----------------------------------------------------------------------------------
  pos          HLA    peptide        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
-----------------------------------------------------------------------------------
    0  HLA-A*03:01  LSAYVGRLS            WT_8         0.063     25369.40    32.00
    1  HLA-A*03:01  SAYVGRLSA            WT_8         0.155      9307.33     8.00
    2  HLA-A*03:01  AYVGRLSAR            WT_8         0.134     11697.07    10.00
    3  HLA-A*03:01  YVGRLSARP            WT_8         0.025     37976.32    50.00
    4  HLA-A*03:01  VGRLSARPK            WT_8         0.182      7008.24     7.00
    5  HLA-A*03:01  GRLSARPKL            WT_8         0.026     37897.49    50.00
    6  HLA-A*03:01  RLSARPKLK            WT_8         0.667        36.68     0.12 <= SB
    7  HLA-A*03:01  LSARPKLKA            WT_8         0.077     21675.38    32.00
    8  HLA-A*03:01  SARPKLKAF            WT_8         0.048     29628.37    50.00
    9  HLA-A*03:01  ARPKLKAFL            WT_8         0.028     37050.94    50.00
   10  HLA-A*03:01  RPKLKAFLA            WT_8         0.031     35648.11    50.00
   11  HLA-A*03:01  PKLKAFLAS            WT_8         0.007     46351.27    50.00
   12  HLA-A*03:01  KLKAFLASP            WT_8         0.114     14623.44    15.00
-----------------------------------------------------------------------------------

Protein WT_8. Allele HLA-A*03:01. Number of high binders 1. Number of weak binders 0. Number of peptides 13

-----------------------------------------------------------------------------------
# Threshold for Strong binding peptides  50.000
# Threshold for Weak binding peptides 500.000
-----------------------------------------------------------------------------------
  pos          HLA    peptide        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
-----------------------------------------------------------------------------------
    0  HLA-A*03:01  LSAYVGRLS            MT_8         0.063     25369.40    32.00
    1  HLA-A*03:01  SAYVGRLSA            MT_8         0.155      9307.33     8.00
    2  HLA-A*03:01  AYVGRLSAW            MT_8         0.037     33539.71    50.00
    3  HLA-A*03:01  YVGRLSAWP            MT_8         0.025     38080.70    50.00
    4  HLA-A*03:01  VGRLSAWPK            MT_8         0.284      2309.18     4.00
    5  HLA-A*03:01  GRLSAWPKL            MT_8         0.029     36637.15    50.00
    6  HLA-A*03:01  RLSAWPKLK            MT_8         0.656        41.45     0.12 <= SB
    7  HLA-A*03:01  LSAWPKLKA            MT_8         0.093     18267.70    16.00
    8  HLA-A*03:01  SAWPKLKAF            MT_8         0.051     28906.50    50.00
    9  HLA-A*03:01  AWPKLKAFL            MT_8         0.041     32064.48    50.00
   10  HLA-A*03:01  WPKLKAFLA            MT_8         0.018     41299.86    50.00
   11  HLA-A*03:01  PKLKAFLAS            MT_8         0.007     46351.27    50.00
   12  HLA-A*03:01  KLKAFLASP            MT_8         0.114     14623.44    15.00
-----------------------------------------------------------------------------------






