package Genome::Model::Tools::Annotate::ReviseMaf;


# By the way here is the MAF dbSNP_Val_Status equivalent of 1-1-1-1-1-1:
#  by2Hit2Allele;byCluster;byFrequency;byHapMap;byOtherPop
#
# You think you could modify this to update column 14 and 15 of external
# maf files?
#
# If the result is a 1 report the method, if zero report nothing,
# semicolon separated.
#
# Better yet would be if you could modify the tool and add an option for
# this to run on any delimited file, with or without headers, and update
# any two columns of the users choice.


use strict;
use warnings;
use Genome;
use Genome::Info::IUB;


class Genome::Model::Tools::Annotate::ReviseMaf {
    is  => 'Genome::Model::Tools::Annotate',
    has => [
        maf_file => {
            type     => 'Text',
            is_input => 1,
            doc      => "MAF file."
        },
	output_file => {
	    type      => 'Text',
	    default   => 0,
	    is_optional => 1,
	    doc       => "Provide a name for your revised maf."
	},
	revise_Entrez_Gene_Id=> {
	    is => 'Boolean',
	    is_optional => 1,
	    default => 0,
	    doc => 'Will append/overwrite the locus_link_id from Genome::Site::TGI::Gene in the Entrez_Gene_Id column.'
	},
	order => {
	    type  =>  'String',
	    is_optional  => 1,
	    default => '1,2,5,6,7,11,12,13,14,15,16',
	    doc   =>  "provide the column number from 1 to 1+n for (Hugo_Symbol Entrez_Gene_Id Chromosome Start_position End_position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode) repectively."
	},
	delimiter => {
	    type  =>  'String',
	    is_optional  => 1,
	    default => "\t",
	    doc   =>  "provide the column seperator if other than a tab"
	},
	
	],
};

sub help_synopsis { 
    "gmt annotate revised-maf --maf-file maf --output-file revised_maf -revise-dbsnp -revise-Entrez-Gene-Id"
}

sub help_detail {
    return <<EOS
	Reads in a maf file and optionally revises the Entrez Gene Id and or dbsnp info. If opting for dbsnp the submitter and population allele frequencies will be appended to the last two columns 

	The defaults are set for a tab delimited with these fields in this order 

Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_position End_position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_file Sequencer

	To run on some variation of this format, see the -order and -delimiter options
	
EOS
}

sub execute { 
    
    my $self = shift;
    
    my ($maf,$dbsnp,$entrez_gene_ids,$column_n) = &parse_maf($self);
    my $dbsnp_submitters;
    if ($self->output_file) {
	my $output_file = $self->output_file;
	open(OUT,">$output_file")  || $self->error_message( "\n\nCouldn't open the output file $output_file\n\n.") && return;
    }
    
    my $order = $self->order;
    my $delimiter = $self->delimiter;
    my ($Hugo_Symbol_n,$Entrez_Gene_Id_n,$Chromosome_n,$Start_position_n,$End_position_n,$Reference_Allele_n,$Tumor_Seq_Allele1_n,$Tumor_Seq_Allele2_n,$dbSNP_RS_n,$dbSNP_Val_Status_n,$Tumor_Sample_Barcode_n) = split(/\,/,$order);

    foreach my $line_number (sort {$a<=>$b} keys %{$maf}) {
	my $line = $maf->{$line_number};
	
	my ($Hugo_Symbol,$Entrez_Gene_Id,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$Tumor_Seq_Allele2,$dbSNP_RS,$dbSNP_Val_Status,$Tumor_Sample_Barcode) = (split(/$delimiter/,$line))[$Hugo_Symbol_n -1 , $Entrez_Gene_Id_n -1 , $Chromosome_n -1 , $Start_position_n -1 , $End_position_n -1 , $Reference_Allele_n -1 , $Tumor_Seq_Allele1_n -1 , $Tumor_Seq_Allele2_n -1 , $dbSNP_RS_n -1 , $dbSNP_Val_Status_n -1 , $Tumor_Sample_Barcode_n -1];
	
	my $revise_dbSNP_RS;
	my $revise_dbSNP_Val_Status;
	my $revise_dbSNP_Submitters;
	my $revise_dbSNP_Allele_freq;
	if ($Chromosome eq "23") {$Chromosome = "X";}
	
	my $revise_Entrez_Gene_Id;
	if ($self->revise_Entrez_Gene_Id) {
	    my $entrez_gene_id = $entrez_gene_ids->{$Hugo_Symbol};
	    if ($entrez_gene_id) {
		$revise_Entrez_Gene_Id = $entrez_gene_id;
	    } else {
		unless ($Entrez_Gene_Id eq "Entrez_Gene_Id") {
		    $revise_Entrez_Gene_Id = '';
		}
	    }
	} else {
	    $revise_Entrez_Gene_Id = $Entrez_Gene_Id;
	}
	if ($self->revise_dbsnp) {
	    my $dbsnp_result = $dbsnp->{$Chromosome}->{$Start_position}->{$End_position}->{$Tumor_Sample_Barcode};
	    if ($dbsnp_result) {
		if ($dbsnp_result eq "no_hit") {
		    $revise_dbSNP_RS = '';
		    $revise_dbSNP_Val_Status = '';
		    $revise_dbSNP_Submitters = '';
		    $revise_dbSNP_Allele_freq = '';
		} else {
		    my ($rs_id,$submitters,$dbsnp_alleles,$validation,$freq) = split(/\,/,$dbsnp_result);
		    
		    if ($rs_id && $submitters && $dbsnp_alleles && $validation) {
			
			$validation =~ s/\;not_validated//gi;
			$validation =~ s/not_validated//;
			$validation =~ s/^\;//;
			$validation =~ s/\;$//;
			
			$revise_dbSNP_Val_Status = $validation;
			$revise_dbSNP_Submitters = $submitters;
			$revise_dbSNP_RS = $rs_id;
			if ($freq) {
			    $revise_dbSNP_Allele_freq = $freq;
			} else {
			    $revise_dbSNP_Allele_freq = '';
			}
		    } else {
			$revise_dbSNP_RS = '';
			$revise_dbSNP_Val_Status = '';
			$revise_dbSNP_Submitters = '';
			$revise_dbSNP_Allele_freq = '';
		    }
		} 
	    } else {
		if ($Hugo_Symbol eq "Hugo_Symbol") {
		#if ($dbSNP_RS eq "dbSNP_RS") {
		    $revise_dbSNP_RS = "dbSNP_RS";
		    $revise_dbSNP_Val_Status = $dbSNP_Val_Status;
		    $revise_dbSNP_Submitters = "dbSNP_Submitters";
		    $revise_Entrez_Gene_Id = $Entrez_Gene_Id;
		    $revise_dbSNP_Allele_freq = 'dbSNP_Population_Allele_Frequencies';
		} else {
		    $revise_dbSNP_RS = '';
		    $revise_dbSNP_Val_Status = '';
		    $revise_dbSNP_Submitters = '';
		    $revise_dbSNP_Allele_freq = '';
		}
	    }
	} else {
	    $revise_dbSNP_RS = "$dbSNP_RS";
	    $revise_dbSNP_Val_Status = "$dbSNP_Val_Status";
	    $revise_dbSNP_Submitters = '';
	    $revise_dbSNP_Allele_freq = '';
	} 
	
	my @line = split(/$delimiter/,$line);
	my @new_line;
	
	for my $n (1..$column_n + 1) {
	    
	    if ($n == $Entrez_Gene_Id_n) {
		push(@new_line,$revise_Entrez_Gene_Id);
	    } elsif ($n == $dbSNP_RS_n) {
		push(@new_line,$revise_dbSNP_RS);
	    } elsif ($n == $dbSNP_Val_Status_n) {
		push(@new_line,$revise_dbSNP_Val_Status);
	    } elsif ($n == $column_n + 1) {
		push(@new_line,$revise_dbSNP_Submitters);
		push(@new_line,$revise_dbSNP_Allele_freq);
	    } else {
		push(@new_line,$line[$n - 1]);
	    }
	}
	my $new_line;
	
	if (@new_line) {
	    no warnings 'uninitialized';
	    $new_line = join "$delimiter" , @new_line;
	}
	
	if ($self->output_file) {
	    print OUT qq($new_line\n);
	} else {
	    print qq($new_line\n);
	}
    }
}


sub parse_maf {
    my ($maf,$dbsnp,$entrez_gene_ids);
    my ($self) = @_;
    
    my $maf_file = $self->maf_file;
    unless (-f $maf_file) {
	$self->error_message( "\n\nCouldn't open the MAF $maf_file\n\n.");
	return;
    }
    
    if ($self->revise_dbsnp) {
	open(TEMP,">temp_out_file_for_look_up_variants.txt") || $self->error_message( "\n\nCouldn't open a file to write a lookup variants file\n\n.") && return;
    }
    
    my $delimiter = $self->delimiter;
    my $order = $self->order;
    my ($Hugo_Symbol_n,$Entrez_Gene_Id_n,$Chromosome_n,$Start_position_n,$End_position_n,$Reference_Allele_n,$Tumor_Seq_Allele1_n,$Tumor_Seq_Allele2_n,$dbSNP_RS_n,$dbSNP_Val_Status_n,$Tumor_Sample_Barcode_n) = split(/\,/,$order);
    
    my $n=0;
    my $column_n=0;

    open(MAF,$maf_file) || $self->error_message( "\n\nCouldn't open the MAF $maf_file\n\n.") && return;
    while (<MAF>) {
	chomp;
	my $line = $_;
	$n++;
	
	my ($Hugo_Symbol,$Entrez_Gene_Id,$Chromosome,$Start_position,$End_position,$Reference_Allele,$Tumor_Seq_Allele1,$Tumor_Seq_Allele2,$dbSNP_RS,$dbSNP_Val_Status,$Tumor_Sample_Barcode) = (split(/$delimiter/,$line))[$Hugo_Symbol_n -1 , $Entrez_Gene_Id_n -1 , $Chromosome_n -1 , $Start_position_n -1 , $End_position_n -1 , $Reference_Allele_n -1 , $Tumor_Seq_Allele1_n -1 , $Tumor_Seq_Allele2_n -1 , $dbSNP_RS_n -1 , $dbSNP_Val_Status_n -1 , $Tumor_Sample_Barcode_n -1];
		
	my @line = split(/$delimiter/,$line);
	my $line_n = @line;
	if ($line_n > $column_n) {
	    $column_n = $line_n;
	}
	
	$maf->{$n}=$line;
	
	if ($self->revise_Entrez_Gene_Id) {
	    my @gene_info = Genome::Site::TGI::Gene->get(gene_name => $Hugo_Symbol);
	    if (@gene_info) {
		for my $info (@gene_info) {
		    my $locus_link_id = $info->locus_link_id;
		    if ($locus_link_id) {
			my $id = $entrez_gene_ids->{$Hugo_Symbol};
			if ($id) {
			    unless ($id =~ /$locus_link_id/) {
				$entrez_gene_ids->{$Hugo_Symbol}="$id:$locus_link_id";
			    }
			} else {
			    $entrez_gene_ids->{$Hugo_Symbol}=$locus_link_id;
			}
		    }
		}
	    }
	}
	
	if ($self->revise_dbsnp) {
	    unless ($Hugo_Symbol eq "Hugo_Symbol") {
		my ($Tumor_Seq_uib);
		if ($Tumor_Seq_Allele1 =~ /^[ACGT]$/ && $Tumor_Seq_Allele2 =~ /^[ACGT]$/) {
		    ($Tumor_Seq_uib) = Genome::Info::IUB->iub_for_alleles($Tumor_Seq_Allele1,$Tumor_Seq_Allele2);
		} elsif ($Reference_Allele =~ /^[ACGT]$/) {
		    $Tumor_Seq_uib = $Reference_Allele;
		} else {
		    $Tumor_Seq_uib = "A";
		    $Reference_Allele = "A";
		}
		unless ($Reference_Allele =~ /^[ACGT]$/) {
		    $Tumor_Seq_uib = "A";
		    $Reference_Allele = "A";
		}
		
		if ($Chromosome eq "23") {$Chromosome = "X";}

		if ($Chromosome =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M)$/ && $Start_position =~ /^[\d]+$/ && $End_position =~ /^[\d]+$/) {
		    print TEMP qq($Chromosome\t$Start_position\t$End_position\t$Reference_Allele\t$Tumor_Seq_uib\t$Tumor_Sample_Barcode\n);
		} else {
		    print qq($Chromosome $Start_position $End_position are not valid coordinates. dbSNP results for this entry will be negative\n);
		}
	    }
	}
    }
    
    return ($maf,$dbsnp,$entrez_gene_ids,$column_n);
    
}

#Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_position End_position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_file Sequencer

