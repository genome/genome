package Genome::Model::Tools::Vcf::MergeWithAnnotation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeWithAnnotation.pm -    Merge variants in a VCF file with annotation and dbSNP information
#               
#   AUTHOR:      Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#   CREATED:    03/13/2012 by D.K.
#   MODIFIED:   03/13/2012 by D.K.
#
#   NOTES:   
#         
#####################################################################################################################################

my %vep_class_rank = ();
$vep_class_rank{'-'} = 				0;
$vep_class_rank{'NMD_TRANSCRIPT'} = 		0;
$vep_class_rank{'PARTIAL_CODON'} = 		0;
$vep_class_rank{'INTERGENIC'} = 		0;
$vep_class_rank{'UPSTREAM'} = 			1;
$vep_class_rank{'DOWNSTREAM'} = 		2;
$vep_class_rank{'INTRONIC'} = 			3;
$vep_class_rank{'5PRIME_UTR'} = 		4;
$vep_class_rank{'3PRIME_UTR'} = 		5;
$vep_class_rank{'WITHIN_NON_CODING_GENE'} = 	6;
$vep_class_rank{'WITHIN_MATURE_miRNA'} = 	7;
$vep_class_rank{'CODING_UNKNOWN'} = 	        7;
$vep_class_rank{'SYNONYMOUS_CODING'} = 		8;
$vep_class_rank{'STOP_LOST'} = 			9;
$vep_class_rank{'SPLICE_SITE'} = 		10;
$vep_class_rank{'ESSENTIAL_SPLICE_SITE'} = 	11;
$vep_class_rank{'NON_SYNONYMOUS_CODING'} = 	12;
$vep_class_rank{'STOP_GAINED'} = 		13;

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Vcf::MergeWithAnnotation {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
    vcf_file   => { is => 'Text', doc => "VCF file in uncompressed format", is_input => 1},
    transcript_annotation_file   => { is => 'Text', doc => "WU Transcript annotation in its native format", is_input => 1},
    vep_annotation_file   => { is => 'Text', doc => "VEP annotation in its native format", is_input => 1},
    regfeat_annotation_file   => { is => 'Text', doc => "Regulatory feature annotation in six-column TSV", is_optional => 1, is_input => 1},
    dbsnp_file => {is => 'String', doc => "dbSNP file ", is_optional => 1, is_input => 1},
    output_file   => { is => 'Text', doc => "Output file to receive limited SNPs", is_input => 1, is_output => 1},
    verbose   => { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merge variants in a VCF file with VEP and dbSNP annotation"                 
}

sub help_synopsis {
    return <<EOS
This command merges a VCF with VEP and dbSNP annotation
EXAMPLE:	gmt vcf merge-with-annotation --vcf myVCF.vcf --vep-annotation-file myVEP.output --dbsnp-file mySites.dbsnp --output-file merged.variants.tsv ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;
    
    ## Get required parameters ##
    my $vcf_file = $self->vcf_file;
    my $vep_file = $self->vep_annotation_file;
    my $output_file = $self->output_file;
    my $verbose = 1 if(defined($self->verbose));
    
    my %stats = ();
    $stats{'vcf_lines'} = $stats{'vcf_lines_pass'} = $stats{'vcf_lines_pass_snp'} = 0;
    
    ## Load the annotation ##
    warn "Loading Transcript annotation...\n";
    my %transcript_annotation = load_annotation($self->transcript_annotation_file);
    warn "Loading VEP annotation...\n";
    my %vep_annotation = load_vep($vep_file);
    warn "Loading dbSNP annotation...\n";
    my %dbsnp_annotation = load_dbsnp($self->dbsnp_file) if($self->dbsnp_file);

    warn "Loading RegFeat annotation...\n";
    my %regfeat_annotation = load_annotation($self->regfeat_annotation_file) if($self->regfeat_annotation_file);
    
    my %vep_class_counts = my %regfeat_class_counts = ();
    
    ## Open the outfile ##
    if($output_file)
    {
        open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
    }
    
    my $input = new FileHandle ($vcf_file);
    my $lineCounter = 0;

    my @sample_names = ();
    my $num_samples = 0;
    
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);
        my $numContents = @lineContents;
        
        if(substr($line, 0, 1) eq '#')
        {
            if(substr($line, 0, 6) eq '#CHROM')
            {
                ## parse out samples ##join("\t", $gene, $transcript, $trv_type, $c_position, $aa_change, $ucsc_cons, $domain);
                print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\t";
                print OUTFILE "dbSNP_rs\tdbSNP_submitters\tdbSNP_maf_alleles\tdbSNP_mafs\t";
                print OUTFILE "vartype\tgene\ttranscript\ttrv_type\tc_position\taa_change\tucsc_cons\tdomain\t";
                print OUTFILE "VEP_transcript\tVEP_gene\tVEP_class\tVEP_codon\tVEP_residue\tVEP_aa_change\tVEP_polyphen\tVEP_sift\tVEP_condel\t";
                print OUTFILE "RegFeat\t";
                print OUTFILE "SamplesRef\tSamplesHet\tSamplesHom\tSamplesMissing\t";
                for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
                {
                    $sample_names[$num_samples] = $lineContents[$colCounter];
                    print OUTFILE "$lineContents[$colCounter]\t";
                    $num_samples++;
                }
                
                print OUTFILE "\n";
                
                warn "$num_samples samples\n";
                $stats{'num_samples'} = $num_samples;
            }
        }
        else
        {
            my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $line);
            $stats{'vcf_lines'}++;
            
            ## Only take sites passing filters ##
            if($filter eq '.' || uc($filter) eq 'PASS')
            {
                $stats{'vcf_lines_pass'}++;
                
                ## Only take sites with a variant allele #
                if($ref && $var && $var ne '.')
                {
                    $stats{'vcf_lines_pass_snp'}++;
                    
                    my $key = join("\t", $chrom, $position, $position, $ref, $var);

                    ## Fix the key ##
                    if(!$vep_annotation{$key})
                    {
                        ## Get shortened key ##
                        
                        $key = join("\t", $chrom, $position, $position, $ref, substr($var, 0, 1));
                    }

                    

                    ## Get dbSNP info if known ##
                    my $dbsnp_info = join("\t", "-", "-", "-", "-");
                    if($dbsnp_annotation{$key})
                    {
                        $stats{'vcf_lines_pass_snp_dbsnp'}++;
#                        my ($rs, $loc, $type, $valstatus, $submitters, $blacklist, $maf_alleles, $mafs) = split(/\t/, $dbsnp_annotation{$key});
                        my ($rs, $submitters, $maf_alleles, $mafs) = split(/\t/, $dbsnp_annotation{$key});
                        $maf_alleles = "-" if(!$maf_alleles);
                        $mafs = "-" if(!$mafs);
#                        $dbsnp_info = join("\t", $rs, $valstatus, $maf_alleles, $mafs);
                        $dbsnp_info = join("\t", $rs, $submitters, $maf_alleles, $mafs);
                    }
                    
                    my $num_gt_ref = my $num_gt_het = my $num_gt_hom = my $num_gt_missing = 0;
                    my $genotypes = "";
                    for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
                    {
                        my ($genotype) = split(/\:/, $lineContents[$colCounter]);
                        $genotype = convert_genotype($ref, $var, $genotype);
                        
                        if($genotype eq $ref . $ref)
                        {
                            $num_gt_ref++;
                        }
                        elsif($genotype eq $ref . $var)
                        {
                            $num_gt_het++;
                        }
                        elsif($genotype eq $var . $var)
                        {
                            $num_gt_hom++;   
                        }
                        else
                        {
                            $num_gt_missing++;
                        }
                        $genotypes .= "\t" if($genotypes);
                        $genotypes .= $genotype;
                    }

                    ## Get VEP annotation info ##

                    my $tx_annot = "-\t-\t-\t-\t-\t-\t-";
                    $tx_annot = $transcript_annotation{$key} if($transcript_annotation{$key});

                    my $vep_annot = "-\t-\t-\t-\t-\t-\t-\t-\t-";

                    if($vep_annotation{$key})
                    {
                        $stats{'vcf_lines_pass_snp_vep'}++;
#                        print join("\t", $key, $dbsnp_info, $vep_annotation{$key}) . "\n";

                        my @vep_annotation = split(/\t/, $vep_annotation{$key});
                        my $vep_class = $vep_annotation[2];
                        $vep_class_counts{$vep_class}++;
                        
                        $vep_annot = $vep_annotation{$key};
                        
                    }
                    else
                    {
#                        warn "No VEP annotation for $key\n";
                    }

                    ## Get regfeats annotation ##
                    
                    my $regfeat_annot = "-";

                    if($regfeat_annotation{$key})
                    {
                        $regfeat_annot = $regfeat_annotation{$key};
                        
                        my @regfeats = split(/\,/, $regfeat_annot);
                        foreach my $feature (@regfeats)
                        {
                            $regfeat_class_counts{$feature}++;                            
                        }
                        
                        $stats{'vcf_lines_pass_snp_regfeat'}++;

                    }


                    print OUTFILE join("\t", $key, $dbsnp_info, $tx_annot, $vep_annot, $regfeat_annot, $num_gt_ref, $num_gt_het, $num_gt_hom, $num_gt_missing, $genotypes) . "\n" if($self->output_file);



                }
            }
        }
        
        if(0)   ## Use to stop as needed ##
        {
            print $stats{'vcf_lines'} . " lines parsed from VCF\n";
            print $stats{'vcf_lines_pass'} . " passed per-site filters\n";
            print $stats{'vcf_lines_pass_snp'} . " passed per-site filters and were variant\n";
            print $stats{'vcf_lines_pass_snp_vep'} . " had VEP annotation\n";
            print $stats{'vcf_lines_pass_snp_dbsnp'} . " were dbSNPs\n";
            return 1;
        }
    }
    
    close($input);

            print $stats{'vcf_lines'} . " lines parsed from VCF\n";
            print $stats{'vcf_lines_pass'} . " passed per-site filters\n";
            print $stats{'vcf_lines_pass_snp'} . " passed per-site filters and were variant\n";
            print $stats{'vcf_lines_pass_snp_vep'} . " had VEP annotation\n";
            print $stats{'vcf_lines_pass_snp_dbsnp'} . " were dbSNPs\n";

    
    close(OUTFILE) if($output_file);
    
    
    return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub convert_genotype
{
	my ($ref, $var, $genotype) = @_;
	
	return("NN") if($genotype eq '.' || length($genotype) < 2);
	
	return($ref . $ref) if($genotype eq '0/0');
	return($ref . $var) if($genotype eq '0/1');
	return($var . $var) if($genotype eq '1/1');
	
	if($var =~ '\,')
	{
		my @vars = split(/\,/, $var);
		
		my ($gt1, $gt2) = split(/\//, $genotype);
		
		if($gt1 == 0)
		{
			$genotype = $ref;
		}
		else
		{
			$genotype = $vars[$gt1 - 1];
		}

		if($gt2 == 0)
		{
			$genotype .= $ref;
		}
		else
		{
			$genotype .= $vars[$gt2 - 1];
		}
	}
	else
	{
		warn "Unable to convert $ref $var $genotype\n";		
	}

	return($genotype);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_vep
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
    my $input = new FileHandle ($annotation_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        
        if(substr($line, 0, 1) eq '#')
        {
            ## Ignore header line ##
        }
        else
        {
            my @lineContents = split(/\t/, $line);
            my ($string) = split(/\t/, $line);
            my ($chrom, $position, $ref, $var) = split(/[\_\/]/, $string);
            my $key = join("\t", $chrom, $position, $position, $ref, $var);

            my $ens_gene = $lineContents[3];
            my $class = $lineContents[6];
            my $cdna_pos = $lineContents[7];
            my $cds_pos = $lineContents[8];
            my $protein_pos = $lineContents[9];
            my $amino_acids = $lineContents[10];
            my $extra = $lineContents[13];

            ## Reset extra variables
            my $gene = my $polyphen = my $sift = my $condel = "";

            my @extraContents = split(/\;/, $extra);
            foreach my $entry (@extraContents)
            {
                    my ($key, $value) = split(/\=/, $entry);

                    $gene = $value if($key eq 'HGNC');
                    $polyphen = $value if($key eq 'PolyPhen');
                    $sift = $value if($key eq 'SIFT');
                    $condel = $value if($key eq 'Condel');	
            }

            my @classes = split(/\,/, $class);
            foreach my $class (@classes)
            {
                    if($polyphen || $sift || $condel)
                    {
                            $annotation{$key} .= "\n" if($annotation{$key});
                            $annotation{$key} .= join("\t", $ens_gene, $gene, $class, $cdna_pos, $protein_pos, $amino_acids, $polyphen, $sift, $condel)
#				print join("\t", $chrom, $position, $alleles, $ens_gene, $gene, $class, $protein_pos, $amino_acids, $polyphen, $sift, $condel) . "\n";
                    }
                    else
                    {
                            $annotation{$key} .= "\n" if($annotation{$key});
                            $annotation{$key} .= join("\t", $ens_gene, $gene, $class, $cdna_pos, "", "", "", "", "")				;
                    }				
            }


        }

    }
    
    close($input);
    
    
    ## Go through each key that has annotation and choose the top result ##
    
    foreach my $key (keys %annotation)
    {
        my @vepResults = split(/\n/, $annotation{$key});
        @vepResults = sort bySeverity @vepResults;
        my $top_result = $vepResults[0];
        $annotation{$key} = $top_result;
    }
    

    return(%annotation);
}



#############################################################
# load_vep_results - parses the file
#
#############################################################

sub bySeverity
{
	my ($ens_gene_a, $gene_a, $class_a, $cdna_pos_a, $protein_pos_a, $amino_acids_a, $polyphen_a, $sift_a, $condel_a) = split(/\t/, $a);

	$polyphen_a = 0 if(!$polyphen_a);
	$sift_a = 0 if(!$sift_a);
	$condel_a = 0 if(!$condel_a);
	
	if($polyphen_a)
	{
		my @temp = split(/[\(\)]/, $polyphen_a);
		$polyphen_a = $temp[1];
	}
	if($sift_a)
	{
		my @temp = split(/[\(\)]/, $sift_a);
		$sift_a = $temp[1];
	}
	if($condel_a)
	{
		my @temp = split(/[\(\)]/, $condel_a);
		$condel_a = $temp[1];
	}
	
	my ($ens_gene_b, $gene_b, $class_b, $cdna_pos_b, $protein_pos_b, $amino_acids_b, $polyphen_b, $sift_b, $condel_b) = split(/\t/, $b);

	$polyphen_b = 0 if(!$polyphen_b);
	$sift_b = 0 if(!$sift_b);
	$condel_b = 0 if(!$condel_b);

	if($polyphen_b)
	{
		my @temp = split(/[\(\)]/, $polyphen_b);
		$polyphen_b = $temp[1];
	}
	if($sift_b)
	{
		my @temp = split(/[\(\)]/, $sift_b);
		$sift_b = $temp[1];
	}
	if($condel_b)
	{
		my @temp = split(/[\(\)]/, $condel_b);
		$condel_b = $temp[1];
	}
	
	my $fxn_code_a = fxn_class_code($class_a);
	my $fxn_code_b = fxn_class_code($class_b);
        
        die "Got no code for $class_a\n" if(!$class_a);
        die "Got no code for $class_b\n" if(!$class_b);

	## Sort by function code severity first ##
	$fxn_code_b <=> $fxn_code_a
	or
	$polyphen_b <=> $polyphen_a
	or
	$sift_b <=> $sift_a
	or
	$condel_b <=> $condel_a
}


#############################################################
# load_vep_results - parses the file
#
#############################################################

sub fxn_class_code
{
	my $class = shift(@_);
	
	my @classes = split(/\,/, $class);
	my $num_classes = @classes;
	
	if($num_classes > 1)
	{
		@classes = sort byCode (@classes);
		$class = $classes[0];
	}
	
	foreach my $test_class (keys %vep_class_rank)
	{
		return($vep_class_rank{$test_class}) if($class eq $test_class);
	}
	
        die "No Rank provided for $class\n";
	return(0);
}


sub byCode
{
	my $code_a = my $code_b = 0;

	foreach my $test_class (keys %vep_class_rank)
	{
		$code_a = $vep_class_rank{$test_class} if($a eq $test_class);
		$code_b = $vep_class_rank{$test_class} if($b eq $test_class);
	}
	
	$code_b <=> $code_a;
}

#############################################################
# load_vep_results - parses the file
#
#############################################################

sub get_code
{
	my $class = shift(@_);
	
	my @classes = split(/\,/, $class);
	my $num_classes = @classes;
	
	if($num_classes > 1)
	{
				
	}
	
	foreach my $test_class (keys %vep_class_rank)
	{
		return($vep_class_rank{$test_class}) if($class eq $test_class);
	}
	
	return(0);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_dbsnp
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
    my $input = new FileHandle ($annotation_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        
        my ($chrom, $chr_start, $chr_stop, $ref, $var, $rs, $genomic, $single, $val_status, $submitters, $zero, $maf_alleles, $mafs) = split(/\t/, $line);

        my $rest_of_line = "";

        if($maf_alleles)
        {
            $rest_of_line = join("\t", $rs, $submitters, $maf_alleles, $mafs);            
        }
        else
        {
            $rest_of_line = join("\t", $rs, $submitters);            
        }

        my @vars = split(/\//, $var);
        foreach my $this_var (@vars)
        {
            my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $this_var);
            $annotation{$key} = $rest_of_line;
        }

    }
    
    close($input);

    return(%annotation);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_annotation
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
    my $input = new FileHandle ($annotation_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        
        my ($chrom, $chr_start, $chr_stop, $ref, $var, $gene, $transcript, $species, $source, $version, $strand, $status, $trv_type, $c_position, $aa_change, $ucsc_cons, $domain) = split(/\t/, $line);

#        my $rest_of_line = "";
 #       
  #      my @lineContents = split(/\t/, $line);
#        my $numContents = @lineContents;
#        for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
#        {
#            $rest_of_line .= "\t" if($colCounter > 5);
#            $rest_of_line .= $lineContents[$colCounter];
#        }

#        my @vars = split(/\//, $var);
#        foreach my $this_var (@vars)
#        {
            my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
            $gene = "-" if(!$gene);
            $transcript = "-" if(!$transcript);
            $c_position = "-" if(!$c_position);
            $aa_change = "-" if(!$aa_change);
            $domain = "-" if(!$domain);
            $ucsc_cons = "-" if(!$ucsc_cons);
            $annotation{$key} = join("\t", $gene, $transcript, $trv_type, $c_position, $aa_change, $ucsc_cons, $domain);
#        }

    }
    
    close($input);

    return(%annotation);
}



1;

