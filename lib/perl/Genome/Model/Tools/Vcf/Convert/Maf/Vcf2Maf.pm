package Genome::Model::Tools::Vcf::Convert::Maf::Vcf2Maf;

#########################################################################
# Vcf2Maf - Converts a VCF file with annotations into a MAF file format #
# 								        #
# AUTHOR: Yanwen You (yyou@genome.wustl.edu)			        #
# 								        #
# CREATED: 6/28/2011 by Yanwen You				        #
# EDITED: 11/11/2011 by William Schierding
#########################################################################

use strict;
use warnings;
use Genome;
use Data::Dumper;

# Debugging
# use diagnostics;
# use Data::Dumper;

class Genome::Model::Tools::Vcf::Convert::Maf::Vcf2Maf {
    is => 'Command::V2',
    has_input => [ # All parameters are required
	    vcf_file => {
	        is => 'Text',
	        doc => 'VCF file to convert -- can only be single sample vcf',
            is_optional => 0,
	    },
	    annotation_file => {
	        is => 'Text',
	        doc => 'Annotation file of all positions within the vcf',
            is_optional => 0,
	    },
	    output_file => {
	        is => 'Text',
	        doc => 'Output MAF file',
            is_optional => 0,
	    },
        remove_silent => {
	        is => 'Boolean',
	        doc => 'Remove silent variants from the maf',
            default_value => 0,
	    },
        annotation_has_header => {
	        is => 'Boolean',
	        doc => 'In the pipeline, annotation files dont have a header',
            default_value => 0,
	    },
    ],
    # add has_optional_input here for optional arguments
};

sub help_brief {
    "Convert VCF file with annotation to a MAF file"
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    my $vcf_file = $self->vcf_file;
    my $annotation_file = $self->annotation_file;
    my $output_file = $self->output_file;

    # Verify existence of files
    unless(-e $vcf_file || $annotation_file) {
        $self->error_message("Error: VCF file or annotation does not exist!\n");
        die $self->error_message;
    }
    # Try to open the files
    unless(open VCF, "<$vcf_file") {
        $self->error_message("Error: Could not open file \"$vcf_file\"\n");
        die $self->error_message;
    }
    unless(open ANNOT, "<$annotation_file") {
        $self->error_message("Error: Could not open file \"$annotation_file\"\n");
        die $self->error_message;
    }
    unless(open MAF, ">$output_file") {
        $self->error_message("Error: Could not open file \"$output_file\" for output\n");
        die $self->error_message;
    }

    # Standard Maf column names
    my @maf_standard_columns = qw(Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_File Sequencer);
    my @maf_nonstandard_columns = qw(chromosome_name_WU start_WU stop_WU reference_WU variant_WU type_WU gene_name_WU transcript_name_WU transcript_species_WU transcript_source_WU transcript_version_WU strand_WU transcript_status_WU trv_type_WU c_position_WU amino_acid_change_WU ucsc_cons_WU domain_WU all_domains_WU deletion_substructures_WU transcript_error_WU);
    my @maf_columns = (@maf_standard_columns, @maf_nonstandard_columns);

    # Make the MAF header
    print MAF join("\t", @maf_columns), "\n";

#	my ($Hugo_Symbol, $Entrez_Gene_Id, $Gsc_Center, $NCBI_Build, $Chromosome, $Start_position, $End_position, $Strand, $Variant_Classification, $Variant_Type, $Reference_Allele, $Variant_Allele1, $Variant_Allele2, $dbSNP_RS, $dbSNP_Val_Status, $Sample_Barcode1, $Sample_Barcode2, $Match_Norm_Seq_Allele1, $Match_Norm_Seq_Allele2, $Validation_Allele1, $Validation_Allele2, $Match_Norm_Validation_Allele1, $Match_Norm_Validation_Allele2, $Verification_Status, $Validation_Status, $Mutation_Status, $Validation_Method, $Sequencing_Phase, $Sequence_Source, $Score, $BAM_file, $Sequencer, $chromosome_name, $start, $stop, $reference, $variant, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error) = split(/\t/, $line);

    # Find annotation file headers
    my $annot_line;
    my @annot_columns;
    if ($self->annotation_has_header) {
        $annot_line = <ANNOT>; chomp $annot_line;
        @annot_columns = split(/\t/, $annot_line);
    }
    else {
        @annot_columns = qw(chromosome_name start stop reference variant type gene_name transcript_name transcript_species transcript_source transcript_version strand transcript_status trv_type c_position amino_acid_change ucsc_cons domain all_domains deletion_substructures transcript_error);
    }

    my %annotation_hash;
    while (my $line = <ANNOT>) {
    	chomp($line);
        my ($chromosome_name, $start, $stop, $reference, $variant, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error) = split(/\t/, $line);
	    #my $variant_name = $gene_name."_".$chromosome_name."_".$start."_".$stop."_".$reference."_".$variant;
	    my $variant_name = $chromosome_name."_".$start."_".$stop."_".$reference."_".$variant; #this only matched vcf format for SNVs
    	$annotation_hash{$variant_name} = "$line";
    }

    my @sample_names;
    my $ref_build;
    while (my $line = <VCF>) {
        chomp($line);
        if ($line =~ /^\#\#/) {
            if ($line =~ /^\#\#reference=/) {
                ($ref_build) = $line =~ /^\#\#reference=(.+$)/;
            }
            next;
        }
        elsif ($line =~ /^\#CHROM/) { #grab sample names off of the header line
            my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);
            @sample_names = @samples;
            next;
        }

        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

        #parse format line to find out where our pattern of interest is in the ::::: system
        my (@format_fields) = split(/:/, $format);
        my $gt_location; #genotype
        my $dp_location; #depth - this is not filled for vasily ucla files - this is currently unused down below
        my $gq_location; #genotype quality - this is only filled in washu vcf for samtools variants, not varscan and not homo ref - this is currrently unused down below
        my $count = 0;
        foreach my $format_info (@format_fields) {
            if ($format_info eq 'GT') {
                $gt_location = $count;
            }
            elsif ($format_info eq 'DP') {
                $dp_location = $count;
            }
            elsif ($format_info eq 'GQ') {
                $gq_location = $count;
            }
            $count++;
        }

        #this file doesn't work if there are unknown genotype locations
        unless ($gt_location || $gt_location == 0) {
            die "Format field doesn't have a GT entry, failed to get genotype for $line\n";
        }

        #check to see if line has 0,1,2,etc as genotype numbering, store those in a hash for future reference

        my %alleles_hash;
        foreach my $sample_info (@samples) {                    
            my (@sample_fields) = split(/:/, $sample_info);
            my $genotype = $sample_fields[$gt_location];
            my $allele1 = my $allele2 = ".";
            ($allele1, $allele2) = split(/\//, $genotype);
            if ($allele1 =~ m/\d+/) {
                $alleles_hash{$allele1}++;
                if ($allele2 =~ m/\d+/) {
                    $alleles_hash{$allele2}++;
                }
            }
        }

        my @allele_options = (sort { $a <=> $b } keys %alleles_hash);
        $count = 0;
        foreach my $sample_info (@samples) {
            my $maf = {};
            my (@sample_fields) = split(/:/, $sample_info);
            my $genotype = $sample_fields[$gt_location];
            if ($genotype eq '0/0' || $genotype eq '.') { #MAF DOESN'T HAVE SPOTS FOR REFERENCE GENOTYPES OR MISSING DATA
                next;
            }
            my $allele1 = my $allele2 = ".";
            my $allele_count;
            ($allele1, $allele2) = split(/\//, $genotype);

            my $variant_name;
            my $pos_stop;
            my $variant_type;
            if ($alt =~ /,/) {
                #multi-allelic variant case
                $variant_type = "Multi-alleleic";
                $pos_stop = $pos;
                $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
            } elsif(length($ref) == 1 and length($alt) == 1) {
                #SNV case
                $variant_type = "SNP";
                $pos_stop = $pos;
                if ($allele1 =~ m/\D+/) {
                    $allele_count = '.';
                    $variant_name = 'missing';
                }
                elsif ($allele1 == $allele2) { #homo
                    if ($allele1 == $allele_options[0]) { #homo first variant
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $ref;
                        $allele_count = 0;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    elsif ($allele1 == $allele_options[1]) { #homo second variant
                    	$maf->{Match_Norm_Seq_Allele1} = $alt;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $allele_count = 2;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                }
                else { #heterozygous
                    if ($alt =~ /,/) { #THIS WON'T WORK
                        my ($alt_ref, $alt_alt) = split(/,/, $alt);
                    	$maf->{Match_Norm_Seq_Allele1} = $alt_ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt_alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_ref";
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_alt";
                    }
                    else {
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    $allele_count = 1;
                }
            } elsif (length($ref) == 1 and length($alt) > 1) {
                #insertion case
                $variant_type = "INS";
                $pos_stop = $pos + 1; #VCF uses 1-based position of base before the insertion (which is the same as 0-based position of first inserted base), insertions have no length -- +1 for annotation format
                $ref = '-';
                $alt = substr($alt, 1);
                if ($allele1 =~ m/\D+/) {
                    $allele_count = '.';
                    $variant_name = 'missing';
                }
                elsif ($allele1 == $allele2) { #homo
                    if ($allele1 == $allele_options[0]) { #homo first variant
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $ref;
                        $allele_count = 0;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    elsif ($allele1 == $allele_options[1]) { #homo second variant
                    	$maf->{Match_Norm_Seq_Allele1} = $alt;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $allele_count = 2;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                }
                else { #heterozygous
                    if ($alt =~ /,/) { #THIS WON'T WORK
                        my ($alt_ref, $alt_alt) = split(/,/, $alt);
                    	$maf->{Match_Norm_Seq_Allele1} = $alt_ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt_alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_ref";
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_alt";
                    }
                    else {
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    $allele_count = 1;
                }
            } elsif (length($ref) > 1 and length($alt) == 1) {
                #deletion case
                $variant_type = "DEL";
                $ref = substr($ref, 1);
                $pos++;
                $pos_stop = $pos + length($ref);
                $alt = '-';
                if ($allele1 =~ m/\D+/) {
                    $allele_count = '.';
                    $variant_name = 'missing';
                }
                elsif ($allele1 == $allele2) { #homo
                    if ($allele1 == $allele_options[0]) { #homo first variant
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $ref;
                        $allele_count = 0;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    elsif ($allele1 == $allele_options[1]) { #homo second variant
                    	$maf->{Match_Norm_Seq_Allele1} = $alt;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $allele_count = 2;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                }
                else { #heterozygous
                    if ($alt =~ /,/) { #THIS WON'T WORK
                        my ($alt_ref, $alt_alt) = split(/,/, $alt);
                    	$maf->{Match_Norm_Seq_Allele1} = $alt_ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt_alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_ref";
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt_alt";
                    }
                    else {
                    	$maf->{Match_Norm_Seq_Allele1} = $ref;
                    	$maf->{Match_Norm_Seq_Allele2} = $alt;
                        $variant_name = "$chr"."_"."$pos"."_"."$pos_stop"."_"."$ref"."_"."$alt";
                    }
                    $allele_count = 1;
                }
            } else {
                die $self->error_message('Unhandled variant type encountered');
            }

            my @alts;
            if ($alt =~ /,/) {
                @alts = split(/,/, $alt);
            }
            else {
                @alts = ($alt);
            }

            foreach my $alt_allele (@alts) { #THIS WONT WORK FOR SPLIT ALTS WHERE THE REF DOESN'T EXIST
            	$maf->{Match_Norm_Seq_Allele2} = $alt_allele;
                # Temporary additions to make the SMG test work for MRSA data (HG Data doesn't have a tumor allele)
                $maf->{Tumor_Seq_Allele1} = $maf->{Match_Norm_Seq_Allele1};
                $maf->{Tumor_Seq_Allele2} = $maf->{Match_Norm_Seq_Allele2};

                my $sample_name = $sample_names[$count];
                my ($chromosome_name, $start, $stop, $reference, $variant, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change,  $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error);
                if ($variant_name eq 'missing') {
                    warn "Skipping missing value for variant, does not have annotation\n";
                    next; #Skip variants with no annotation -- MAYBE RUN ANNOTATION AGAIN???
                }
                elsif (defined $variant_name && $annotation_hash{$variant_name}) {
                    ($chromosome_name, $start, $stop, $reference, $variant, $type, $gene_name, $transcript_name, $transcript_species, $transcript_source, $transcript_version, $strand, $transcript_status, $trv_type, $c_position, $amino_acid_change, $ucsc_cons, $domain, $all_domains, $deletion_substructures, $transcript_error) = split(/\t/,$annotation_hash{$variant_name});
                }
                else { 
                    warn "Variant $variant_name in $annotation_file does not have annotation, skipping this variant\n";
                    next; #Skip variants with no annotation -- MAYBE RUN ANNOTATION AGAIN???
                    $chromosome_name = $chr;
                    $start = $pos;
                    $stop = $pos_stop;
                    $reference = $ref;
                    $variant = $alt;
                    $type = $variant_type;
                    $gene_name = $transcript_name = $transcript_species = $transcript_source = $transcript_version = $strand = $transcript_status = $trv_type = $c_position = $amino_acid_change = $ucsc_cons = $domain = $all_domains = $deletion_substructures = $transcript_error = "-";
                }

            	$maf->{Chromosome} = $chromosome_name;
                $maf->{Start_Position} = $start;
                $maf->{End_Position} = $stop;

                $maf->{Reference_Allele} = $reference;

                $maf->{Variant_Type} = $type; # SNP, INS, DEL, etc.
                $maf->{Hugo_Symbol} = $gene_name;
    #taken from gmt annotate revise-maf
                my $entrez_gene_ids;
                my $Hugo_Symbol = $maf->{Hugo_Symbol};
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
                    $maf->{Entrez_Gene_Id} = $entrez_gene_ids->{$Hugo_Symbol};
                }
                else {
                    $maf->{Entrez_Gene_Id} = "-";
                }

                $maf->{Strand} = $strand;
                $maf->{Variant_Classification} = &trv_to_mutation_type($trv_type);
                warn "Unrecognized trv_type \"$trv_type\" in $sample_name data $sample_info annotation file $annotation_hash{$variant_name}: $maf->{Hugo_Symbol}, chr$maf->{Chromosome}:$maf->{Start_Position}-$maf->{End_Position}\n" if ! $maf->{Variant_Classification};

                $maf->{NCBI_Build} = $ref_build;

                if (defined ($id)) {
                    $maf->{dbSNP_RS} = $id;
                    $maf->{dbSNP_Val_Status} = $id;
                }
                else {
                    $maf->{dbSNP_RS} = "-";
                    $maf->{dbSNP_Val_Status} = "-";
                }

	            # required to match correctly
                $maf->{Tumor_Sample_Barcode} = $sample_name;
                $maf->{Matched_Norm_Sample_Barcode} = $sample_name; # Ex. H_MRS-6201-1025127

                #did we pass filters?
                if ($filter eq "PASS" || $filter eq "-" || !defined($filter)) {
                    $maf->{Verification_Status} = "Strandfilter_Passed";
                }
                else {
                    $maf->{Verification_Status} = $filter;
                }

                #load maf with default names
                $maf->{Center} = 'genome.wustl.edu';
                $maf->{Tumor_Validation_Allele1} = "-";
                $maf->{Tumor_Validation_Allele2} = "-";
                $maf->{Match_Norm_Validation_Allele1} = "-";
                $maf->{Match_Norm_Validation_Allele2} = "-";
                $maf->{Validation_Status} = "Unknown";
                $maf->{Mutation_Status} = "Germline";
                $maf->{Validation_Method} = "4";
                $maf->{Sequencing_Phase} = "Capture";
                $maf->{Sequence_Source} = "-";
                $maf->{Score} = "-";
                $maf->{BAM_File} = "-";
                $maf->{Sequencer} = "GaIIx or HiSeq";

                #fill out non-standard maf fields with our annotation
                $maf->{chromosome_name_WU} = $maf->{Chromosome};
                $maf->{start_WU} = $maf->{Start_Position};
                $maf->{stop_WU} = $maf->{End_Position};
                $maf->{reference_WU} = $maf->{Reference_Allele};
                $maf->{variant_WU} = $alt;
                $maf->{type_WU} = $maf->{Variant_Type};
                $maf->{gene_name_WU} = $maf->{Hugo_Symbol};
                $maf->{transcript_name_WU} = $transcript_name;
                $maf->{transcript_species_WU} = $transcript_species;
                $maf->{transcript_source_WU} = $transcript_source;
                $maf->{transcript_version_WU} = $transcript_version;
                $maf->{strand_WU} = $strand;
                $maf->{transcript_status_WU} = $transcript_status;
                $maf->{trv_type_WU} = $maf->{Variant_Classification};
                $maf->{c_position_WU} = $c_position;
                $maf->{amino_acid_change_WU} = $amino_acid_change;
                $maf->{ucsc_cons_WU} = $ucsc_cons;
                $maf->{domain_WU} = $domain;
                $maf->{all_domains_WU} = $all_domains;
                $maf->{deletion_substructures_WU} = $deletion_substructures;
                $maf->{transcript_error_WU} = $transcript_error;

                # print it out to the file
                if ($self->remove_silent && $maf->{Variant_Classification} eq 'Silent') {
                    next;
                }
                # First one without tab
                print MAF (defined $maf->{$maf_columns[0]} ? $maf->{$maf_columns[0]} : "--");
                foreach my $column (@maf_columns[1 .. $#maf_columns]) {
                	print MAF "\t", (defined $maf->{$column} ? $maf->{$column} : "--");
                }
                print MAF "\n";
            }
        }
    }
    return 1; # No error
}


################################################################################
# subroutines                                                                  #
################################################################################

#############################################################
# trv_to_mutation_type - Converts WU var types to MAF variant classifications
#
#############################################################
sub trv_to_mutation_type {
    my $trv_type = shift;
  
    return( "Missense_Mutation" ) if( $trv_type eq "missense" );
    return( "Nonsense_Mutation" ) if( $trv_type eq "nonsense" || $trv_type eq "nonstop" );
    return( "Silent" ) if( $trv_type eq "silent" );
    return( "Splice_Site" ) if( $trv_type eq "splice_site" || $trv_type eq "splice_site_del" || $trv_type eq "splice_site_ins" );
    return( "Frame_Shift_Del" ) if( $trv_type eq "frame_shift_del" );
    return( "Frame_Shift_Ins" ) if( $trv_type eq "frame_shift_ins" );
    return( "In_Frame_Del" ) if( $trv_type eq "in_frame_del" );
    return( "In_Frame_Ins" ) if( $trv_type eq "in_frame_ins" );
    return( "RNA" ) if( $trv_type eq "rna" );
    return( "3'UTR" ) if( $trv_type eq "3_prime_untranslated_region" );
    return( "5'UTR" ) if( $trv_type eq "5_prime_untranslated_region" );
    return( "3'Flank" ) if( $trv_type eq "3_prime_flanking_region" );
    return( "5'Flank" ) if( $trv_type eq "5_prime_flanking_region" );
  
    return( "Intron" ) if( $trv_type eq "intronic" || $trv_type eq "splice_region" || $trv_type eq "splice_region_ins" || $trv_type eq "splice_region_del");
    return( "Targeted_Region" ) if( $trv_type eq "-" );
  
    return( "" );
}

1;
