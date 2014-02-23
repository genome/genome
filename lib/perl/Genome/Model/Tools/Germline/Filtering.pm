package Genome::Model::Tools::Germline::Filtering;

use warnings;
use strict;
use IO::File;
use Genome;

class Genome::Model::Tools::Germline::Filtering {
    is => 'Command',
    has_input => [
    variant_file => {
        is => 'String',
        doc => "MAF, VCF, or WU file of germline variants that require filtering (Valid file extensions: maf, vcf, var). WU variant files must use the standard 5-column format or 21-column annotation format, with an extra column specifying the sample containing the variant (or a comma-delimited list of samples)",
    },
    output_file => {
        is => 'String',
        doc => "WU annotation file of germline variants that passed all filtering steps. Contains extra column with sample names",
    },
    ],
    has_optional_input => [
    reference_transcripts => {
        is => 'String', example_values => ["NCBI-human.ensembl/67_37l_v2"],
        doc => "The annotation build to use with the WU annotator, and to find transcript lengths",
    },
    problem_genes => {
        is => 'String', default => "PDE4DIP,CDC27,MUC4,DUX4",
        doc => "Comma-delimited list of genes whose variants should be filtered out (Olfactory receptor genes are always filtered out)",
    },
    liftover_hg18_to_hg19 => {
        is => 'Boolean', default => 0,
        doc => "If the input variants are on Build36, liftOver and re-annotate to Build37 before applying filters",
    },
    exclude_3prime_vars => {
        is => 'Boolean', default => 1,
        doc => "Sites near the 3' end of a transcript (within 5% of transcript length) are filtered out (unless in a functional domain or splice site)",
    },
    max_mut_freq => {
        is => 'Number', default => 0.02,
        doc => "Variants reported in more than this fraction of cases in the cohort, are removed",
    },
    max_allele_freq => {
        is => 'Number', default => 0.01,
        doc => "Variants reported in more than this fraction of alleles in NHLBI or 1000G, are removed",
    },
    num_cases_in_cohort => {
        is => 'Number',
        doc => "Number of cases in the cohort. This is a required input if --variant-file is a MAF or WU variant file",
    },
    ],
    doc => "A germline variant filtering tool",
};

sub help_synopsis {
    return <<EOS
  gmt germline filtering --variant-file 248_ucec_cases_germline.maf --output-file 248_ucec_cases_germline.anno
EOS
}

sub help_detail {
    return <<EOS
  This tool takes a germline variant list and filters out variants that are likely non-functional.
  Since filters are designed for Build37 loci only, use --liftover-hg18-to-hg19 if necessary.

  The recommended workflow prior to using this tool:
  1) Create a model-group of succeeded somatic-variation models using the TN-pairs under study
  2) Run 'gmt capture somatic-variation-group --output-germline-calls' on the model-group (this
     pulls germline variants from the models, runs standard false-positive-filters, and generates
     WU annotation files per sample)
  3) Concatenate all the resulting WU annotation files into a single file, retaining sample names
     in a 22nd column, and use it as the input to this tool
EOS
}

sub _doc_authors {
    return <<EOS
  Cyriac Kandoth, Ph.D.
  Nathan D. Dees, Ph.D.
EOS
}

sub execute {

    my $self = shift;
    $DB::single = 1;

    my ( $variant_file, $output_file ) = ( $self->variant_file, $self->output_file );
    my %problem_genes = map{($_,1)} split( /,/, $self->problem_genes );
    my ( $max_mut_freq, $max_allele_freq ) = ( $self->max_mut_freq, $self->max_allele_freq );
    my $exclude_3prime_vars = $self->exclude_3prime_vars;

    # Parse out the variants into a hash with WU annotations
    $self->debug_message( "\nParsing the variant file...\n" );
    my %vars = my %all_samples;
    ( -s $variant_file ) or die "Input file does not exist or has zero size\n";
    if( $variant_file =~ m/\.maf$/ ) {
        ( defined $self->num_cases_in_cohort ) or die "Please specify '--num-cases-in-cohort' when input is a MAF\n";
        %vars = $self->parse_maf( $variant_file, \%all_samples );
    }
    elsif( $variant_file =~ m/\.vcf$/ ) {
        %vars = $self->parse_vcf( $variant_file, \%all_samples );
    }
    elsif( $variant_file =~ m/\.var$/ ) {
        ( defined $self->num_cases_in_cohort ) or die "Please specify '--num-cases-in-cohort' when input is a WU variant file\n";
        %vars = $self->parse_var( $variant_file, \%all_samples );
    }
    else {
        die "The input file must use an extension: .maf, .vcf, or .var\n";
    }

    # We need to count the total number of samples in the variant list, unless the user provided it
    my $total_num_samples = scalar( keys %all_samples );
    $total_num_samples = $self->num_cases_in_cohort if( defined $self->num_cases_in_cohort );

    # Check if there are enough samples to apply the max_mut_freq threshold
    if( $max_mut_freq * $total_num_samples < 1 ) {
        die "With max-mut-freq = $max_mut_freq, and only $total_num_samples samples, zero variants will pass the filters\n";
    }
    elsif( $max_mut_freq * $total_num_samples < 2 ) {
        $self->status_message( "WARNING: With max-mut-freq = $max_mut_freq, and only $total_num_samples samples, no recurrent variants will pass the filters\n" );
    }

    $self->debug_message( "\nLoading transcript lengths based on " . $self->reference_transcripts . "...\n" );
    my %tr_95pc_of_length;
    my $anno_dir = Genome::Model::Build::ImportedAnnotation->get( name=>$self->reference_transcripts )->data_directory;
    map{my @c = split(/\t|,/); if($c[8] eq "utr_exon" and $c[17] ne "0"){$tr_95pc_of_length{$c[30]} = ($c[17]*0.95);}}`cat $anno_dir/annotation_data/substructures/*.csv`;

    $self->debug_message( "\nLoading 1000G and NHLBI data...\n" );
    # Locate the tabix-indexed file containing variant allele frequencies and dbSNP IDs from 1000G
    my $onekg_file = '/gscmnt/gc6132/info/medseq/1000_genomes/downloads/2012-03-27/ALL.wgs.phase1_release_v3.20101123.snps_indels.sites.vars.gz';
    # Load the NHLBI-exome variant allele frequencies and dbSNP IDs into a hash
    my $nhlbi_file = '/gscmnt/sata170/info/medseq/NHLBI/ESP5400/NHLBI_ESP5400.snv.vars';
    my %nhlbi_af;
    my $nhlbi_fh = IO::File->new( $nhlbi_file ) or die "Couldn't open $nhlbi_file.";
    $nhlbi_fh->getline; # Skip header line
    while( my $line = $nhlbi_fh->getline ) {
        my @col = split( /\t/, $line );
        my $key = join( "\t", @col[0..4] );
        $nhlbi_af{$key} = $col[7];
    }
    $nhlbi_fh->close;

    # For each variant, run the filters and store pass/fail status
    $self->debug_message( "\nRunning filters on each variant...\n" );
    foreach my $key ( keys %vars ) {
        my @wu_anno = @{$vars{$key}{wu_anno}};
        my ( $var_type, $gene, $tr_name, $tr_source, $var_class, $c_pos, $domain, $tr_errors ) = ( @wu_anno[5..7,9,13,14,17,20] );

        $vars{$key}{filter} = "";

        # Tag variants annotated to transcripts that are not entirely error-free
        $vars{$key}{filter} .= "transcript_error;" unless( $tr_errors eq 'no_errors' );

        # Tag all LOC/ENSG/Corf genes and XM/XR transcripts
        $vars{$key}{filter} .= "LOC_ENSG_Corf_gene;" if( $gene =~ m/^(LOC\d+|ENSG\d+|C\w+orf\d+)$/ );
        $vars{$key}{filter} .= "XM_XR_transcript;" if( $tr_name =~ m/^(XM|XR)\_\d+/ );

        # Tag variants on user-defined problematic genes and Olfactory receptor (OR*) genes
        $vars{$key}{filter} .= "problem_gene;" if( defined $problem_genes{$gene} or $gene =~ m/^OR\d+\w+\d+$/ );

        # Tag variants in non-coding regions
        $vars{$key}{filter} .= "non_coding;" unless( $var_class =~ m/^(frame_shift|in_frame|splice_site|missense$|nonsense$|nonstop$|silent$)/ );

        # Tag variants near the 3' end of the transcript, as long as it's not in a known functional domain or splice site
        if( $exclude_3prime_vars ) {
            if( $domain eq 'NULL' and $var_class =~ m/^(frame_shift|in_frame|missense$|nonsense$|nonstop$|silent)/  ) {
                ( my $c_pos_value ) = $c_pos =~ /[._](\d+)$/;
                if( defined $tr_95pc_of_length{$tr_name} ) {
                    ( $c_pos_value =~ m/^\d+$/ ) or die "Cannot parse out nucleotide position from $c_pos";
                    $vars{$key}{filter} .= "near_3prime_end;" if( $c_pos_value > $tr_95pc_of_length{$tr_name} );
                }
            }
        }

        # Tag common variants - that are seen in a sufficiently large subset of cases in the cohort
        my $num_cases = scalar( @{$vars{$key}{samples}} );
        $vars{$key}{filter} .= "common_in_cohort;" if( $num_cases / $total_num_samples > $max_mut_freq );

        # Tag variants sufficiently common in the 1000G dataset
        my ( $chr, $start, $stop, $ref, $var ) = split( /\t/, $key );
        my $var_freq_1000g;
        open( TABIX_PIPE, "tabix $onekg_file $chr:$start-$stop |" );
        while( my $line = <TABIX_PIPE> ) {
            my @col = split( /\t/, $line );
            $var_freq_1000g = $col[7] if( $key eq join( "\t", @col[0..4] ));
        }
        close( TABIX_PIPE );
        $vars{$key}{filter} .= "common_in_1000G;" if( defined $var_freq_1000g && $var_freq_1000g > $max_allele_freq );

        # Tag variants sufficiently common in the NHLBI dataset
        $vars{$key}{filter} .= "common_in_NHLBI;" if( defined $nhlbi_af{$key} && $nhlbi_af{$key} > $max_allele_freq );

        # ::TODO:: Run VEP annotation, and apply additional filters based on that

        $vars{$key}{filter} = 'PASS' if( $vars{$key}{filter} eq "" );
    }

    # Summarize the most common reasons why variants were filtered out
    $self->status_message( scalar( keys %vars ) . " unique variants across " . $total_num_samples . " cases ran through filters as follows:\n" );
    my %counts;
    ++$counts{$vars{$_}{filter}} foreach( keys %vars );
    foreach my $filter ( sort {$counts{$b} <=> $counts{$a}} keys %counts ) {
        $self->status_message( "$counts{$filter}" . "\t$filter\n" );
    }

    # Create a WU anno file of variants that passed all filters, and the samples they were seen in
    my $out_fh = IO::File->new( $output_file, ">" ) or die "Couldn't open $output_file.";
    $out_fh->print( "chromosome_name\tstart\tstop\treference\tvariant\ttype\tgene_name\t" .
    "transcript_name\ttranscript_species\ttranscript_source\ttranscript_version\tstrand\t" .
    "transcript_status\ttrv_type\tc_position\tamino_acid_change\tucsc_cons\tdomain\t" .
    "all_domains\tdeletion_substructures\ttranscript_error\tsample\n" );
    foreach my $key ( sort keys %vars ) {
        my @wu_anno = @{$vars{$key}{wu_anno}};
        if( $vars{$key}{filter} eq 'PASS' ) {
            foreach my $sample ( @{$vars{$key}{samples}} ) {
                $out_fh->print( join( "\t", @wu_anno ) . "\t$sample\n" );
            }
        }
    }
    $out_fh->close;

    return 1;
}

# Parse a MAF and return a hashful of WU annotations per variant
sub parse_maf {
    my ( $self, $maf_file, $all_samples_ref ) = @_;
    my ( %vars, %dedup_samples );
    my $line = "";
    my $need_wu_anno = 0;

    # Open the MAF and handle the headers before looping through each variant
    my $maf_fh = IO::File->new( $maf_file ) or die "Couldn't open $maf_file.";
    while( $line = $maf_fh->getline ) { # Skip all comment lines
        last if( $line !~ m/^(#|$)/ );
    }
    chomp( $line );
    my @col = split( /\t/, $line );
    if( $col[0] eq 'Hugo_Symbol' ) {
        $need_wu_anno = 1 unless( defined $col[52] and $col[52] eq 'transcript_error' );
    }
    else {
        die "Provided MAF file does not appear to have a header line!\n";
    }

    # Now we can assume the rest of the file is just variants, and parse it quicker
    while( $line = $maf_fh->getline ) {
        chomp( $line );
        @col = split( /\t/, $line );

        # Parse the line, and fill the hash indexed by locus & ref/var
        my $sample = $col[15];
        my $key = join( "\t", @col[4..6,10] ) . "\t" . ( $col[10] eq $col[11] ? $col[12] : $col[11] );
        @{$vars{$key}{wu_anno}} = @col[32..52] unless( defined $vars{$key} or $need_wu_anno );

        # Keep track of all the samples each variant is reported in
        unless( defined $dedup_samples{$key}{$sample} ) {
            push( @{$vars{$key}{samples}}, $sample );
            $dedup_samples{$key}{$sample} = 1;
            $all_samples_ref->{$sample} = 1;
        }
    }
    $maf_fh->close;

    # Before we proceed, make sure the user isn't pulling our leg
    my $sample_count = scalar( keys %{$all_samples_ref} );
    if( defined $self->num_cases_in_cohort && $sample_count > $self->num_cases_in_cohort ) {
        die "$sample_count samples found in MAF, more that what's specified in '--num-cases-in-cohort'!\n";
    }

    # Write the variants in a 5-col format, in case we need to liftOver and/or annotate
    my $var_5col_list = Genome::Sys->create_temp_file_path;
    my $var_5col_list_sorted = Genome::Sys->create_temp_file_path;
    ( $var_5col_list and $var_5col_list_sorted ) or die "Couldn't create a temp file. $!";
    if( $self->liftover_hg18_to_hg19 or $need_wu_anno ) {
        my $outFh = IO::File->new( $var_5col_list, ">" ) or die "Couldn't open $var_5col_list.";
        foreach my $key ( sort keys %vars ) {
            $outFh->print( "$key\t" . join( ",", @{$vars{$key}{samples}} ) . "\n" );
        }
        $outFh->close;

        # Sort the variant list by genomic loci, so that the WU annotator can run quicker
        Genome::Sys->shellcmd( cmd => "joinx sort --stable --input-file $var_5col_list --output-file $var_5col_list_sorted" );
    }

    if( $self->liftover_hg18_to_hg19 ) {
        %vars = $self->liftover_to_hg19( $var_5col_list_sorted );
    }
    elsif( $need_wu_anno ) {
        %vars = $self->annotate_to_build37( $var_5col_list_sorted );
    }
    return %vars;
}

# Parse a VCF, run WU annotator, and return a hashful of WU annotations per variant
sub parse_vcf {
    my ( $self, $vcf_file, $all_samples_ref ) = @_;
    my ( @sample_names, @vcf_header );

    # Reformat each variant for the WU annotator
    my $var_5col_list = Genome::Sys->create_temp_file_path;
    ( $var_5col_list ) or die "Couldn't create a temp file. $!";
    my $outFh = IO::File->new( $var_5col_list, ">" ) or die "Couldn't open $var_5col_list.";
    my $vcf_fh = IO::File->new( $vcf_file ) or die "Couldn't open $vcf_file.";
    while( my $line = $vcf_fh->getline ) {

        chomp( $line );
        my ( $chr, $pos, $rsid, $ref, $alt_line, $qual, $filter, $info_line, $format_line, @rest ) = split( /\t/, $line );

        # Skip headers, but if this is the line containing sample names, parse those out for later
        next if( $line =~ m/^(##|$)/ );
        if( $line =~ m/^#CHROM/ ) {
            @sample_names = @rest;
            my $sample_count = scalar( @sample_names );
            # Before we proceed, make sure the user isn't pulling our leg
            if( defined $self->num_cases_in_cohort && $sample_count != $self->num_cases_in_cohort ) {
                die "$sample_count samples found in VCF, differs from what's specified in '--num-cases-in-cohort'!\n";
            }
            %{$all_samples_ref} = map{( $_, 1 )} @sample_names;
            next;
        }

        # If there are multiple alt alleles, set them up separately for annotation
        my @alt_alleles = split( /,/, $alt_line );
        foreach my $alt( @alt_alleles ) {

            # Parse out allele depths from the formatted columns, and find samples with this variant allele
            my @format_keys = split( /\:/, $format_line );
            my @var_samples;
            for( my $i = 0; $i < scalar( @rest ); ++$i ) {
                my ( $j, $k ) = ( 0, 0 );
                my %format_info = map {( $format_keys[$j++], $_ )} split( /\:/, $rest[$i] );

                # If allele depth is not provided, we can't do anything. Just skip the variant
                my %allele_depth;
                if( !defined $format_info{AD} or $format_info{AD} eq '.' ) {
                    %allele_depth = map {( $_, 0 )} @alt_alleles;
                }
                else {
                    %allele_depth = map {( $alt_alleles[$k++], $_ )} split( /,/, $format_info{AD} );
                }
                if( $allele_depth{$alt} > 0 ) {
                    push( @var_samples, $sample_names[$i] );
                }
            }

            # Parse the info column and reformat the variant for the WU annotator
            my %info = map {(m/=/ ? (split(/=/)) : ($_,1))} split( /\;/, $info_line );
            my $start = my $stop = my $var = "";
            my ( $ref_length, $alt_length, $indel_size ) = ( length( $ref ), length( $alt ), 0 );

            if( !defined $info{VT} && $ref_length == 1 && $alt_length == 1 ) { # Handle a SNP-only VCF
                ( $start, $stop ) = ( $pos, $pos );
                $var = $alt;
            }
            elsif( $info{VT} eq 'SNP' ) { # Handle SNPs
                ( $start, $stop ) = ( $pos, $pos );
                $var = $alt;
            }
            elsif( $info{VT} eq 'INDEL' && $ref_length < $alt_length ) { # Handle insertions
                $indel_size = $alt_length - $ref_length;
                ( $ref, $var ) = ( "-", substr( $alt, $ref_length, $indel_size ));
                ( $start, $stop ) = ( $pos, $pos + 1 );
            }
            elsif( $info{VT} eq 'INDEL' && $ref_length > $alt_length ) { # Handle deletion
                $indel_size = $ref_length - $alt_length;
                ( $ref, $var ) = ( substr( $ref, $alt_length, $indel_size ), "-" );
                ( $start, $stop ) = ( $pos + 1, $pos + $indel_size );
            }
            elsif( $info{VT} eq 'SV' ) { # Skip SVs
                next;
            }
            else { # Poop-out on unknown variant types
                die "Unhandled variant type in VCF:\n$line\nPlease check parser!\n";
            }

            # Print variants to a temporary file, keeping track of the samples containing them
            $outFh->print( "$chr\t$start\t$stop\t$ref\t$var\t" . join( ",", @var_samples ) . "\n" );
        }
    }
    $vcf_fh->close;
    $outFh->close;

    # Sort the variant list by genomic loci, so that the WU annotator can run quicker
    my $var_5col_list_sorted = Genome::Sys->create_temp_file_path;
    ( $var_5col_list_sorted ) or die "Couldn't create a temp file. $!";
    Genome::Sys->shellcmd( cmd => "joinx sort --stable --input-file $var_5col_list --output-file $var_5col_list_sorted" );

    my %vars;
    if( $self->liftover_hg18_to_hg19 ) {
        %vars = $self->liftover_to_hg19( $var_5col_list_sorted );
    }
    else {
        %vars = $self->annotate_to_build37( $var_5col_list_sorted );
    }
    return %vars;
}

# Parse a WU variant file and return a hashful of WU annotations per variant
sub parse_var {
    my ( $self, $var_file, $all_samples_ref ) = @_;
    my ( %dedup_samples, %vars );
    my ( $need_wu_anno, $sample_name_col_idx ) = ( 0, 21 );

    # Check for problems in the WU variant file, and store it into a hash
    my $var_fh = IO::File->new( $var_file ) or die "Couldn't open $var_file.";
    while( my $line = $var_fh->getline ) {
        chomp( $line );
        next if( $line =~ m/^(#|chromosome_name|$)/ );
        my @col = split( /\t/, $line );

        # Using the first line in the file, find out if it is missing WU annotation columns
        unless( %vars ) {
            if( scalar( @col ) == 6 ) {
                ( $need_wu_anno, $sample_name_col_idx ) = ( 1, 5 );
            }
            elsif( scalar( @col ) != 22 ) {
                die "Unrecognized WU variant file format! Use the standard 5-column format or 21-column annotation format, with an extra column specifying the sample containing that variant (or a comma-delimited list of samples)\n";
            }
        }

        # Fill the hash indexed by locus & ref/var
        my $key = join( "\t", @col[0..4] );
        @{$vars{$key}{wu_anno}} = @col[0..20] unless( defined $vars{$key} or $need_wu_anno );

        # Keep track of all the samples each variant is reported in
        my @samples = split( /,/, $col[$sample_name_col_idx] );
        foreach my $sample ( @samples ) {
            unless( defined $dedup_samples{$key}{$sample} ) {
                push( @{$vars{$key}{samples}}, $sample );
                $dedup_samples{$key}{$sample} = 1;
                $all_samples_ref->{$sample} = 1;
            }
        }
    }
    $var_fh->close;

    # Before we proceed, make sure the user isn't pulling our leg
    my $sample_count = scalar( keys %{$all_samples_ref} );
    if( defined $self->num_cases_in_cohort && $sample_count > $self->num_cases_in_cohort ) {
        die "$sample_count samples found in MAF, more that what's specified in '--num-cases-in-cohort'!\n";
    }

    # Write the variants in a 5-col format, in case we need to liftOver and/or annotate
    my $var_5col_list = Genome::Sys->create_temp_file_path;
    my $var_5col_list_sorted = Genome::Sys->create_temp_file_path;
    ( $var_5col_list and $var_5col_list_sorted ) or die "Couldn't create a temp file. $!";
    if( $self->liftover_hg18_to_hg19 or $need_wu_anno ) {
        my $outFh = IO::File->new( $var_5col_list, ">" ) or die "Couldn't open $var_5col_list.";
        foreach my $key ( sort keys %vars ) {
            $outFh->print( "$key\t" . join( ",", @{$vars{$key}{samples}} ) . "\n" );
        }
        $outFh->close;

        # Sort the variant list by genomic loci, so that the WU annotator can run quicker
        Genome::Sys->shellcmd( cmd => "joinx sort --stable --input-file $var_5col_list --output-file $var_5col_list_sorted" );
    }

    if( $self->liftover_hg18_to_hg19 ) {
        %vars = $self->liftover_to_hg19( $var_5col_list_sorted );
    }
    elsif( $need_wu_anno ) {
        %vars = $self->annotate_to_build37( $var_5col_list_sorted );
    }
    return %vars;
}

# Given a hashful of Build36 variants, return a similar hash remapped and reannotated to Build37
sub liftover_to_hg19 {
    my ( $self, $build36_file ) = @_;

    # Write all unique variant loci to a temporary file
    my $build37_file = Genome::Sys->create_temp_file_path;
    ( $build37_file ) or die "Couldn't create a temp file. $!";

    # Run liftOver using the wrapper that can handle annotation format files
    $self->debug_message( "\nRunning liftOver to remap the provided Build36 variants to Build37...\n" );
    my $lift_cmd = Genome::Model::Tools::LiftOver->create(
        source_file => $build36_file, destination_file => $build37_file,
        input_is_annoformat => 1, lift_direction => "hg18ToHg19",
    );
    ( $lift_cmd->execute ) or die "Failed to run 'gmt lift-over'";
    $lift_cmd->delete;

    # Run the Build37 WU annotator, and return a hashful of annotated Build37 variants
    return $self->annotate_to_build37( $build37_file );
}

# Return a hashful of WU annotations, given a list of Build37 variants (with comma-delimited sample names) like so:
# 19    51507544    51507547    TATT    -    AML1,AML23,TCGA-AB-1543
sub annotate_to_build37 {
    my ( $self, $build37_file ) = @_;

    # Run the Build37 WU annotator
    $self->debug_message( "\nRunning WU Annotator based on build: ". $self->reference_transcripts . "...\n" );
    my $build37_anno = Genome::Sys->create_temp_file_path;
    ( $build37_anno ) or die "Couldn't create a temp file. $!";
    my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        variant_file => $build37_file, output_file => $build37_anno,
        annotation_filter => "top", reference_transcripts => $self->reference_transcripts,
        extra_columns => "samples", use_version => 2,
    );
    ( $anno_cmd->execute ) or die "Failed to run 'gmt annotate transcript-variants'";
    $anno_cmd->delete;

    # Parse the remapped and reannotated output, and create a new hash to return
    my %b37_vars;
    my $inFh = IO::File->new( $build37_anno ) or die "Couldn't open $build37_anno.";
    while( my $line = $inFh->getline ) {
        next if( $line =~ m/^(#|chromosome_name|$)/ );
        chomp( $line );
        my @col = split( /\t/, $line );

        # Parse the line, and fill the hash indexed by locus & ref/var
        my @wu_anno = @col[0..5,7..21];
        my $samples = $col[6];
        my $key = join( "\t", @col[0..4] );
        @{$b37_vars{$key}{wu_anno}} = @wu_anno;
        @{$b37_vars{$key}{samples}} = split( /,/, $samples );
    }
    $inFh->close;

    return %b37_vars;
}

1;
