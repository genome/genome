package Genome::Model::Tools::Music::Bmr::CalcBmr;

use warnings;
use strict;
use Genome;
use IO::File;
use Bit::Vector;
use List::Util qw( min sum );

our $VERSION = $Genome::Model::Tools::Music::VERSION;

# These constants let us use space-efficient arrays instead of hashes, while keeping the code fairly readable
use constant { AT_Transitions => 0, AT_Transversions => 1, CG_Transitions => 2, CG_Transversions => 3,
               CpG_Transitions => 4, CpG_Transversions => 5, Indels => 6, Truncations => 7, Overall => 8,
               covd_bases => 0, mutations => 1, bmr => 2 };

class Genome::Model::Tools::Music::Bmr::CalcBmr {
    is => 'Genome::Model::Tools::Music::Bmr::Base',
    has_input => [
        roi_file => { is => 'Text', doc => "Tab delimited list of ROIs [chr start stop gene_name] (See DESCRIPTION)" },
        reference_sequence => { is => 'Text', doc => "Path to reference sequence in FASTA format" },
        bam_list => { is => 'Text', doc => "Tab delimited list of BAM files [sample_name normal_bam tumor_bam] (See DESCRIPTION)" },
        output_dir => { is => 'Text', doc => "Directory where output files will be written (Use the same one used with calc-covg)" },
        maf_file => { is => 'Text', doc => "List of mutations using TCGA MAF specification v2.3" },
        bmr_groups => { is => 'Integer', doc => "Number of clusters of samples with comparable BMRs (See DESCRIPTION)", is_optional => 1, default => 1 },
        show_skipped => { is => 'Boolean', doc => "Report each skipped mutation, not just how many", is_optional => 1, default => 0 },
        separate_truncations => { is => 'Boolean', doc => "Group truncational mutations as a separate category", is_optional => 1, default => 0 },
        merge_concurrent_muts => { is => 'Boolean', doc => "Multiple mutations of a gene in the same sample are treated as 1", is_optional => 1, default => 0 },
        genes_to_ignore => { is => 'Text', doc => "Comma-delimited list of genes to ignore for background mutation rates", is_optional => 1 },
        skip_non_coding => { is => 'Boolean', doc => "Skip non-coding mutations from the provided MAF file", is_optional => 1, default => 1 },
        skip_silent => { is => 'Boolean', doc => "Skip silent mutations from the provided MAF file", is_optional => 1, default => 1 },
    ],
    has_output => [
        bmr_output => { is => 'Number', is_optional => 1, doc => "TODO" },
        gene_mr_file => { is => 'Text', is_optional => 1, doc => "TODO" },
    ],
    doc => "Calculates mutation rates given per-gene coverage (from \"music bmr calc-covg\"), and a mutation list",
};

sub help_synopsis {
    return <<HELP
 ... music bmr calc-bmr \\
    --bam-list input_dir/bam_list \\
    --maf-file input_dir/myMAF.tsv \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv

 ... music bmr calc-bmr \\
    --bam-list input_dir/bam_list \\
    --maf-file input_dir/myMAF.tsv \\
    --output-dir output_dir/ \\
    --reference-sequence input_dir/all_sequences.fa \\
    --roi-file input_dir/all_coding_exons.tsv \\
    --genes-to-ignore GENE1,GENE2
HELP
}

sub help_detail {
    return <<HELP;
Given a mutation list (MAF), and per-gene coverage data calculated using \"music bmr calc-covg\"),
this script calculates overall Background Mutation Rate (BMR) and BMRs in the categories of
AT/CG/CpG Transitions, AT/CG/CpG Transversions, and Indels. An optional category for truncational
mutations can also be specified. The script generates a file with per-gene mutation rates that can
be used with the tool that tests for significantly mutated genes (music smg).
HELP
}

sub _additional_help_sections {
    return (
    "ARGUMENTS",
<<EOS

=over 4

=item --roi-file

=over 8

=item The regions of interest (ROIs) of each gene are typically regions targeted for sequencing or
  are merged exon loci (from multiple transcripts) of genes with 2-bp flanks (splice junctions).
  ROIs from the same chromosome must be listed adjacent to each other in this file. This allows
  the underlying C-based code to run much more efficiently and avoid re-counting bases seen in
  overlapping ROIs (for overall covered base counts). For per-gene base counts, an overlapping
  base will be counted each time it appears in an ROI of the same gene. To avoid this, be sure to
  merge together overlapping ROIs of the same gene. BEDtools' mergeBed can help if used per gene.

=back

=item --reference-sequence

=over 8

=item The reference sequence in FASTA format. If a reference sequence index is not found next to
  this file (a .fai file), it will be created.

=back

=item --bam-list

=over 8

=item Provide a file containing sample names and normal/tumor BAM locations for each. Use the tab-
  delimited format [sample_name normal_bam tumor_bam] per line. Additional columns like clinical
  data are allowed, but ignored. The sample_name must be the same as the tumor sample names used
  in the MAF file (16th column, with the header Tumor_Sample_Barcode).

=back

=item --bmr-groups

=over 8

=item Ideally, we want to test the mutation rate (MR) of a gene in a sample against the background
  mutation rate (BMR) across that sample. But if the BMRs of some samples are comparable, we can
  instead test the MR of a gene across a group of samples with comparable BMR, against the overall
  BMR of that group. This argument specifies how many such groups you want to cluster samples
  into. By default, it is assumed that all samples have comparable BMRs (bmr-groups = 1).

=back

=item --output-dir

=over 8

=item This should be the same output directory used when running "music bmr calc-covg". The
  following outputs of this script will also be created/written:
  overall_bmrs: File containing categorized overall background mutation rates.
  gene_mrs: File containing categorized per-gene mutation rates.

=back

=item --genes-to-ignore

=over 8

=item A comma-delimited list of genes to ignore for overall BMR calculations. List genes that are
  known factors in this disease and whose mutations should not be classified as background.

=back

=back

EOS
  );
}

sub _doc_authors {
    return " Cyriac Kandoth, Ph.D.";
}

sub _doc_see_also {
    return <<EOS
B<genome-music-bmr>(1),
B<genome-music>(1),
B<genome>(1)
EOS
}

sub execute {
    my $self = shift;
    my $roi_file = $self->roi_file;
    my $ref_seq = $self->reference_sequence;
    my $bam_list = $self->bam_list;
    my $output_dir = $self->output_dir;
    my $maf_file = $self->maf_file;
    my $show_skipped = $self->show_skipped;
    my $bmr_groups = $self->bmr_groups;
    my $separate_truncations = $self->separate_truncations;
    my $merge_concurrent_muts = $self->merge_concurrent_muts;
    my $genes_to_ignore = $self->genes_to_ignore;
    my $skip_non_coding = $self->skip_non_coding;
    my $skip_silent = $self->skip_silent;

    # Check on all the input data before starting work
    print STDERR "ROI file not found or is empty: $roi_file\n" unless( -s $roi_file );
    print STDERR "Reference sequence file not found: $ref_seq\n" unless( -e $ref_seq );
    print STDERR "List of BAMs not found or is empty: $bam_list\n" unless( -s $bam_list );
    print STDERR "Output directory not found: $output_dir\n" unless( -e $output_dir );
    print STDERR "MAF file not found or is empty: $maf_file\n" unless( -s $maf_file );
    return undef unless( -s $roi_file && -e $ref_seq && -s $bam_list && -e $output_dir && -s $maf_file );

    # Check on the files we expect to find within the provided output directory
    $output_dir =~ s/(\/)+$//; # Remove trailing forward slashes if any
    my $gene_covg_dir = "$output_dir/gene_covgs"; # Should contain per-gene coverage files per sample
    my $total_covgs_file = "$output_dir/total_covgs"; # Should contain overall coverages per sample
    print STDERR "Directory with per-gene coverages not found: $gene_covg_dir\n" unless( -e $gene_covg_dir );
    print STDERR "Total coverages file not found or is empty: $total_covgs_file\n" unless( -s $total_covgs_file );
    return undef unless( -e $gene_covg_dir && -s $total_covgs_file );

    # Outputs of this script will be written to these locations in the output directory
    my $overall_bmr_file = "$output_dir/overall_bmrs";
    my $gene_mr_file = "$output_dir/gene_mrs";
    $self->gene_mr_file( $gene_mr_file );

    # Build a hash to quickly lookup the genes to be ignored for overall BMRs
    my %ignored_genes = ();
    if( defined $genes_to_ignore ) {
        %ignored_genes = map { $_ => 1 } split( /,/, $genes_to_ignore );
    }

    # Parse out the names of the samples which should match the names of the coverage files needed
    my ( @all_sample_names, %sample_idx );
    my $idx = 0;
    my $sampleFh = IO::File->new( $bam_list ) or die "Couldn't open $bam_list. $!";
    while( my $line = $sampleFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample ) = split( /\t/, $line );
        push( @all_sample_names, $sample );
        $sample_idx{$sample} = $idx++;
    }
    $sampleFh->close;

    # If the reference sequence FASTA file hasn't been indexed, do it
    my $ref_seq_idx = "$ref_seq.fai";
    Genome::Sys->shellcmd( cmd => "samtools faidx $ref_seq" ) unless( -e $ref_seq_idx );

    # Parse gene names and ROIs. Mutations outside these ROIs will be skipped
    my ( @all_gene_names, %gene_idx );
    $idx = 0;
    my $roi_bitmask = $self->create_empty_genome_bitmask( $ref_seq_idx );
    my $roiFh = IO::File->new( $roi_file ) or die "Couldn't open $roi_file. $!";
    while( my $line = $roiFh->getline ) {
        next if( $line =~ m/^#/ );
        chomp $line;
        my ( $chr, $start, $stop, $gene ) = split( /\t/, $line );
        if( !$roi_bitmask->{$chr} or $start > $roi_bitmask->{$chr}->Size ) {
            print STDERR "Skipping invalid ROI bitmask $chr:$start-$stop\n";
            next;
        }
        $roi_bitmask->{$chr}->Interval_Fill( $start, $stop );
        unless( defined $gene_idx{$gene} ) {
            push( @all_gene_names, $gene );
            $gene_idx{$gene} = $idx++;
        }
    }
    $roiFh->close;

    # These are the various categories that each mutation will be classified into
    my @mut_classes = ( AT_Transitions, AT_Transversions, CG_Transitions, CG_Transversions, CpG_Transitions, CpG_Transversions, Indels );
    push( @mut_classes, Truncations ) if( $separate_truncations );
    # Save the actual class names for reporting purposes, because the elements above are really just numerical constants
    my @mut_class_names = qw( AT_Transitions AT_Transversions CG_Transitions CG_Transversions CpG_Transitions CpG_Transversions Indels );
    push( @mut_class_names, 'Truncations' ) if( $separate_truncations );

    my @sample_mr; # Stores per sample covg and mutation information
    foreach my $sample ( @all_sample_names ) {
        $sample_mr[$sample_idx{$sample}][$_][mutations] = 0 foreach( @mut_classes );
        $sample_mr[$sample_idx{$sample}][$_][covd_bases] = 0 foreach( @mut_classes );
    }

    # Load the covered base-counts per sample from the output of "music bmr calc-covg"
    print STDERR "Loading per-sample coverages stored in $total_covgs_file\n";
    my $sample_cnt_in_total_covgs_file = 0;
    my $totCovgFh = IO::File->new( $total_covgs_file ) or die "Couldn't open $total_covgs_file. $!";
    while( my $line = $totCovgFh->getline ) {
        next unless( $line =~ m/^\S+\t\d+\t\d+\t\d+\t\d+$/ and $line !~ m/^#/ );
        chomp( $line );
        ++$sample_cnt_in_total_covgs_file;
        my ( $sample, $covd_bases, $covd_at_bases, $covd_cg_bases, $covd_cpg_bases ) = split( /\t/, $line );
        $sample_mr[$sample_idx{$sample}][AT_Transitions][covd_bases] = $covd_at_bases;
        $sample_mr[$sample_idx{$sample}][AT_Transversions][covd_bases] = $covd_at_bases;
        $sample_mr[$sample_idx{$sample}][CG_Transitions][covd_bases] = $covd_cg_bases;
        $sample_mr[$sample_idx{$sample}][CG_Transversions][covd_bases] = $covd_cg_bases;
        $sample_mr[$sample_idx{$sample}][CpG_Transitions][covd_bases] = $covd_cpg_bases;
        $sample_mr[$sample_idx{$sample}][CpG_Transversions][covd_bases] = $covd_cpg_bases;
        $sample_mr[$sample_idx{$sample}][Indels][covd_bases] = $covd_bases;
        $sample_mr[$sample_idx{$sample}][Truncations][covd_bases] = $covd_bases if( $separate_truncations );
    }
    $totCovgFh->close;

    unless( $sample_cnt_in_total_covgs_file == scalar( @all_sample_names )) {
        print STDERR "Mismatching number of samples in $total_covgs_file and $bam_list\n";
        return undef;
    }

    my @gene_mr; # Stores per gene covg and mutation information
    foreach my $gene ( @all_gene_names ) {
        foreach my $sample ( @all_sample_names ) {
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$_][mutations] = 0 foreach( @mut_classes );
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$_][covd_bases] = 0 foreach( @mut_classes );
        }
    }

    # Sum up the per-gene covered base-counts across samples from the output of "music bmr calc-covg"
    print STDERR "Loading per-gene coverage files stored under $gene_covg_dir/\n";
    foreach my $sample ( @all_sample_names ) {
        my $sample_covg_file = "$gene_covg_dir/$sample.covg";
        my $sampleCovgFh = IO::File->new( $sample_covg_file ) or die "Couldn't open $sample_covg_file. $!";
        while( my $line = $sampleCovgFh->getline ) {
            next unless( $line =~ m/^\S+\t\d+\t\d+\t\d+\t\d+\t\d+$/ and $line !~ m/^#/ );
            chomp( $line );
            my ( $gene, undef, $covd_bases, $covd_at_bases, $covd_cg_bases, $covd_cpg_bases ) = split( /\t/, $line );
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][AT_Transitions][covd_bases] += $covd_at_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][AT_Transversions][covd_bases] += $covd_at_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][CG_Transitions][covd_bases] += $covd_cg_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][CG_Transversions][covd_bases] += $covd_cg_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][CpG_Transitions][covd_bases] += $covd_cpg_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][CpG_Transversions][covd_bases] += $covd_cpg_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][Indels][covd_bases] += $covd_bases;
            $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][Truncations][covd_bases] += $covd_bases if( $separate_truncations );
        }
        $sampleCovgFh->close;
    }

    # Run "joinx ref-stats" to classify SNVs as being at AT, CG, or CpG sites in the reference
    print STDERR "Running 'joinx1.7 ref-stats' to read reference FASTA and identify SNVs at AT, CG, CpG sites\n";
    my $maf_bed = Genome::Sys->create_temp_file_path();
    my $mafBedFh = IO::File->new( $maf_bed, ">" ) or die "Temporary file could not be created. $!";
    my $mafFh = IO::File->new( $maf_file ) or die "Couldn't open $maf_file. $!";
    while( my $line = $mafFh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $chr, $start, $stop, $mutation_type, $ref, $var1, $var2 ) = @cols[4..6,9..12];
        if( $mutation_type =~ m/^(SNP|DNP|ONP|TNP)$/ ) {
            $mafBedFh->print( "$chr\t" . ( $start - 2 ) . "\t" . ( $start + 1 ) . "\n" );
        }
    }
    $mafFh->close;
    $mafBedFh->close;
    my $refstats_file = Genome::Sys->create_temp_file_path();
    `joinx1.7 ref-stats --ref-bases --bed $maf_bed --fasta $ref_seq --output $refstats_file`;

    # Parse through the ref-stats output and load it into hashes for quick lookup later
    my ( %ref_base, %cpg_site );
    my $refStatsFh = IO::File->new( $refstats_file ) or die "Couldn't open $refstats_file. $!";
    while( my $line = $refStatsFh->getline ) {
        next if( $line =~ m/^#/ );
        chomp $line;
        my ( $chr, undef, $pos, undef, undef, undef, $ref ) = split( /\t/, $line );
        my $locus = "$chr\t" . ( $pos - 1 );
        $ref_base{$locus} = substr( $ref, 1, 1 );
        $cpg_site{$locus} = 1 if( $ref =~ m/CG/i );
    }
    $refStatsFh->close;

    # Create a hash to help classify SNVs
    my %classify;
    $classify{$_} = AT_Transitions foreach( qw( AG TC ));
    $classify{$_} = AT_Transversions foreach( qw( AC AT TA TG ));
    $classify{$_} = CG_Transitions foreach( qw( CT GA ));
    $classify{$_} = CG_Transversions foreach( qw( CA CG GC GT ));

    # Parse through the MAF file and categorize each somatic mutation
    print STDERR "Parsing MAF file to classify mutations\n";
    my %skip_cnts;
    $mafFh = IO::File->new( $maf_file ) or die "Couldn't open $maf_file. $!";
    while( my $line = $mafFh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref, $var1, $var2, $sample ) = @cols[0,4..6,8..12,15];
        # Skip mutations in samples that are not in the provided bam list
        unless( defined $sample_idx{$sample} ) {
            $skip_cnts{"belong to unrecognized samples"}++;
            print STDERR "Skipping unrecognized sample ($sample not in BAM list): $gene, $chr:$start-$stop\n" if( $show_skipped );
            next;
        }

        # If the mutation classification is odd, quit with error
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            print STDERR "Unrecognized Variant_Classification \"$mutation_class\" in MAF file: $gene, $chr:$start-$stop\n";
            print STDERR "Please use TCGA MAF Specification v2.3.\n";
            return undef;
        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
          ( $skip_silent && $mutation_class =~ m/^Silent$/ )) {
            $skip_cnts{"are classified as $mutation_class"}++;
            print STDERR "Skipping $mutation_class mutation: $gene, $chr:$start-$stop\n" if( $show_skipped );
            next;
        }

        # If the mutation type is odd, quit with error
        if( $mutation_type !~ m/^(SNP|DNP|TNP|ONP|INS|DEL|Consolidated)$/ ) {
            print STDERR "Unrecognized Variant_Type \"$mutation_type\" in MAF file: $gene, $chr:$start-$stop\n";
            print STDERR "Please use TCGA MAF Specification v2.3.\n";
            return undef;
        }

        # Skip mutations that were consolidated into others (E.g. SNP consolidated into a TNP)
        if( $mutation_type =~ m/^Consolidated$/ ) {
            $skip_cnts{"are consolidated into another"}++;
            print STDERR "Skipping consolidated mutation: $gene, $chr:$start-$stop\n" if( $show_skipped );
            next;
        }

        # Skip mutations that fall completely outside any of the provided regions of interest
        if( $self->count_bits( $roi_bitmask->{$chr}, $start, $stop ) == 0 ) {
            $skip_cnts{"are outside any ROIs"}++;
            print STDERR "Skipping mutation that falls outside ROIs: $gene, $chr:$start-$stop\n" if( $show_skipped );
            next;
        }

        # Skip mutations whose gene names don't match any of those in the ROI list
        unless( defined $gene_idx{$gene} ) {
            $skip_cnts{"have unrecognized gene names"}++;
            print STDERR "Skipping unrecognized gene name (not in ROI file): $gene, $chr:$start-$stop\n" if( $show_skipped );
            next;
        }

        my $class = '';
        # Check if the mutation is the truncating type, if the user wanted a separate category of those
        if( $separate_truncations && $mutation_class =~ m/^(Nonsense_Mutation|Splice_Site|Frame_Shift_Del|Frame_Shift_Ins)/ ) {
            $class = Truncations;
        }
        # Classify the mutation as AT/CG/CpG Transition, AT/CG/CpG Transversion
        elsif( $mutation_type =~ m/^(SNP|DNP|ONP|TNP)$/ ) {
            # ::TBD:: For DNPs and TNPs, we use only the first base for mutation classification
            $ref = substr( $ref, 0, 1 );
            $var1 = substr( $var1, 0, 1 );
            $var2 = substr( $var2, 0, 1 );

            # If the alleles are anything but A, C, G, or T then quit with error
            if( $ref !~ m/[ACGT]/ || $var1 !~ m/[ACGT]/ || $var2 !~ m/[ACGT]/ ) {
                print STDERR "Unrecognized allele in column Reference_Allele, Tumor_Seq_Allele1, or Tumor_Seq_Allele2: $gene, $chr:$start-$stop\n";
                print STDERR "Please use TCGA MAF Specification v2.3.\n";
                return undef;
            }

            # Use the classify hash to find whether this SNV is an AT/CG Transition/Transversion
            $class = $classify{ "$ref$var1" } if( defined $classify{ "$ref$var1" } );
            $class = $classify{ "$ref$var2" } if( defined $classify{ "$ref$var2" } );

            # Check if the ref base in the MAF matched what we fetched from the ref-seq
            my $locus = "$chr\t$start";
            if( defined $ref_base{$locus} && $ref_base{$locus} ne $ref ) {
                print STDERR "Reference allele $ref for $gene variant at $chr:$start-$stop is " . $ref_base{$locus} . " in the FASTA. Using it anyway.\n";
            }

            # Check if a C or G reference allele belongs to a CpG pair in the refseq
            if(( $ref eq 'C' || $ref eq 'G' ) && defined $cpg_site{$locus} ) {
                $class = (( $class == CG_Transitions ) ? CpG_Transitions : CpG_Transversions );
            }
        }
        # Classify it as an indel (excludes splice-site and frame-shift if user wanted truncations separately)
        elsif( $mutation_type =~ m/^(INS|DEL)$/ ) {
            $class = Indels;
        }

        # The user's gene exclusion list affects only the overall BMR calculations
        $sample_mr[$sample_idx{$sample}][$class][mutations]++ unless( defined $ignored_genes{$gene} );
        $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$class][mutations]++;
    }
    $mafFh->close;

    # Display statistics related to parsing the MAF
    print STDERR "Finished Parsing the MAF file to classify mutations\n";
    foreach my $skip_type ( sort {$skip_cnts{$b} <=> $skip_cnts{$a}} keys %skip_cnts ) {
        print STDERR "Skipped " . $skip_cnts{$skip_type} . " mutation(s) that $skip_type\n" if( defined $skip_cnts{$skip_type} );
    }

    # If the user wants, merge together concurrent mutations of a gene in the same sample
    if( $merge_concurrent_muts ) {
        foreach my $sample ( @all_sample_names ) {
            foreach my $gene ( @all_gene_names ) {
                next unless( defined $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}] );
                my $num_muts = 0;
                $num_muts += $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$_][mutations] foreach( @mut_classes );
                if( $num_muts > 1 ) {
                    foreach my $class ( @mut_classes ) {
                        my $muts_in_class = $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$class][mutations]; # Num of muts of gene in this class
                        $sample_mr[$sample_idx{$sample}][$class][mutations] -= $muts_in_class; # Take it out of the sample total
                        $muts_in_class /= $num_muts; # Turn it into a fraction of the total number of muts in this gene
                        $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$class][mutations] = $muts_in_class; # Use the fraction as the num muts of gene in this class
                        $sample_mr[$sample_idx{$sample}][$class][mutations] += $muts_in_class; # Add the same fraction to the sample total
                    }
                }
            }
        }
    }

    # Calculate per-sample BMRs, and also subtract out covered bases in genes the user wants ignored
    foreach my $sample ( @all_sample_names ) {
        my $tot_muts = 0;
        foreach my $class ( @mut_classes ) {
            # Subtract the covered bases in this class that belong to the genes to be ignored
            # ::TBD:: Some of these bases may also belong to another gene (on the other strand maybe?), and those should not be subtracted
            foreach my $ignored_gene ( keys %ignored_genes ) {
                $sample_mr[$sample_idx{$sample}][$class][covd_bases] -= $gene_mr[$sample_idx{$sample}][$gene_idx{$ignored_gene}][$class][covd_bases] if( defined $gene_mr[$sample_idx{$sample}][$gene_idx{$ignored_gene}] );
            }
            $tot_muts += $sample_mr[$sample_idx{$sample}][$class][mutations];
        }
        $sample_mr[$sample_idx{$sample}][Overall][bmr] = $tot_muts / $sample_mr[$sample_idx{$sample}][Indels][covd_bases];
    }

    # Cluster samples into bmr-groups using k-means clustering
    my @sample_bmrs = map { $sample_mr[$sample_idx{$_}][Overall][bmr] } @all_sample_names;
    my @bmr_clusters = k_means( $bmr_groups, \@sample_bmrs );

    # Calculate overall BMRs for each cluster of samples, and print them to file
    my %cluster_bmr; # Stores per cluster categorized BMR
    my $totBmrFh = IO::File->new( $overall_bmr_file, ">" ) or die "Couldn't open $overall_bmr_file. $!";
    $totBmrFh->print( "#User-specified genes skipped in these calculations: $genes_to_ignore\n" ) if( defined $genes_to_ignore );
    my ( $covered_bases_sum, $mutations_sum ) = ( 0, 0 );
    for( my $i = 0; $i < scalar( @bmr_clusters ); ++$i ) {
        my @samples_in_cluster = map { $all_sample_names[$_] } @{$bmr_clusters[$i]};
        unless( $bmr_groups == 1 ) {
            $totBmrFh->print( "#BMR sub-group ", $i + 1, " (", scalar( @{$bmr_clusters[$i]} ), " samples)\n" );
            $totBmrFh->print( "#Samples: ", join( ",", @samples_in_cluster ), "\n" );
        }
        $totBmrFh->print( "#Mutation_Class\tCovered_Bases\tMutations\tOverall_BMR\n" );

        my ( $tot_covd_bases, $tot_muts ) = ( 0, 0 );
        foreach my $class ( @mut_classes ) {
            my ( $covd_bases, $mutations ) = ( 0, 0 );
            foreach my $sample ( @samples_in_cluster ) {
                $covd_bases += $sample_mr[$sample_idx{$sample}][$class][covd_bases];
                $mutations += $sample_mr[$sample_idx{$sample}][$class][mutations];
            }
            $tot_covd_bases = $covd_bases if( $class == Indels ); # Save this to calculate overall BMR below
            # Calculate overall BMR for this mutation class and print it to file
            $cluster_bmr{$i}[$class][bmr] = ( $covd_bases == 0 ? 0 : ( $mutations / $covd_bases ));
            $totBmrFh->print( join( "\t", $mut_class_names[$class], $covd_bases, $mutations, $cluster_bmr{$i}[$class][bmr] ), "\n" );
            $tot_muts += $mutations;
        }
        $totBmrFh->print( join( "\t", "Overall_BMR", $tot_covd_bases, $tot_muts, $tot_muts / $tot_covd_bases ), "\n\n" );
        $covered_bases_sum += $tot_covd_bases;
        $mutations_sum += $tot_muts;
    }
    $totBmrFh->close;
    $self->bmr_output( $mutations_sum / $covered_bases_sum );

    # Print out a file containing per-gene mutation counts and covered bases for use by "music smg"
    my $geneBmrFh = IO::File->new( $gene_mr_file, ">" ) or die "Couldn't open $gene_mr_file. $!";
    $geneBmrFh->print( "#Gene\tMutation_Class\tCovered_Bases\tMutations\tBMR\n" );
    foreach my $gene ( sort @all_gene_names ) {
        my ( $tot_covd_bases, $tot_muts ) = ( 0, 0 );
        for( my $i = 0; $i < scalar( @bmr_clusters ); ++$i ) {
            my @samples_in_cluster = map { $all_sample_names[$_] } @{$bmr_clusters[$i]};
            foreach my $class ( @mut_classes ) {
                my ( $covd_bases, $mutations ) = ( 0, 0 );
                foreach my $sample( @samples_in_cluster ) {
                    if( defined $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}] ) {
                        $covd_bases += $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$class][covd_bases];
                        $mutations += $gene_mr[$sample_idx{$sample}][$gene_idx{$gene}][$class][mutations];
                    }
                }
                my $rename_class = $mut_class_names[$class];
                $rename_class = ( $rename_class . "_SubGroup" . ( $i + 1 )) if( $bmr_groups > 1 );
                $geneBmrFh->print( join( "\t", $gene, $rename_class, $covd_bases, $mutations, $cluster_bmr{$i}[$class][bmr] ), "\n" );
                $tot_muts += $mutations;
                $tot_covd_bases += $covd_bases if( $class == Indels );
            }
        }
        $geneBmrFh->print( join( "\t", $gene, "Overall", $tot_covd_bases, $tot_muts, $self->bmr_output ), "\n" );
    }
    $geneBmrFh->close;

    return 1;
}

# Creates an empty whole genome bitmask based on the given reference sequence index
sub create_empty_genome_bitmask {
    my ( $self, $ref_seq_idx_file ) = @_;
    my %genome;
    my $refFh = IO::File->new( $ref_seq_idx_file ) or die "Couldn't open $ref_seq_idx_file. $!";
    while( my $line = $refFh->getline ) {
        my ( $chr, $length ) = split( /\t/, $line );
        $genome{$chr} = Bit::Vector->new( $length + 1 ); # Adding a base for 1-based coordinates
    }
    $refFh->close;
    return \%genome;
}

# Counts the number of bits that are set in the given region of a Bit:Vector
sub count_bits {
    my ( $self, $vector, $start, $stop ) = @_;
    my $count = 0;
    for my $pos ( $start..$stop ) {
        ++$count if( $vector->bit_test( $pos ));
    }
    return $count;
}

# Given a list of numerical values, returns k clusters based on k-means clustering
sub k_means {
    my ( $k, $list_ref ) = @_;
    my @vals = @{$list_ref};
    my $num_vals = scalar( @vals );

    # Start with the first k values as the centroids
    my @centroids = @vals[0..($k-1)];
    my @prev_centroids = map { 0 } @centroids;
    my @groups = ();

    my $diff_means = 1; # Arbitrary non-zero value
    # Repeat until the difference between these centroids and the previous ones, converges to zero
    while( $diff_means > 0 ) {
        @groups = ();
        # Group values into clusters based on closest centroid
        for( my $i = 0; $i < $num_vals; ++$i ) {
            my @distances = map { abs( $vals[$i] - $_ ) } @centroids;
            my $closest = min( @distances );
            for( my $j = 0; $j < $k; ++$j ) {
                if( $distances[$j] == $closest ) { push( @{$groups[$j]}, $i ); last; }
            }
        }

        # Calculate means to be the new centroids, and the sum of differences
        $diff_means = 0;
        for( my $i = 0; $i < $k; ++$i ) {
            $centroids[$i] = sum( map {$vals[$_]} @{$groups[$i]} );
            $centroids[$i] /= scalar( @{$groups[$i]} );
            $diff_means += abs( $centroids[$i] - $prev_centroids[$i] );
        }

        # Save the current centroids for comparisons with those in the next iteration
        @prev_centroids = @centroids;
    }
    return @groups;
}

1;
