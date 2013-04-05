package Genome::Model::Tools::Music::Data::ParseCosmic;

use warnings;
use strict;
use IO::File;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Data::ParseCosmic {
    is => 'Command::V2',
    has_input => [
        input_complete_file => { is => 'String', doc => "CosmicCompleteExport downloaded from cosmic" },
        input_ins_mut_file => { is => 'String', doc => "CosmicInsMutExport downloaded from cosmic" },
        output_file => { is => 'String', doc => "Cleaned up and standardized variant list with sample IDs" },
    ],
    has_optional_input => [
        hg18_fasta => { is => 'String', doc => "Reference FASTA file for HG18 (Build36)",
                        example_values => ["/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa"] },
        hg19_fasta => { is => 'String', doc => "Reference FASTA file for HG19 (Build37)",
                        example_values => ["/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa"] },
        max_indel_length => { is => 'Integer', doc => "Skip indels longer than this many bps", default => 100 },
        genome_wide_only => { is => 'Boolean', doc => "Pull only variants from genome-wide or exome-wide screens", default => 0 },
    ],
    doc => "Parses and standardizes nucleotide-level variants downloaded from COSMIC",
};

sub help_detail {
    return <<HELP;
The two expected input files are downloadable from Sanger's Public FTP at:
ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export

"input_complete_file" should be the complete mutation list from COSMIC:
CosmicCompleteExport_v62_291112.tsv

"input_ins_mut_file" should be the separate file with nucleotide sequences for somatic insertions:
CosmicInsMutExport_v62_291112.tsv

Shown above are the filenames from COSMIC release v62. You can find other versions of these files
on Sanger's FTP. This tool should work fine with most of them.

Note that this tool will skip the following based on data available in COSMIC:

=over

 - variants that don't have both nucleotide change and genomic loci available
 - variants with reference alleles that don't match the provided hg18/hg19 reference sequence
 - variants annotated as germline, or from samples annotated as normals
 - variants annotated as belonging to metastases, adjacent tissue, recurrent tumors, etc.
 - insertions whose sequence don't match the separate cosmic file for somatic insertions
 - deletions whose length don't match the genomic loci
 - TNPs, ONPs and complex indels (DNPs are not skipped)
HELP
}

#########################################################################

sub execute {
    my $self = shift;
    my $input_complete_file = $self->input_complete_file;
    my $input_ins_mut_file = $self->input_ins_mut_file;
    my $hg18_fasta = $self->hg18_fasta;
    my $hg19_fasta = $self->hg19_fasta;
    my $max_indel_length = $self->max_indel_length;
    my $genome_wide_only = $self->genome_wide_only;
    my $output_file = $self->output_file;

    # Some hashes for quick lookup
    my %valid_chrom = map {($_,1)} ( 1..22, qw( X Y MT ));
    my %complement = qw( A T T A C G G C N N );

    # Load up the sequences of insertions, which COSMIC happens to store separately for whatever reason
    my %ins_nucs = map{chomp; split(/\t/)} grep {m/^\w+\t[ACGTacgt]+$/} `cut -f 9,11 $input_ins_mut_file`;

    # From the main cosmic file, parse out genomic loci and ref/var for SNVs and indels
    my %vars = ();
    my %tissue_site_counts = ();
    my $cosmicFh = IO::File->new( $input_complete_file );
    warn "Parsing COSMIC DB for usable SNVs and indels...\n";
    while( my $line = $cosmicFh->getline ) {
        next if( $line =~ m/^Gene name/ );
        chomp( $line );
        my @col = split( /\t/, $line );
        my ( $samp_name, $samp_id, $samp_site, $samp_hist, $samp_subtype, $wgs, $origin ) = @col[3,4,6,8,9,10,23];
        my ( $mut_id, $mut_nuc, $mut_aa, $mut_type, $zygos, $hg18_locus, $hg18_strand, $hg19_locus, $hg19_strand, $status ) = @col[11..20];

        # Skip variants of samples that were not "Genome-wide screened" (if the flag is turned on)
        next if( $genome_wide_only and $wgs ne 'y' );

        # Skip variants belonging to recurrent tumors, metastases, adjacent, etc. - We don't want to report the same variant from the same patient more than once
        next unless( $origin =~ m/^(primary|secondary|surgery)/ );

        # Skip known germline variants or variants from samples annotated as normals
        next if( $status =~ m/germline/ or $samp_hist =~ m/NORMAL/ or $samp_subtype =~ m/normal/);

        # Skip mutation types that we cannot deal with
        next if( $mut_type =~ m/^Whole gene deletion/ );

        # Skip variants that don't have both nucleotide change and loci available
        next unless( $mut_nuc and ( $hg18_locus or $hg19_locus ));

        # Pull b36 or b37 loci of the variant giving preference to b36, which we will liftOver later
        my $build = ( $hg18_locus ? "b36" : "b37" );
        my $locus = ( $hg18_locus ? $hg18_locus : $hg19_locus );
        my ( $chr, $start, $stop ) = $locus =~ m/^(\w+):(\d+)-(\d+)$/;
        die "Cannot identify locus for:\n$line\n" unless( $chr and $start and $stop );
        $chr =~ s/^(23|x)$/X/; $chr =~ s/^(24|y)$/Y/; $chr =~ s/^(25|M|m|mt)$/MT/;
        die "Invalid chrom name in:\n$line\n" unless( $valid_chrom{$chr} );

        # Fetch the reference sequence at this locus to do some QC
        my $ref_seq = ( $hg18_locus ? $hg18_fasta : $hg19_fasta );
        my $decr_start = $start - 1;
        my $fetched_ref = `echo "$chr\t$decr_start\t$stop" | joinx ref-stats --ref-bases - $ref_seq | grep -v ^# | cut -f 7`;
        chomp( $fetched_ref );

        # Try to find out what kind of variant this is, and pull ref/var if possible
        my ( $ref, $var, @tmp );
        if( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*([ACGTacgt]+)>([ACGTacgt]+)$/ ) { # SNV
            ( $ref, $var ) = ( uc( $tmp[0] ), uc( $tmp[1] ));
        }
        elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*ins([ACGTacgt]+)$/ ) { # insertion
            ( $ref, $var ) = ( "-", uc( $tmp[0] ));
            my $ins_nuc = uc( $ins_nucs{$mut_id} );
            if( length( $var ) > $max_indel_length ) {
                warn "Skipped: Long insertion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
                next;
            }
            elsif( $ins_nuc and $var ne $ins_nuc ) {
                warn "Skipped: Inserted sequence per COSMIC ($ins_nuc) differs from $mut_nuc at $build locus: $locus\n";
                next;
            }
        }
        elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del([ACGTacgt]+)$/ ) { # deletion
            ( $ref, $var ) = ( uc( $tmp[0] ), "-" );
            if( length( $ref ) > $max_indel_length ) {
                warn "Skipped: Long deletion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
                next;
            }
            elsif( length( $ref ) != ( $stop - $start + 1 )) {
                warn "Skipped: Length of deleted sequence in $mut_nuc doesn't match $build locus: $locus\n";
                next;
            }
        }
        elsif( @tmp = $mut_nuc =~ m/^c\.\d+[+-]*\d*_*\d*[+-]*\d*del(\d+)$/ ) { # deletion without sequence specified
            my $del_length = $tmp[0];
            if( $del_length > $max_indel_length ) {
                warn "Skipped: Long deletion >$max_indel_length bps in $mut_nuc at $build locus: $locus\n";
                next;
            }

            if( length( $fetched_ref ) != $del_length ) {
                warn "Skipped: Length of deleted sequence in $mut_nuc doesn't match $build locus: $locus\n";
                next;
            }
            ( $ref, $var ) = ( uc( $fetched_ref ), "-" );
        }

        unless( $ref and $var ) {
            warn "Skipped: Cannot identify ref/var from $mut_nuc at $build locus: $locus\n";
            next;
        }

        unless(( length($ref)==1 && length($var)==1 ) || ( length($ref)==2 && length($var)==2 ) ||
               ( $ref eq "-" && $var=~m/[ACGT]+/ ) || ( $ref=~m/[ACGT]+/ && $var eq "-" )) {
            warn "Skipped: Cannot handle $ref/$var at $build locus: $locus\n";
            next;
        }

        # Reverse complement the ref/var if the variant is from the negative strand
        my $locus_strand = ( $hg18_strand ? $hg18_strand : $hg19_strand );
        my $rc_ref = reverse map {$complement{$_}} split( //, $ref );
        if( $locus_strand eq "-" and $rc_ref eq $fetched_ref ) {
            if( $ref eq "-" && $var =~ m/[ACGT]+/ ) {
                $var = reverse map {$complement{$_}} split( //, $var );
            }
            elsif( $ref =~ m/[ACGT]+/ && $var eq "-" ) {
                $ref = $rc_ref;
            }
            else {
                $ref = $rc_ref;
                $var = reverse map {$complement{$_}} split( //, $var );
            }
        }

        # Skip SNVs and deletions with reference alleles that don't match the reference fasta
        if( $ref ne $fetched_ref and $ref ne "-" ) {
            warn "Skipped: Variant's ref allele $ref is $fetched_ref in the $build fasta at locus: $locus\n";
            next;
        }

        # If we get to this point, then we have all the information we need. Save it into a hash
        $samp_site = $samp_hist if( $samp_site eq "NS" );
        ++$tissue_site_counts{$samp_site} unless( defined $vars{$samp_id} );
        $vars{$samp_id}{$build}{"$chr\t$start\t$stop\t$ref\t$var"} = 1;
    }
    $cosmicFh->close;

    # Remove output file if it already exists
    unlink $output_file if( -e $output_file );

    # De-deduplicate variants by sample ID, and dump it into the same output file
    my ( $samp_count, $tot_mut_count ) = ( 0, 0 );
    warn "\nRunning liftOver and prepping Build37 files in a more familiar format...\n";
    foreach my $samp_id ( keys %vars ) {

        # Make 2 temporary files to use with liftover
        my $build36_file = Genome::Sys->create_temp_file_path;
        my $build37_file = Genome::Sys->create_temp_file_path;
        ( $build36_file and $build37_file ) or die "Couldn't create a temp file. $!";

        # Write any build36 variants into a 5 column WU format for use with liftOver
        my $mut_count = 0;
        if( defined $vars{$samp_id}{b36} && scalar( keys %{$vars{$samp_id}{b36}} ) > 0 ) {
            $mut_count += scalar( keys %{$vars{$samp_id}{b36}} );
            my $outFh = IO::File->new( $build36_file, ">" );
            foreach my $line ( keys %{$vars{$samp_id}{b36}} ) {
                $outFh->print( "$line\t$samp_id\n" );
            }
            $outFh->close;

            # Run liftOver using the GMT wrapper (::TODO:: use liftOver binary instead)
            my $lift_cmd = Genome::Model::Tools::LiftOver->create(
                source_file => $build36_file, destination_file => $build37_file,
                input_is_annoformat => 1, lift_direction => "hg18ToHg19",
            );
            ( $lift_cmd->execute ) or die "Failed to run 'gmt lift-over' for variants in $samp_id\n";
            $lift_cmd->delete;
        }

        # Append additional b37 loci to the liftOver'ed file, and remove any duplicates
        if( defined $vars{$samp_id}{b37} && scalar( keys %{$vars{$samp_id}{b37}} ) > 0 ) {
            $mut_count += scalar( keys %{$vars{$samp_id}{b37}} );
            my $outFh = IO::File->new( $build37_file, ">>" );
            foreach my $line ( keys %{$vars{$samp_id}{b37}} ) {
                $outFh->print( "$line\t$samp_id\n" );
            }
            $outFh->close;
        }
        print `sort -u $build37_file -o $build37_file`;

        # Sort by loci, and write to output file
        if( -s $build37_file ) {
            ++$samp_count;
            $tot_mut_count += $mut_count;
            print `joinx sort -s -i $build37_file >> $output_file`;
        }
    }

    print "\nFetched $tot_mut_count SNVs and indels across $samp_count samples from the following tissue types:\n";
    foreach my $key ( sort {$tissue_site_counts{$b} <=> $tissue_site_counts{$a}} keys %tissue_site_counts ) {
        print $tissue_site_counts{$key} . "\t$key\n";
    }

    return 1;
}

1;
