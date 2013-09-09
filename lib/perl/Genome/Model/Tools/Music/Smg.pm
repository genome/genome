package Genome::Model::Tools::Music::Smg;

use warnings;
use strict;
use Genome;
use IO::File;
use Statistics::Descriptive;
use Carp;
use POSIX qw( WIFEXITED );

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Smg {
    is => 'Command::V2',
    has_input => [
        gene_mr_file => { is => 'Text', doc => "File with per-gene mutation rates (Created using \"music bmr calc-bmr\")" },
        output_file => { is => 'Text', is_output => 1, doc => "Output file that will list significantly mutated genes and their p-values" },
    ],
    has_optional_input => [
        max_fdr => { is => 'Number', default => 0.20, doc => "The maximum allowed false discovery rate for a gene to be considered an SMG" },
        skip_low_mr_genes => { is => 'Boolean', default => 1, doc => "Skip testing genes with MRs lower than the background MR" },
        downsample_large_genes => { is => 'Boolean', default => 0, doc => "Downscale #bps in large genes, and the #muts proportionally" },
        bmr_modifier_file => { is => 'Text', doc => "Tab delimited multipliers per gene that modify BMR before testing [gene_name bmr_modifier]" },
        processors => { is => 'Integer', default => 1, doc => "Number of cores to use (requires 'foreach' and 'doMC' R packages)" },
    ],
    doc => "Identify significantly mutated genes."
};

sub help_synopsis {
    return <<HELP
 ... music smg \\
      --gene-mr-file output_dir/gene_mrs \\
      --output-file output_dir/smgs

(A "gene-mr-file" can be generated using the tool "music bmr calc-bmr".)
HELP
}

sub help_detail {
    return <<HELP
This script runs R-based statistical tools to identify Significantly Mutated Genes (SMGs), when
given per-gene mutation rates categorized by mutation type, and the overall background mutation
rates (BMRs) for each of those categories (gene_mr_file, created using "music bmr calc-bmr").

P-values and false discovery rates (FDRs) for each gene in gene_mr_file is calculated using three
tests: Fisher\'s Combined P-value test (FCPT), Likelihood Ratio test (LRT), and the Convolution
test (CT). For a gene, if its FDR for at least 2 of these tests is <= max_fdr, it will be output
as an SMG. Another output file with prefix "_detailed" will have p-values and FDRs for all genes.
HELP
}

sub _additional_help_sections {
    return (
    "ARGUMENTS",
<<EOS

=over 4

=item --skip-low-mr-genes

=over 8

=item Genes with consistently lower MRs than the BMRs across mutation categories, may show up in
  the results as an SMG (by CT or LRT). If such genes are not of interest, they may be assigned a
  p-value of 1. This should also speed things up. Genes with higher Indel or Truncation rates than
  the background will not be skipped even if the gene's overall MR is lower than the BMR. If
  bmr-modifiers are applied, this step uses the modified BMRs instead.

=back

=item --downsample-large-genes

=over 8

=item Large genes have a higher power for detecting mutational significance than smaller genes. So
  a workaround is to downsample their covd bps and mutation count, keeing MR unchanged. A gene with
  more covd bps than Q3+(Q3-Q1)*4.5 is considered an outlier, where Q1 & Q3 are the lower and upper
  quartiles respectively, of #covd_bps in the gene-mr-file. For each outlier, the #covd bps in each
  mutation category is brought closer to Q3+(Q3-Q1)*4.5, and the #muts reduced proportionally, to
  keep the gene's MR mostly unchanged.

=back

=item --bmr-modifier-file

=over 8

=item The user can provide a BMR modifier for each gene in the ROI file, which is a multiplier for
  the categorized background mutation rates, before testing them against the gene's categorized
  mutation rates. Such a file can be used to correct for regional or systematic bias in mutation
  rates across the genome that may be correlated to CpG deamination or DNA repair processes like
  transcription-coupled repair or mismatch repair. Mutation rates have also been associated with
  DNA replication timing, where higher mutation rates are seen in late replicating regions. Note
  that the same per-gene multiplier is used on each mutation category of BMR. Any genes from the
  ROI file that are not in the BMR modifier file will be tested against unmodified overall BMRs
  per mutation category. BMR modifiers of <=0 are not permitted, because that's just silly.

=back

=back

=cut

EOS
  );
}

sub _doc_authors {
    return <<EOS
 Qunyuan Zhang, Ph.D.
 Cyriac Kandoth, Ph.D.
 Nathan D. Dees, Ph.D.
EOS
}

sub execute {
    my $self = shift;
    my $gene_mr_file = $self->gene_mr_file;
    my $output_file = $self->output_file;
    my $output_file_detailed = $output_file . "_detailed";
    my $max_fdr = $self->max_fdr;
    my $skip_low_mr_genes = $self->skip_low_mr_genes;
    my $downsample_large_genes = $self->downsample_large_genes;
    my $bmr_modifier_file = $self->bmr_modifier_file;
    my $processors = $self->processors;

    # Check on all the input data before starting work
    print STDERR "Gene mutation rate file not found or is empty: $gene_mr_file\n" unless( -s $gene_mr_file );
    print STDERR "BMR modifier file not found or is empty: $bmr_modifier_file\n" unless( !defined $bmr_modifier_file || -s $bmr_modifier_file );
    return undef unless( -s $gene_mr_file && ( !defined $bmr_modifier_file || -s $bmr_modifier_file ));

    # Collect per-gene mutation rates for reporting in results later
    my ( %gene_muts, %gene_bps, %mut_categ_hash );
    my $inMrFh = IO::File->new( $gene_mr_file ) or die "Couldn't open $gene_mr_file. $!\n";
    while( my $line = $inMrFh->getline ) {
        next if( $line =~ m/^#/ );
        my ( $gene, $type, $covd_bps, $mut_cnt, undef ) = split( /\t/, $line );

        # Warn user about cases where there could be fewer covered bps than mutations detected
        if( $type ne "Overall" and $covd_bps < $mut_cnt ) {
            warn "#More $type seen in $gene than there are bps with sufficient coverage!\n";
        }

        if( $type eq "Overall" ) {
            $gene_muts{$gene}{Overall} = $mut_cnt;
            $gene_bps{$gene} = $covd_bps;
        }
        elsif( $type =~ m/^Truncations/ ) {
            $gene_muts{$gene}{Truncations} += $mut_cnt;
            $mut_categ_hash{Truncations} = 1;
        }
        elsif( $type =~ m/^Indels/ ) {
            $gene_muts{$gene}{Indels} += $mut_cnt;
            $mut_categ_hash{Indels} = 1;
        }
        elsif( $type =~ m/^(AT_|CG_|CpG_)(Transitions|Transversions)/ ) {
            $gene_muts{$gene}{SNVs} += $mut_cnt;
            $mut_categ_hash{SNVs} = 1;
        }
        else {
            die "Unrecognized mutation category in gene-mr-file. $!\n";
        }
    }
    $inMrFh->close;
    my @mut_categs = sort keys %mut_categ_hash;

    # If requested, downsample large genes and create another gene_mr_file
    if( $downsample_large_genes ) {

        # Collect covd_bps across genes per mutation category
        my ( %gene_categ_bps );
        my $inMrFh = IO::File->new( $gene_mr_file ) or die "Couldn't open $gene_mr_file. $!\n";
        while( my $line = $inMrFh->getline ) {
            next if( $line =~ m/^#/ );
            my ( $gene, $type, $covd_bps, undef, undef ) = split( /\t/, $line );
            push( @{$gene_categ_bps{$type}}, $covd_bps );
        }
        $inMrFh->close;

        # For each mutation category, find the Q1 and Q3 quartiles, and the median
        my ( %q1, %q2, %q3 );
        foreach my $type ( keys %gene_categ_bps ) {
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data( @{$gene_categ_bps{$type}} );
            $q1{$type} = $stat->quantile( 1 );
            $q2{$type} = $stat->quantile( 2 );
            $q3{$type} = $stat->quantile( 3 );
        }

        # For each mutation category, calculate the outlier cutoff
        my %cutoff;
        print "#Selecting large genes to downsample based on these cutoffs...\n";
        print "#Mutation_Type\tMedian\tIQR\tCutoff\n";
        foreach my $type ( sort keys %gene_categ_bps ) {
            unless( $type eq "Overall" ) {
                my $IQR = $q3{$type} - $q1{$type}; # Interquartile range
                $cutoff{$type} = sprintf( "%.0f", $q3{$type} + ( $IQR * 4.5 ));
                print "$type\t" . $q2{$type} . "\t$IQR\t" . $cutoff{$type} . "\n";
            }
        }

        # Create a new gene_mr_file with large genes downsampled in #covd_bps and #muts
        print "\n#Changed lines in $gene_mr_file as follows...\n";
        my $new_gene_mr_file = $gene_mr_file . "_downsampled_large_genes";
        my $outMrFh = IO::File->new( $new_gene_mr_file, ">" ) or die "Couldn't create $new_gene_mr_file. $!\n";
        $inMrFh = IO::File->new( $gene_mr_file ) or die "Couldn't open $gene_mr_file. $!\n";
        while( my $line = $inMrFh->getline ) {
            if( $line =~ m/^#/ ) {
                $outMrFh->print( $line );
                next;
            }
            my ( $gene, $type, $covd_bps, $mut_cnt, $bmr ) = split( /\t/, $line );

            #If the #covd bps in this line is an outlier for its category, then downsample it
            #The category "Overall" is only for reporting purposes, so we can leave it unchanged
            if( $type ne "Overall" and $covd_bps > $cutoff{$type} ) {
                # If this line has zero mutations, then simply set $covd_bps to the cutoff
                if( $mut_cnt == 0 ) {
                    $covd_bps = $cutoff{$type};
                }
                # Otherwise, iteratively reduce $mut_cnt, no lower than 1 mutation, bringing
                # $covd_bps as close to the cutoff as possible, but keeping the MR unchanged
                else {
                    for( my $i = $mut_cnt - 1; $i >= 0; --$i ) {
                        if(( $covd_bps * $i / $mut_cnt ) <= $cutoff{$type} || $i == 0 ) {
                            ++$i;
                            $covd_bps = sprintf( "%.0f", $covd_bps * $i / $mut_cnt );
                            $mut_cnt = $i;
                            last;
                        }
                    }
                }
                # Print out the changed #bps and #muts, for users' benefit
                if( $line ne "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr" ) {
                    print "< $line";
                    print "> $gene\t$type\t$covd_bps\t$mut_cnt\t$bmr";
                }
            }

            $outMrFh->print( "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr" );
        }
        $outMrFh->close;
        $inMrFh->close;

        $gene_mr_file = $new_gene_mr_file;
    }

    # If BMR modifiers were provided, then create another gene_mr_file with modified BMRs
    if( defined $bmr_modifier_file ) {
        my $inBmrModFh = IO::File->new( $bmr_modifier_file ) or die "Couldn't open $bmr_modifier_file. $!\n";
        my %bmr_modifier = ();
        while( my $line = $inBmrModFh->getline ) {
            next if( $line =~ m/^#/ );
            chomp( $line );
            my ( $gene, $modifier ) = split( /\t/, $line );
            ( $modifier > 0 ) or die "$modifier is an invalid bmr-modifier. Please fix values in $bmr_modifier_file.\n";
            $bmr_modifier{$gene} = $modifier;
        }
        $inBmrModFh->close;

        my $tmp_gene_mr_file = Genome::Sys->create_temp_file_path();
        my $tmpFh = IO::File->new( $tmp_gene_mr_file, ">" ) or die "Temporary file could not be created. $!";
        my $inMrFh = IO::File->new( $gene_mr_file ) or die "Couldn't open $gene_mr_file. $!\n";
        while( my $line = $inMrFh->getline ) {
            if( $line =~ m/^#/ ) {
                $tmpFh->print( $line );
                next;
            }
            chomp( $line );
            my ( $gene, $type, $covd_bps, $mut_cnt, $bmr ) = split( /\t/, $line );
            $bmr = $bmr * $bmr_modifier{$gene} if( defined $bmr_modifier{$gene} );
            $tmpFh->print( "$gene\t$type\t$covd_bps\t$mut_cnt\t$bmr\n" );
        }
        $tmpFh->close;
        $inMrFh->close;

        $gene_mr_file = $tmp_gene_mr_file;
    }

    # Create a temporary intermediate file to hold the p-values
    my $pval_file = Genome::Sys->create_temp_file_path();
    my $pvalFh = IO::File->new( $pval_file, ">" ) or die "Temporary file could not be created. $!";

    # Call R for Fisher combined test, Likelihood ratio test, and convolution test on each gene
    my $smg_cmd = "R --slave --args < " . __FILE__ . ".R $gene_mr_file $pval_file smg_test $processors $skip_low_mr_genes";
    WIFEXITED( system $smg_cmd ) or croak "Couldn't run: $smg_cmd ($?)";

    # Call R for calculating FDR on the p-values calculated in the SMG test
    my $fdr_cmd = "R --slave --args < " . __FILE__ . ".R $pval_file $output_file_detailed calc_fdr $processors $skip_low_mr_genes";
    WIFEXITED( system $fdr_cmd ) or croak "Couldn't run: $fdr_cmd ($?)";

    # Parse the R output to identify the SMGs (significant by at least 2 of 3 tests)
    my $smgFh = IO::File->new( $output_file_detailed ) or die "Couldn't open $output_file_detailed. $!\n";
    my ( @newLines, @smgLines );
    my $header = "#Gene\t" . join( "\t", @mut_categs );
    $header .= "\tTot Muts\tCovd Bps\tMuts pMbp\tP-value FCPT\tP-value LRT\tP-value CT\tFDR FCPT\tFDR LRT\tFDR CT\n";
    while( my $line = $smgFh->getline ) {
        chomp( $line );
        if( $line =~ m/^Gene\tp.fisher\tp.lr\tp.convol\tfdr.fisher\tfdr.lr\tfdr.convol$/ ) {
            push( @newLines, $header );
            push( @smgLines, $header );
        }
        else {
            my ( $gene, @pq_vals ) = split( /\t/, $line );
            my ( $p_fcpt, $p_lrt, $p_ct, $q_fcpt, $q_lrt, $q_ct ) = @pq_vals;
            my @mut_cnts;
            foreach( @mut_categs ) {
                # If a mutation count is a fraction, round down the digits after the decimal point
                push( @mut_cnts, (( $gene_muts{$gene}{$_} =~ m/\./ ) ? sprintf( "%.2f", $gene_muts{$gene}{$_} ) : $gene_muts{$gene}{$_} ));
            }
            my $mut_per_mbp = ( $gene_bps{$gene} ? sprintf( "%.2f", ( $gene_muts{$gene}{Overall} / $gene_bps{$gene} * 1000000 )) : 0 );
            push( @newLines, join( "\t", $gene, @mut_cnts, $gene_muts{$gene}{Overall}, $gene_bps{$gene}, $mut_per_mbp, @pq_vals ) . "\n" );

            # If the FDR of at least two of these tests is less than the maximum allowed, we consider it an SMG
            if(( $q_fcpt <= $max_fdr && $q_lrt <= $max_fdr ) || ( $q_fcpt <= $max_fdr && $q_ct <= $max_fdr ) ||
               ( $q_lrt <= $max_fdr && $q_ct <= $max_fdr )) {
                push( @smgLines, join( "\t", $gene, @mut_cnts, $gene_muts{$gene}{Overall}, $gene_bps{$gene}, $mut_per_mbp, @pq_vals ) . "\n" );
            }
        }
    }
    $smgFh->close;

    # Add per-gene SNV and Indel counts to the detailed R output, and make the header friendlier
    my $outDetFh = IO::File->new( $output_file_detailed, ">" ) or die "Couldn't open $output_file_detailed. $!\n";
    $outDetFh->print( @newLines );
    $outDetFh->close;

    # Do the same for only the genes that we consider SMGs
    my $outFh = IO::File->new( $output_file, ">" ) or die "Couldn't open $output_file. $!\n";
    $outFh->print( @smgLines );
    $outFh->close;

    return 1;
}

1;
