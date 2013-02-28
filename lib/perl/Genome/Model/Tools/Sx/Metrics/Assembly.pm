package Genome::Model::Tools::Sx::Metrics::Assembly;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Metrics::Assembly {
    is => 'Genome::Model::Tools::Sx::Metrics',
    has => [
        supercontigs => {
            is => 'Hash',
            is_optional => 1,
            default_value => {},
        },
        contigs => {
            is => 'Hash',
            is_optional => 1,
            default_value => {}, 
        },
        ( map { 
                $_ => {
                    #is => 'UR::Value', 
                    is_optional => 1,
                    default_value => ( ( /_percent$/ or /_success$/ or /^coverage_/ or /^core_/ ) ? 'NA' : 0 ), 
                },
            } __PACKAGE__->metric_names
        ),
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->tier_two( $self->tier_one ) if $self->tier_one and not $self->tier_two;

    return $self;
}

sub metric_names {
    return (qw/
        tier_one tier_two major_contig_threshold
        assembler_name assembler_version assembler_params assembler_kmer
        assembly_length
        coverage genome_size insert_size read_processor subject_name
        chaff_rate
        content_at content_gc content_nx
        contigs_average_length
        contigs_count
        contigs_length
        contigs_length_5k
        contigs_minor_count
        contigs_minor_length
        contigs_minor_average_length
        contigs_minor_n50_count
        contigs_minor_n50_length
        contigs_minor_read_count
        contigs_minor_read_percent
        contigs_major_average_length
        contigs_major_count
        contigs_major_length
        contigs_major_n50_count
        contigs_major_n50_length
        contigs_major_read_count
        contigs_major_read_percent
        contigs_major_percent
        contigs_maximum_length
        contigs_n50_count
        contigs_n50_length
        contigs_t1_count
        contigs_t1_length
        contigs_t1_average_length
        contigs_t1_maximum_length
        contigs_t1_n50_count
        contigs_t1_n50_length
        contigs_t1_n50_not_reached
        contigs_t2_count
        contigs_t2_length
        contigs_t2_average_length
        contigs_t2_maximum_length
        contigs_t2_n50_count
        contigs_t2_n50_length
        contigs_t2_n50_not_reached
        contigs_t3_count
        contigs_t3_length
        contigs_t3_average_length
        contigs_t3_maximum_length
        contigs_t3_n50_count
        contigs_t3_n50_length
        contigs_t3_n50_not_reached

        core_gene_present_percent
        core_gene_group_present_count
        core_gene_survey_result

        coverage_5x
        coverage_4x
        coverage_3x
        coverage_2x
        coverage_1x
        coverage_0x
        reads_attempted
        reads_processed
        reads_processed_length
        reads_processed_success
        reads_processed_average_length
        reads_processed_length_q20
        reads_processed_length_q20_per_read
        reads_processed_length_q20_redundancy
        reads_assembled
        reads_assembled_success
        reads_assembled_success_percent
        reads_assembled_chaff_rate
        reads_assembled_duplicate
        reads_assembled_in_scaffolds
        reads_assembled_unique
        reads_not_assembled
        reads_not_assembled_percent
        scaffolds_1M
        scaffolds_250K_1M
        scaffolds_100K_250K
        scaffolds_10K_100K
        scaffolds_5K_10K
        scaffolds_2K_5K
        scaffolds_0K_2K
        supercontigs_average_length
        supercontigs_count
        supercontigs_length
        supercontigs_minor_count
        supercontigs_minor_length
        supercontigs_minor_average_length
        supercontigs_minor_n50_count
        supercontigs_minor_n50_length
        supercontigs_minor_read_count
        supercontigs_minor_read_percent
        supercontigs_major_read_count
        supercontigs_major_read_percent
        supercontigs_major_average_length
        supercontigs_major_count
        supercontigs_major_length
        supercontigs_major_n50_count
        supercontigs_major_n50_length
        supercontigs_major_percent
        supercontigs_maximum_length
        supercontigs_n50_count
        supercontigs_n50_length
        supercontigs_t1_count
        supercontigs_t1_length
        supercontigs_t1_average_length
        supercontigs_t1_maximum_length
        supercontigs_t1_n50_count
        supercontigs_t1_n50_length
        supercontigs_t1_n50_not_reached
        supercontigs_t2_count
        supercontigs_t2_length
        supercontigs_t2_average_length
        supercontigs_t2_maximum_length
        supercontigs_t2_n50_count
        supercontigs_t2_n50_length
        supercontigs_t2_n50_not_reached
        supercontigs_t3_count
        supercontigs_t3_length
        supercontigs_t3_average_length
        supercontigs_t3_maximum_length
        supercontigs_t3_n50_count
        supercontigs_t3_n50_length
        supercontigs_t3_n50_not_reached
        /)
};

sub add_scaffolds_file {
    my $self = shift;
    return $self->_add_file('add_scaffold', @_);
}

sub add_contigs_file {
    my $self = shift;
    return $self->_add_file('add_contig', @_);
}

sub add_contigs_file_with_contents {
    my $self = shift;
    return $self->_add_file('add_contig_with_contents', @_);
}

sub add_reads_file {
    my $self = shift;
    return $self->_add_file('add_read', @_);
}

sub add_reads_metrics_file {
    my ($self, $file) = @_;

    Carp::confess("No reads metrics file given!") if not $file;
    Carp::confess("Reads metrics file does not exist! $file") if not -s $file;

    my $reads_metrics = Genome::Model::Tools::Sx::Metrics->from_file($file);
    if ( not $reads_metrics ) {
        $self->error_message('Failed to get reads metrics from file! '.$file);
        return;
    }
    $self->reads_processed( $self->reads_processed + $reads_metrics->count );
    $self->reads_processed_length( $self->reads_processed_length + $reads_metrics->bases );

    return 1;
}

sub add_reads_file_with_q20 {
    my $self = shift;
    return $self->_add_file('add_read_with_q20', @_);
}

sub _add_file {
    my ($self, $method, $file_config) = @_;

    Carp::confess('No method for sequence to add file!') if not $method;
    Carp::confess("No file config to add file!") if not $file_config;

    my $reader = Genome::Model::Tools::Sx::Reader->create(config => [ $file_config ]);
    if ( not $reader ) {
        $self->error_message('Failed to open read file: '.$file_config);
        return;
    }

    while ( my $seqs = $reader->read ) {
        for my $seq ( @$seqs ) {
            $self->$method($seq);
        }
    }

    return 1;
}

sub add_scaffold {
    my ($self, $scaffold) = @_;

    Carp::confess('No scaffold to add!') if not $scaffold;

    my $contig_id = 1;
    for my $seq ( split(/n+/i, $scaffold->{seq}) ) {
        $self->add_contig({
                id => $scaffold->{id}.'.'.$contig_id,
                seq => $seq,
            });
        $contig_id++;
    }

    return 1;
}
sub add_contig {
    my ($self, $contig) = @_;

    my $id = $contig->{id};
    $id =~ s/^contig//i;
    my ($supercontig_number, $contig_number) = split(/\./, $id);
    $contig_number = 0 if not defined $contig_number;
    $contig_number = $supercontig_number.'.'.$contig_number;

    my $contig_length = length($contig->{seq});
    $self->supercontigs_count( $self->supercontigs_count + 1 ) if not exists $self->supercontigs->{$supercontig_number};
    $self->supercontigs_length( $self->supercontigs_length + $contig_length );
    $self->supercontigs->{$supercontig_number} += $contig_length;

    $self->contigs_count( $self->contigs_count + 1 );
    $self->assembly_length( $self->assembly_length + $contig_length );
    $self->contigs_length( $self->contigs_length + $contig_length );
    $self->contigs->{$contig_number} = $contig_length;

    if ( $contig_length >= $self->major_contig_threshold ) {
        $self->contigs_major_length( $self->contigs_major_length + $contig_length );
    }
    else {
        $self->contigs_minor_length( $self->contigs_minor_length + $contig_length );
    }

    return 1;
}

sub add_contig_with_contents {
    my ($self, $contig) = @_;

    my $add_contig = $self->add_contig($contig);
    return if not $add_contig;

    my %base_counts = ( a => 0, t => 0, g => 0, C => 0, n => 0, x => 0 );
    foreach my $base ( split('', $contig->{seq}) ) {
        $base_counts{lc $base}++;
    }

    $self->content_at( $self->content_at + $base_counts{a} + $base_counts{t} );
    $self->content_gc( $self->content_gc + $base_counts{g} + $base_counts{c} );
    $self->content_nx( $self->content_nx + $base_counts{n} + $base_counts{x} );

    $self->contigs_length_5k( $self->contigs_length_5k + length($contig->{seq}) ) if length $contig->{seq} >= 5000;

    return 1;
}

sub add_read {
    my ($self, $read) = @_;

    $self->reads_processed( $self->reads_processed + 1 );
    $self->reads_processed_length( $self->reads_processed_length + length($read->{seq}) );

    return 1;
}

sub add_read_with_q20 {
    my ($self, $read) = @_;

    $self->reads_processed( $self->reads_processed + 1 );
    $self->reads_processed_length( $self->reads_processed_length + length($read->{seq}) );
    $self->reads_processed_length_q20( $self->reads_processed_length_q20 + Genome::Model::Tools::Sx::Functions->calculate_qualities_over_minumum($read->{qual}, 20));

    return 1;
}

sub calculate {
    my $self = shift;

    # Supercontigs major/minor lengths
    $self->supercontigs_major_length(0);
    $self->supercontigs_minor_length(0);
    for my $length ( values %{$self->supercontigs} ) {
        if ( $length >= $self->major_contig_threshold ) {
            $self->supercontigs_major_length( $self->supercontigs_major_length + $length );
        }
        else {
            $self->supercontigs_minor_length( $self->supercontigs_minor_length + $length );
        }
    }

    # Reads
    $self->reads_processed_average_length( int($self->reads_processed_length / $self->reads_processed) ) if $self->reads_processed;

    # Tiers
    my $t1 = $self->tier_one;
    my $t2 = $self->tier_two;
    my $t3 = $self->contigs_length - ($t1 + $t2);
    for my $type (qw/ contigs supercontigs /) {
        my $metrics = $self->$type;

        #TOTAL CONTIG VARIABLES
        my $total_contig_number = 0;    my $cumulative_length = 0;
        my $maximum_contig_length = 0;
        my $n50_contig_number = 0;      my $n50_contig_length = 0;
        my $not_reached_n50 = 1;        
        #MAJOR CONTIGS VARIABLES
        my $major_contig_bases = 0;     my $major_contig_number = 0;
        my $major_contig_not_reached_n50 = 1;
        my $major_n50_contig_length = 0;
        my $major_n50_contig_number = 0;
        my $major_contigs_n50_length;
        if ( $type eq 'contigs' ) {
            $major_contigs_n50_length = $self->contigs_major_length * 0.50;
        } else {
            $major_contigs_n50_length = $self->supercontigs_major_length * 0.50;
        }
        #MINOR CONTIG VARIABLES
        my $minor_contig_number = 0;    my $minor_contig_bases = 0;
        my $minor_n50_contig_number = 0;my $minor_n50_contig_length = 0;
        my $minor_n50_not_reached = 1;
        my $minor_contigs_n50_length;
        if ( $type eq 'contigs' ) {
            $minor_contigs_n50_length = $self->contigs_minor_length * 0.50;
        } else {
            $minor_contigs_n50_length = $self->supercontigs_minor_length * 0.50;
        }
        #TIER 1 VARIABLES
        my $total_t1_bases = 0;
        my $t1_n50_contig_number = 0;   my $t1_n50_contig_length = 0;
        my $t1_not_reached_n50 = 1;     my $t1_max_length = 0;
        my $t1_count = 0;
        #TIER 2 VARIABLES
        my $total_t2_bases = 0;
        my $t2_n50_contig_number = 0;   my $t2_n50_contig_length = 0;
        my $t2_not_reached_n50 = 1;     my $t2_max_length = 0;
        my $t2_count = 0;
        #TIER 3 VARIABLES
        my $total_t3_bases = 0;
        my $t3_n50_contig_number = 0;   my $t3_n50_contig_length = 0;
        my $t3_not_reached_n50 = 1;     my $t3_max_length = 0;
        my $t3_count = 0;
        #ASSESS CONTIG / SUPERCONTIG SIZE VARIABLES
        my $larger_than_1M = 0;         my $larger_than_250K = 0;
        my $larger_than_100K = 0;       my $larger_than_10K = 0;
        my $larger_than_5K = 0;         my $larger_than_2K = 0;
        my $larger_than_0K = 0;

        foreach my $c (sort {$metrics->{$b} <=> $metrics->{$a}} keys %{$metrics}) {
            $total_contig_number++;
            $cumulative_length += $metrics->{$c};
            #LONGEST CONTIG
            if ($metrics->{$c} > $maximum_contig_length) {
                $maximum_contig_length = $metrics->{$c};
            }
            #ALL CONTIGS
            if ($not_reached_n50) {
                $n50_contig_number++;
                if ($cumulative_length >= ($self->contigs_length * 0.50)) {
                    $n50_contig_length = $metrics->{$c};
                    $not_reached_n50 = 0;
                }
            }
            #MAJOR CONTIG/SUPERCONTIGS
            if ($metrics->{$c} >= $self->major_contig_threshold) {
                $major_contig_bases += $metrics->{$c};
                $major_contig_number++;
                if ( $major_contig_not_reached_n50 ) {
                    $major_n50_contig_number++;
                    if ( $major_contig_bases >= $major_contigs_n50_length ) { 
                        $major_contig_not_reached_n50 = 0;
                        $major_n50_contig_length = $metrics->{$c};
                    }
                }
            }
            #MINOR CONTIGS/SUPERCONTIGS
            if ($metrics->{$c} < $self->major_contig_threshold) {
                $minor_contig_number++;
                $minor_contig_bases += $metrics->{$c};
                if ( $minor_n50_not_reached ) {
                    $minor_n50_contig_number++;
                    if ( $minor_contig_bases >= $minor_contigs_n50_length ) {
                        $minor_n50_not_reached = 0;
                        $minor_n50_contig_length = $metrics->{$c};
                    }
                }
            }
            #TIER 1
            if ($total_t1_bases < $t1) {
                $total_t1_bases += $metrics->{$c};
                if ($t1_not_reached_n50) {
                    $t1_n50_contig_number++;
                    if ($cumulative_length >= ($t1 * 0.50)) {
                        $t1_n50_contig_length = $metrics->{$c};
                        $t1_not_reached_n50 = 0;
                    }
                }
                $t1_count++;
                if ($t1_max_length == 0) {
                    $t1_max_length = $metrics->{$c}
                }
            }
            #TIER 2
            elsif ($total_t2_bases < $t2) {
                $total_t2_bases += $metrics->{$c};
                if ($t2_not_reached_n50) {
                    $t2_n50_contig_number++;
                    if ($cumulative_length >= ($t2 * 0.50)) {
                        $t2_n50_contig_length = $metrics->{$c};
                        $t2_not_reached_n50 = 0;
                    }
                }
                $t2_count++;
                if ($t2_max_length == 0) {
                    $t2_max_length = $metrics->{$c}
                }
            }
            #TIER 3
            else {
                $total_t3_bases += $metrics->{$c};
                if ($t3_not_reached_n50) {
                    $t3_n50_contig_number++;
                    if ($cumulative_length >= ($t3 * 0.50)) {
                        $t3_n50_contig_length = $metrics->{$c};
                        $t3_not_reached_n50 = 0;
                    }
                }
                $t3_count++;
                if ($t3_max_length == 0) {
                    $t3_max_length = $metrics->{$c}
                }
            }

            #FOR SUPERCONTIGS CONTIGUITY METRICS .. calculated number of contigs > 1M, 250K, 100-250K etc
            if ($metrics->{$c} > 1000000) {
                $larger_than_1M++;
            }
            elsif ($metrics->{$c} > 250000) {
                $larger_than_250K++;
            }
            elsif ($metrics->{$c} > 100000) {
                $larger_than_100K++;
            }
            elsif ($metrics->{$c} > 10000) {
                $larger_than_10K++;
            }
            elsif ($metrics->{$c} > 5000) {
                $larger_than_5K++;
            }
            elsif ($metrics->{$c} > 2000) {
                $larger_than_2K++;
            }
            else {
                $larger_than_0K++;
            }
        }

        $self->{$type.'_average_length'} = int ($cumulative_length / $total_contig_number + 0.50);
        $self->{$type.'_maximum_length'} = $maximum_contig_length;
        $self->{$type.'_n50_length'} = $n50_contig_length;
        $self->{$type.'_n50_count'} = $n50_contig_number;
        $self->{$type.'_major_count'} = $major_contig_number;
        $self->{$type.'_major_length'} = $major_contig_bases;
        $self->{$type.'_major_average_length'} = ($major_contig_number > 0) 
        ? int ($major_contig_bases / $major_contig_number + 0.50) 
        : 0;
        $self->{$type.'_major_percent'} = ( $major_contig_bases > 0 )
        ? sprintf( '%.1f', $major_contig_bases / $cumulative_length * 100 )
        : 0;
        $self->{$type.'_major_n50_length'} = $major_n50_contig_length;
        $self->{$type.'_major_n50_count'} = $major_n50_contig_number;

        $self->{$type.'_t1_length'} = $total_t1_bases;
        $self->{$type.'_t1_count'} = $t1_count;
        $self->{$type.'_t1_average_length'} = ($t1_count > 0) ? int ($total_t1_bases/$t1_count + 0.5) : 0;
        $self->{$type.'_t1_n50_count'} = $t1_n50_contig_number;
        $self->{$type.'_t1_n50_length'} = $t1_n50_contig_length;
        $self->{$type.'_t1_n50_not_reached'} = $t1_not_reached_n50; 
        $self->{$type.'_t1_maximum_length'} = $t1_max_length;

        $self->{$type.'_t2_length'} = $total_t2_bases;
        $self->{$type.'_t2_count'} = $t2_count;
        $self->{$type.'_t2_average_length'} = ($t2_count > 0) ? int ($total_t2_bases/$t2_count + 0.5) : 0;
        $self->{$type.'_t2_n50_count'} = $t2_n50_contig_number;
        $self->{$type.'_t2_n50_length'} = $t2_n50_contig_length;
        $self->{$type.'_t2_n50_not_reached'} = $t2_not_reached_n50; 
        $self->{$type.'_t2_maximum_length'} = $t2_max_length;

        $self->{$type.'_t3_length'} = $total_t3_bases;
        $self->{$type.'_t3_count'} = $t3_count;
        $self->{$type.'_t3_average_length'} = ($t3_count > 0) ? int ($total_t3_bases/$t3_count + 0.5) : 0;
        $self->{$type.'_t3_n50_count'} = $t3_n50_contig_number;
        $self->{$type.'_t3_n50_length'} = $t3_n50_contig_length;
        $self->{$type.'_t3_n50_not_reached'} = $t3_not_reached_n50; 
        $self->{$type.'_t3_maximum_length'} = $t3_max_length;

        $self->{$type.'_minor_count'} = $minor_contig_number;
        $self->{$type.'_minor_length'} = $minor_contig_bases;
        $self->{$type.'_minor_average_length'} = ( $minor_contig_bases > 0 ) ?
        int($minor_contig_bases/$minor_contig_number + 0.5)
        : 0;
        $self->{$type.'_minor_n50_count'} = $minor_n50_contig_number;
        $self->{$type.'_minor_n50_length'} = $minor_n50_contig_length;

        if ($type eq 'supercontigs') {
            $self->{'scaffolds_1M'} = $larger_than_1M;
            $self->{'scaffolds_250K_1M'} = $larger_than_250K;
            $self->{'scaffolds_100K_250K'} = $larger_than_100K;
            $self->{'scaffolds_10K_100K'} = $larger_than_10K;
            $self->{'scaffolds_5K_10K'} = $larger_than_5K;
            $self->{'scaffolds_2K_5K'} = $larger_than_2K;
            $self->{'scaffolds_0K_2K'} = $larger_than_0K;
        }
    }

    if ( $self->reads_processed and $self->reads_attempted ) {
        $self->reads_processed_success( sprintf('%0.3f', $self->reads_processed / $self->reads_attempted) );
    }

    if ( $self->reads_assembled != 0 and $self->reads_processed != 0 ) {
        #if ( $self->reads_assembled ne 'NA' and $self->reads_assembled != 0 and $self->reads_processed ne 'NA' and $self->reads_processed != 0 ) {
        $self->reads_assembled_success(
            sprintf('%0.3f', $self->reads_assembled / $self->reads_processed)
        );
        $self->reads_assembled_success_percent(
             sprintf('%.1f', $self->{reads_assembled_unique} / $self->{reads_processed} * 100)
        );
        $self->reads_not_assembled( $self->reads_assembled - $self->reads_processed );
        $self->reads_not_assembled_percent(
            sprintf('%.1f', ( $self->{reads_processed} - $self->{reads_assembled_unique} ) / $self->{reads_processed} * 100)
        );
    }

    return 1;
}

1;

