package Genome::Model::ReferenceAlignment::Report::ReferenceCoverage;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Report::ReferenceCoverage{
    is => 'Genome::Model::Report',
    has => [
            name => {
                default_value =>'Reference Coverage',
            },
            description => {
                calculate_from => [qw/ model_name build_id /],
                calculate => q|
                    return sprintf(
                        'Reference Coverage for model (Name <%s> Build Id <%s>)',
                        $self->model_name,
                        $self->build_id,
                    );
                |,
            },
        ],
    has_optional => [
        _coverage_hash_ref => { },
    ],
};

sub _add_to_report_xml {
    my $self = shift;

    #Currently no need for 44k nodes
    #unless ($self->_generate_stats) {
    #    $self->error_message('Failed to generate coverage stats report data!');
    #    return;
    #}
    unless ($self->_generate_summary) {
        $self->error_message('Failed to generate coverage summary report data!');
        return;
    }
    unless ($self->_generate_breakdown) {
        $self->error_message('Failed to generate breakdown report data!');
        return;
    }
    unless ($self->_generate_progression) {
        $self->error_message('Failed to generate coverage progression report data!');
        return;
    }

    unless ($self->_generate_relative_coverage) {
        $self->error_message('Failed to generate relative coverage report data!');
        return;
    }
    unless ($self->_generate_breadth) {
        $self->error_message('Failed to generate breadth of coverage report data!');
        return;
    }
    unless ($self->_generate_size_histos) {
        $self->error_message('Failed to generate size histograms of coverage report data!');
        return;
    }
    return 1;
}

sub _generate_summary {
    my $self = shift;

    my $ref = $self->_get_coverage_summary_hash_ref;
    my $total_references;
    my $touched_references;
    foreach my $ref_name (keys %{$ref}) {
        my $cov    = $ref->{$ref_name}->{cov};
        if ($cov > 0) {
            $touched_references++;
        }
        $total_references++;
    }
    my $pc_touched = sprintf("%.02f",( ($touched_references / $total_references) * 100 ));
    my @headers = qw(category-name category-items total percent-category);
    my @rows;

    my @touched_row = ('references-touched', $touched_references, $total_references, $pc_touched);
    push @rows, \@touched_row;

    my $total_reads = 0;
    my $aligned_reads = 0;
    my @instrument_data = $self->build->instrument_data;
    for my $instrument_data (@instrument_data) {
        my @alignments = $self->build->alignment_results_for_instrument_data($instrument_data);
        for my $alignment (@alignments) {
            $total_reads += $alignment->total_read_count;
            $aligned_reads += $alignment->total_aligned_read_count;
            #May need to handle PE differently
        }
    }
    my $pc_aligned;
    if ($total_reads > 0) {
        $pc_aligned = sprintf("%.02f",( ($aligned_reads / $total_reads) * 100 ));
    } else {
        $pc_aligned = '0.00';
    }
    my @align_row = ('aligned-reads', $aligned_reads, $total_reads, $pc_aligned);
    push @rows, \@align_row;
    unless ($self->_add_dataset(
        name => 'summary',
        title => 'cDNA Alignment Summary',
        'display-type' => 'table',
        row_name => 'category',
        headers => \@headers,
        rows => \@rows,
    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}

sub _generate_stats {
    my $self = shift;

    my $stats_tsv = $self->build->coverage_stats_file;
    my @headers = (
        'reference-name','percent-covered','reference-bases','covered-bases','missing-bases',
        'average-coverage-depth','average-coverage-depth-sd','median-coverage-depth','gaps','average-gap-length',
        'average-gap-length','median-gap-length','minimum-depth-filter','discarded-bases','percent-discarded-bases'
    );
    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $stats_tsv,
        separator => "\t",
        headers => \@headers,
    );
    unless ($parser) {
        $self->error_message('Failed to create tab delimited parser for file '. $stats_tsv);
        die($self->error_message);
    }
    my @rows;
    while (my $stats = $parser->next) {
        my @row;
        for my $header (@headers) {
            push @row, $stats->{$header};
        }
        push @rows, \@row;
    }
    unless ($self->_add_dataset(
                        name => 'coverage-stats',
                        title => 'Coverage Statisticts Per Reference',
                        'display-type' => 'table',
                        row_name => 'reference',
                        headers => \@headers,
                        rows => \@rows,
                    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}

sub _generate_breakdown {
    my $self = shift;
    my $breakdown_tsv = $self->build->breakdown_file;

    my @headers = ('category-name','reads','percent-total-reads','percent-category','poly-at',
                   'percent-poly-at','poly-a','percent-poly-a','poly-t','percent-poly-t');

    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $breakdown_tsv,
        separator => "\t",
        headers => \@headers,
    );
    unless ($parser) {
        $self->error_message('Failed to create tab delimited parser for file '. $breakdown_tsv);
        die($self->error_message);
    }
    my @rows;
    while (my $breakdown = $parser->next) {
        my @row;
        for my $header (@headers) {
            push @row, $breakdown->{$header};
        }
        push @rows, \@row;
    }
    unless ($self->_add_dataset(
                        name => 'breakdown',
                        title => 'Alignment Breakdown',
                        'display-type' => 'table',
                        row_name => 'category',
                        headers => \@headers,
                        rows => \@rows,
                    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}

sub _generate_progression {
    my $self = shift;

    my $progression_tsv = $self->build->coverage_progression_file;

    my @headers = ('lane-name','bases-covered','reference-bases','percent-coverage','percent-gain','lane-id','lane-error-average','lane-error-standard-deviation');

    my $parser = Genome::Utility::IO::SeparatedValueReader->create(
        input => $progression_tsv,
        separator => "\t",
        headers => \@headers,
    );
    unless ($parser) {
        $self->error_message('Failed to create tab delimited parser for file '. $progression_tsv);
        die($self->error_message);
    }
    my @rows;
    while (my $stats = $parser->next) {
        my @row;
        for my $header (@headers) {
            push @row, $stats->{$header};
        }
        push @rows, \@row;
    }
    unless ($self->_add_dataset(
                        name => 'progression',
                        title => 'Coverage Progression By Lane',
                        'display-type' => 'graph-line',
                        'x-axis' => 'lane-name',
                        'y1-axis' => 'percent-coverage',
                        'y1-axis-title' => 'Total Coverage(%)',
                        'y2-axis' => 'percent-gain',
                        'y2-axis-title' => 'Coverage Gain(%)',
                        row_name => 'lane',
                        headers => \@headers,
                        rows => \@rows,
                    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}

sub _generate_relative_coverage {
    my $self = shift;

    my $relative_coverage_node = $self->_xml->createElement('relative-coverage')
        or return;
    $self->_datasets_node->addChild($relative_coverage_node)
        or return;
    my %params = (
        title => 'Depth By Relative Position',
        'display-type' => 'graph-line',
        'x-axis-title' => 'Relative Position(5\'->3\')',
        'x-axis' => 'relative-position',
        'y-axis-title' => 'Depth',
        'y-axis' => 'depth',
    );
    for my $attr (keys %params) {
        $relative_coverage_node->addChild( $self->_xml->createAttribute($attr, $params{$attr}) )
            or return;
    }

    my @headers = ('relative-position','depth');
    my @size_files = $self->build->relative_coverage_files;
    for my $size_file ( @size_files ) {
        unless ($size_file =~ /bias_(\w+)$/) {
            die('Failed to parse bias file path '. $size_file);
        }
        my $size = lc($1);
        my $size_bin_node = $self->_xml->createElement('relative-coverage-bin')
            or return;
        $size_bin_node->addChild( $self->_xml->createAttribute('size', $size) )
            or return;
        $relative_coverage_node->addChild($size_bin_node)
            or return;
        my $parser = Genome::Utility::IO::SeparatedValueReader->create(
            input => $size_file,
            separator => "\t",
            headers => \@headers,
        );
        unless ($parser) {
            $self->error_message('Failed to create tab delimited parser for file '. $size_file);
            die($self->error_message);
        }
        while (my $rel_cov = $parser->next) {
            my $reading_node = $size_bin_node->addChild( $self->_xml->createElement('reading') );
            for my $header (@headers) {
                my $element = $reading_node->addChild( $self->_xml->createElement($header) )
                    or return;
                $element->appendTextNode($rel_cov->{$header});
            }
        }
    }
    return 1;
}

sub _generate_breadth {
    my $self = shift;

    my $ref = $self->_get_coverage_summary_hash_ref;

    # TODO: make the ranges a set on input params
    my %bins;
    foreach my $ref_name (keys %{$ref}) {
        my $size = $ref->{$ref_name}->{size};
        my $cov    = $ref->{$ref_name}->{cov};
        if (($size >= 100) && ($size <= 2_999)) {
            if ($cov >= 90) {
                $bins{SMALL_COVERED}++;
            }
            $bins{SMALL_TOTAL}++;
        }
        elsif (($size >= 3_000) && ($size <= 6_999)) {
            if ($cov >= 50) {
                $bins{MEDIUM_COVERED}++;
            }
            $bins{MEDIUM_TOTAL}++;
        }
        elsif (($size >= 7_000)) {
            if ($cov >= 30) {
                $bins{LARGE_COVERED}++;
            }
            $bins{LARGE_TOTAL}++;
        }
    }
    my %desc = (
                SMALL => {
                    minimum_basepair => '100',
                    maximum_basepair => '<=2999',
                    minimum_depth => '1X',
                    minimum_percent_coverage =>'90',
                },
                MEDIUM => {
                    minimum_basepair => '>=3000',
                    maximum_basepair => '<=6999',
                    minimum_depth => '1X',
                    minimum_percent_coverage =>'50',
                },
                LARGE => {
                    minimum_basepair => '>=7000',
                    maximum_basepair => 'inf',
                    minimum_depth => '1X',
                    minimum_percent_coverage =>'30',
                },
            );
    my @bin_headers= qw(reference-size minimum-basepair maximum-basepair minimum-depth minimum-percent-coverage covered-references total-references percent-covered);
    # Hard coded to retain order
    my @rows;
    for my $key ('SMALL', 'MEDIUM', 'LARGE') {
        my $covered = $bins{$key .'_COVERED'};
        my $total = $bins{$key .'_TOTAL'};
        my $pc = sprintf("%.02f",(($covered / $total) * 100));
        my @row = ($key,$desc{$key}->{minimum_basepair},$desc{$key}->{maximum_basepair},$desc{$key}->{minimum_depth},$desc{$key}->{minimum_percent_coverage},$covered,$total,$pc);
        push @rows, \@row;
    }
    unless ($self->_add_dataset(
        name => 'coverage-bins',
        title => 'Breadth and Depth of Coverage Summary',
        row_name => 'coverage-bin',
        'display-type' => 'table',
        headers => \@bin_headers,
        rows => \@rows,
    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}


sub _generate_size_histos {
    my $self = shift;

    my $ref = $self->_get_coverage_summary_hash_ref;

    # Walk through reference BIN sizes and act on each bin independently.
    my $bin_sizes = {
                     1  => [1,         500],
                     2  => [501,       1_000],
                     3  => [1_001,     2_000],
                     4  => [2_001,     3_000],
                     5  => [3_001,     4_000],
                     6  => [4_001,     5_000],
                     7  => [5_001,     6_000],
                     8  => [6_001,     7_000],
                     9  => [7_001,     8_000],
                     10 => [8_001,     9_000],
                     11 => [9_001,     10_000],
                     12 => [10_001,    11_000],
                     13 => [11_001,    12_000],
                     14 => [12_001,    13_000],
                     15 => [13_001,    14_000],
                     16 => [14_001,    15_000],
                     17 => [15_000, 1_000_000],  # catch all
                 };
    my $range_sizes = {
                       1  => [0,     0],
                       2  => [0.01,  10],
                       3  => [10.01, 20],
                       4  => [20.01, 30],
                       5  => [30.01, 40],
                       6  => [40.01, 50],
                       7  => [50.01, 60],
                       8  => [60.01, 70],
                       9  => [70.01, 80],
                       10 => [80.01, 90],
                       11 => [90.01, 99.99],
                       12 => [100,  100],
                   };
    my @headers = qw(minimum-basepair maximum-basepair minimum-coverage-percent maximum-coverage-percent references-covered percent-references-covered);
    my @rows;
    foreach my $bin (sort {$a <=> $b} keys %{$bin_sizes}) {
        # Gather information based on BIN size.
        my $bininfo = {};

        # Check for reference inclusion in BIN.
      BININCLUSION:
        foreach my $ref_name (keys %{$ref}) {
            if ( $ref->{$ref_name}->{size} >= $bin_sizes->{$bin}[0] &&
                 $ref->{$ref_name}->{size} <= $bin_sizes->{$bin}[1] ) {
                # Inclusion.
                $bininfo->{total_num_genes}++;
              RANGEINCLUSION:
                foreach my $range (sort {$a <=> $b} keys %{$range_sizes}) {
                    if ( $ref->{$ref_name}->{cov} >= $range_sizes->{$range}[0] &&
                         $ref->{$ref_name}->{cov} <= $range_sizes->{$range}[1] ) {
                        $bininfo->{perc_cov}->{$range}++;
                        $bininfo->{depth_cov}->{$range}++;
                    } else {
                        next RANGEINCLUSION;
                    }
                }
            } else {
                # Exclusion.
                next BININCLUSION;
            }
        }

        foreach my $range_val (sort {$a <=> $b} keys %{$range_sizes}) {
            my @row;
            if ($bininfo->{perc_cov}->{$range_val}) {
                my $percent = ($bininfo->{perc_cov}->{$range_val} / $bininfo->{total_num_genes}) * 100;
                $percent    = sprintf( "%.2f", $percent );
                @row = ($bin_sizes->{$bin}[0],$bin_sizes->{$bin}[1],$range_sizes->{$range_val}[0],$range_sizes->{$range_val}[1],$bininfo->{perc_cov}->{$range_val},$percent);
            } else {
                @row = ($bin_sizes->{$bin}[0],$bin_sizes->{$bin}[1],$range_sizes->{$range_val}[0],$range_sizes->{$range_val}[1],'0','0');
            }
            push @rows, \@row;
        }
    }
    unless ($self->_add_dataset(
        name => 'size-bins',
        row_name => 'size-bin',
        title => 'Breadth and Depth of Coverage',
        headers => \@headers,
        rows => \@rows,
    )) {
        $self->error_message('Failed to add dataset.');
        return;
    }
    return 1;
}


sub _get_coverage_summary_hash_ref {
    my $self = shift;
    unless ( $self->_coverage_hash_ref ) {
        my $stats_tsv = $self->build->coverage_stats_file;

        my @stats_headers = (
            'reference_name','pc_covered','reference_bases','covered_bases','missing_bases',
            'avg_coverage_depth','sd_avg_coverage_depth','med_coverage_depth','gaps','avg_gap_length',
            'sd_avg_gap_length','med_gap_length','min_depth','discarded_bases','pc_discarded_bases'
        );
        my $parser = Genome::Utility::IO::SeparatedValueReader->create(
            input => $stats_tsv,
            separator => "\t",
            headers => \@stats_headers,
        );
        unless ($parser) {
            $self->error_message('Failed to create tab delimited parser for file '. $stats_tsv);
            die($self->error_message);
        }
        my %hash;
        while (my $stats = $parser->next) {
            my $ref             = $stats->{'reference_name'};
            $hash{$ref}->{cov}   = $stats->{'pc_covered'};
            $hash{$ref}->{size}  = $stats->{'reference_bases'};
            $hash{$ref}->{depth} = $stats->{'avg_coverage_depth'};
        }
        $self->_coverage_hash_ref(\%hash);
    }
    return $self->_coverage_hash_ref;
}

1;
