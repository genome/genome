package Genome::Model::Tools::RefCov::Stats;

use strict;
use warnings;

use Statistics::Descriptive;
use PDL;
use PDL::NiceSlice;

my @HEADERS = qw/
                    name
                    percent_ref_bases_covered
                    total_ref_bases
                    total_covered_bases
                    missing_bases
                    ave_cov_depth
                    sdev_ave_cov_depth
                    med_cov_depth
                    gap_number
                    ave_gap_length
                    sdev_ave_gap_length
                    med_gap_length
                    min_depth_filter
                    min_depth_discarded_bases
                    percent_min_depth_discarded
/;

my %HEADER_DESCRIPTIONS = (
    name                        => 'Name of Coverage Region',
    percent_ref_bases_covered   => 'Percent of Reference Bases Covered',
    total_ref_bases             => 'Total Number of Reference Bases',
    total_covered_bases         => 'Total Number of Covered Bases',
    missing_bases               => 'Number of Missing Bases',
    ave_cov_depth               => 'Average Coverage Depth',
    sdev_ave_cov_depth          => 'Standard Deviation Average Coverage Depth',
    med_cov_depth               => 'Median Coverage Depth',
    gap_number                  => 'Number of Gaps',
    ave_gap_length              => 'Average Gap Length',
    sdev_ave_gap_length         => 'Standard Deviation Average Gap Length',
    med_gap_length              => 'Median Gap Length',
    min_depth_filter            => 'Min. Depth Filter',
    min_depth_discarded_bases   => 'Discarded Bases (Min. Depth Filter)',
    percent_min_depth_discarded => 'Percent Discarded Bases (Min. Depth Filter)',

);

class Genome::Model::Tools::RefCov::Stats {
    has => [
        #coverage => {
        #    is => 'ArrayRef',
        #    doc => 'An array of integers representing the depth of coverage at each base position',
        #},
        name => {
            is => 'String',
            doc => 'A unique name for the coverage region',
        },
        min_depth => {
            is => 'Integer',
            doc => 'The minimum depth to consider a base position covered',
            is_optional => 1,
            default_value => 1,
        },
        stats_mode => {
            default_value => 'normal',
            valid_values => ['brief','normal'],
            is_optional => 1,
        },
    ],
    has_optional => {
        percent_ref_bases_covered   => {},
        #Redundant with ref_length
        total_ref_bases             => {},
        total_covered_bases         => {},
        missing_bases               => {},
        ave_cov_depth               => {},
        sdev_ave_cov_depth          => {},
        med_cov_depth               => {},
        gap_number                  => {},
        ave_gap_length              => {},
        sdev_ave_gap_length         => {},
        med_gap_length              => {},
        min_depth_filter            => {},
        min_depth_discarded_bases   => {},
        percent_min_depth_discarded => {},
        #OPTIONAL
        nonzero_minimum_coverage_depth => { },
        minimum_coverage_depth => { },
        maximum_coverage_depth => { },
    },
};

sub create {
    my $class = shift;
    my %params = @_;

    my $coverage = delete($params{coverage});
    my $pdl = delete($params{pdl});

    my $self = $class->SUPER::create(%params);
    $PDL::BIGPDL=1;

    if ($coverage) {
        $pdl = pdl @{$coverage};
    }
    if (defined $pdl) {
        $self->{_coverage_pdl} = $pdl;
        $self->_main_calculation_code;
    }
    return $self;
}

sub coverage {
    my $self = shift;
    my $coverage_pdl = $self->{_coverage_pdl};
    my @coverage = $coverage_pdl->list;
    return \@coverage;
}

sub ref_length {
    my $self = shift;
    my $coverage_pdl = $self->{_coverage_pdl};
    return $coverage_pdl->nelem;
}

# This method is poorly named, it's really a "reload" method with new params.
sub calculate_coverage_stats {
    my $self = shift;
    my %params = @_;
    my $coverage = delete($params{coverage});
    my $pdl = delete($params{pdl});
    unless ($coverage || defined($pdl)) {
        die('A coverage array ref or PDL is required!');
    }
    my $min_depth = delete($params{min_depth});
    my $name = delete($params{name});
    if ($coverage) {
        $pdl = pdl @{$coverage};
    }
    $self->{_coverage_pdl} = $pdl;
    $self->min_depth($min_depth);
    $self->name($name);
    $self->_main_calculation_code;
    return $self;
}

sub add_coverage {
    my $self = shift;
    my %params = @_;
    my $coverage = delete($params{coverage});
    my $pdl = delete($params{pdl});
    unless ($coverage || defined($pdl)) {
        die('A coverage array ref or PDL is required!');
    }
    my $min_depth = delete($params{min_depth});
    my $name = delete($params{name});
    if ($coverage) {
        $pdl = pdl @{$coverage};
    }
    my $coverage_pdl = $self->{_coverage_pdl};
    unless ($pdl->nelem eq $coverage_pdl->nelem) {
        die('The new coverage length is '. $pdl->nelem .' but the existing length is '. $coverage_pdl->nelem);
    }
    $coverage_pdl += $pdl;
    if (defined($min_depth)) {
        $self->min_depth($min_depth);
    }
    if (defined($name)) {
        $self->name($name);
    }
    $self->_main_calculation_code;
    return $self;
}

# This seems to work on up to 30-40Mb PDLs

sub append_coverage {
    my $self = shift;

    my %params = @_;
    my $coverage = delete($params{coverage});
    my $pdl = delete($params{pdl});
    unless ($coverage || defined($pdl)) {
        die('A coverage array ref or PDL is required!');
    }
    my $coverage_pdl = $self->{_coverage_pdl};
    unless (defined($coverage_pdl)) {
        die('In order to append_coverage there must be an existing coverage PDL!');
    }
    if ($coverage) {
        $pdl = pdl @{$coverage};
    }
    $coverage_pdl = $coverage_pdl->append($pdl);
    $self->{_coverage_pdl} = $coverage_pdl;
    $self->_main_calculation_code;
    return $self;
}

sub _main_calculation_code {
    my $self = shift;
    my $calculation_method = '_'. $self->stats_mode .'_calculation_code';
    return $self->$calculation_method;
}

# This generates basic statistics about coverage.
# This does not track gaps.
# But it works on much larger PDL sizes
sub _brief_calculation_code {
    my $self = shift;
    my $coverage_pdl = $self->{_coverage_pdl};

    my ($non_zero_bp,$zero_bp) = $self->zero_filtered;

    my ($total_covered, $lt_min_depth_bp, $min_depth_pdl) = $self->min_depth_filtered;
    my $discarded_bases = ($lt_min_depth_bp - $zero_bp);
    $self->total_covered_bases( $total_covered );
    
    my ($mean_pos_depth,$prms_pos_depth,$med_pos_depth,$coverage_min,$coverage_max,$adev_pos_depth,$coverage_rms) = $min_depth_pdl->stats;
    
    # --------------------------------------------------
    # F O R M A T
    # --------------------------------------------------
    # [0]   Percent of Reference Bases Covered
    # [1]   Total Number of Reference Bases
    # [2]   Total Number of Covered Bases
    # [3]   Number of Missing Bases
    # [4]   Average Coverage Depth
    # [5]   Standard Deviation Average Coverage Depth
    # [6]   Median Coverage Depth
    # (NOT SUPPORTED) [7]   Number of Gaps
    # (NOT SUPPORTED) [8]   Average Gap Length
    # (NOT SUPPORTED) [9]   Standard Deviation Average Gap Length
    # (NOT SUPPORTED) [10]  Median Gap Length
    # [11]  Min. Depth Filter
    # [12]  Discarded Bases (Min. Depth Filter)
    # [13]  Percent Discarded Bases (Min. Depth Filter)
    # (DEPRECATED)  [14]  Max. Unique Filter
    # (DEPRECATED)  [15]  Total Number of Reads Layered
    # (DEPRECATED)  [16]  Total Number of Unique Start Site Reads Layered
    # (DEPRECATED)  [17]  Percent Redundancy of Read Layers
    # (DEPRECATED)  [18]  Zenith
    # (DEPRECATED)  [19]  Nadir
    # --------------------------------------------------

    # [0] Percent of Reference Bases Covered
    $self->percent_ref_bases_covered( _round( ($self->total_covered_bases / $self->ref_length()) * 100 ) );

    # [1] Total Number of Reference Bases
    $self->total_ref_bases( $self->ref_length() );

    # [2] Total Number of Covered Bases
    # (... set above.)

    # [3] Number of Missing Bases
    $self->missing_bases( $self->ref_length() - $self->total_covered_bases() );

    # [4] Average Coverage Depth
    $self->ave_cov_depth( _round( $mean_pos_depth ) );

    # [5] Standard Deviation Average Coverage Depth
    $self->sdev_ave_cov_depth( _round( $prms_pos_depth ) );

    # [6] Median Coverage Depth
    $self->med_cov_depth( _round( $med_pos_depth ) );

    # [7] Number of Gaps
    $self->gap_number( 'n/a' );

    # [8] Average Gap Length
    $self->ave_gap_length( 'n/a' );

    # [9] Standard Deviation Average Gap Length
    $self->sdev_ave_gap_length( 'n/a' );

    # [10] Median Gap Length
    $self->med_gap_length( 'n/a' );

    # [11] Min. Depth Filter
    $self->min_depth_filter( $self->min_depth );

    # [12] Discarded Bases (Min. Depth Filter)
    $self->min_depth_discarded_bases( $discarded_bases );

    # [13] Percent Discarded Bases (Min. Depth Filter)
    $self->percent_min_depth_discarded( _round( ($self->min_depth_discarded_bases() / $self->total_ref_bases()) * 100 ) );

    # OPTIONAL COLUMNS
    $self->minimum_coverage_depth($coverage_min->sclr);
    $self->maximum_coverage_depth($coverage_max->sclr);

    return $self;
}

sub zero_filtered {
    my $self = shift;
    my $coverage_pdl = $self->{_coverage_pdl};

    my ($non_zero_idx,$zero_idx) = which_both($coverage_pdl);
    my $non_zero_bp = $non_zero_idx->nelem;
    my $zero_bp = $zero_idx->nelem;
    return ($non_zero_bp,$zero_bp);
}

sub min_depth_filtered {
    my $self = shift;
    my $coverage_pdl = $self->{_coverage_pdl};

    my ($gt_min_depth_idx,$lt_min_depth_idx) = which_both($coverage_pdl >= $self->min_depth);
    my $min_depth_pdl = $coverage_pdl->copy;
    my $lt_pdl = $min_depth_pdl->index($lt_min_depth_idx);
    $lt_pdl .= 0;
    # TODO: When a pdl is switched out there should be a method that clears the boject values.
    # These objects are reused because of the overhead in creating millions of objects
    $self->{_min_depth_filtered_pdl} = $min_depth_pdl;
    return ($gt_min_depth_idx->nelem,$lt_min_depth_idx->nelem,$min_depth_pdl);
}

sub min_depth_filtered_coverage {
    my $self = shift;
    my $min_depth_pdl = $self->{_min_depth_filtered_pdl};
    unless (defined($min_depth_pdl)) {
        (undef,undef,$min_depth_pdl) = $self->min_depth_filtered;
    }
    my @filtered_coverage = $min_depth_pdl->list;
    return \@filtered_coverage;
}

sub _normal_calculation_code {
    my $self = shift;

    # ** NOTE **
    # Thu Jul  9 21:35:35 CDT 2009
    # The min_depth filter permanently changes the incoming coverage-depth
    # string; all downstream interactions will be based on the revised
    # coverage-depth string values, where values not matching the minimum
    # criteria are set to "0".

    # ** NOTE **
    # Thu Jul  9 21:35:46 CDT 2009
    # The functions related to "redundancy" of start-sites has been deprecated
    # in this version. The values needed for these calcualtions was previously
    # provided by the RefCov layering engine. We may possibly be able to
    # re-introduce this functionality by delving into the pileup information
    # provided by the Bio::DB::Bam package--however, for now, we will not be
    # attempting redundancy calculations.

    # MAIN CALCULATIONS FOLLOW THIS LINE
    # _____________________________________________________________________________

    my $min_depth = $self->min_depth;

    my ($non_zero_bp, $zero_bp) = $self->zero_filtered;

    my ($total_covered, $lt_min_depth_bp, $min_depth_pdl) = $self->min_depth_filtered;
    my $discarded_bases = ($lt_min_depth_bp - $zero_bp);
    $self->total_covered_bases( $total_covered );


    # COVERAGE
    my ($mean_pos_depth,$prms_pos_depth,$med_pos_depth,$coverage_min,$coverage_max,$adev_pos_depth,$coverage_rms) = $min_depth_pdl->stats;

    # GAPS
    my ($quantity,$uniq_vals) = rle($min_depth_pdl);
    my $gaps_pdl_zero_padded = $quantity->where($uniq_vals < $min_depth);
    my $gaps_pdl = $gaps_pdl_zero_padded->where($gaps_pdl_zero_padded > 0);
    my $gaps = $gaps_pdl->nelem;
    my $mean_gap_size = 0;
    my $prms_gap_size = 0;
    my $med_gap_size = 0;
    my $gap_min = 0;
    my $gap_max = 0;
    my $adev_gap_size = 0;
    my $gap_rms = 0;
    if ($gaps) {
        ($mean_gap_size,$prms_gap_size,$med_gap_size,$gap_min,$gap_max,$adev_gap_size,$gap_rms) = $gaps_pdl->stats;
    }
    # --------------------------------------------------
    # F O R M A T
    # --------------------------------------------------
    # [0]   Percent of Reference Bases Covered
    # [1]   Total Number of Reference Bases
    # [2]   Total Number of Covered Bases
    # [3]   Number of Missing Bases
    # [4]   Average Coverage Depth
    # [5]   Standard Deviation Average Coverage Depth
    # [6]   Median Coverage Depth
    # [7]   Number of Gaps
    # [8]   Average Gap Length
    # [9]   Standard Deviation Average Gap Length
    # [10]  Median Gap Length
    # [11]  Min. Depth Filter
    # [12]  Discarded Bases (Min. Depth Filter)
    # [13]  Percent Discarded Bases (Min. Depth Filter)
    # (DEPRECATED)  [14]  Max. Unique Filter
    # (DEPRECATED)  [15]  Total Number of Reads Layered
    # (DEPRECATED)  [16]  Total Number of Unique Start Site Reads Layered
    # (DEPRECATED)  [17]  Percent Redundancy of Read Layers
    # (DEPRECATED)  [18]  Zenith
    # (DEPRECATED)  [19]  Nadir
    # --------------------------------------------------

    # [0] Percent of Reference Bases Covered
    $self->percent_ref_bases_covered( _round( ($self->total_covered_bases / $self->ref_length()) * 100 ) );

    # [1] Total Number of Reference Bases
    $self->total_ref_bases( $self->ref_length() );

    # [2] Total Number of Covered Bases
    # (... set above.)

    # [3] Number of Missing Bases
    $self->missing_bases( $self->ref_length() - $self->total_covered_bases() );

    # [4] Average Coverage Depth
    $self->ave_cov_depth( _round( $mean_pos_depth ) );

    # [5] Standard Deviation Average Coverage Depth
    $self->sdev_ave_cov_depth( _round( $prms_pos_depth ) );

    # [6] Median Coverage Depth
    $self->med_cov_depth( _round( $med_pos_depth ) );

    # [7] Number of Gaps
    if ($gaps) {
        $self->gap_number( $gaps );
    } else {
        $self->gap_number( '0' );
    }

    # [8] Average Gap Length
    $self->ave_gap_length( _round( $mean_gap_size ) );

    # [9] Standard Deviation Average Gap Length
    $self->sdev_ave_gap_length( _round( $prms_gap_size ) );

    # [10] Median Gap Length
    $self->med_gap_length( _round( $med_gap_size ) );

    # [11] Min. Depth Filter
    $self->min_depth_filter( $min_depth );

    # [12] Discarded Bases (Min. Depth Filter)
    $self->min_depth_discarded_bases( $discarded_bases );

    # [13] Percent Discarded Bases (Min. Depth Filter)
    $self->percent_min_depth_discarded( _round( ($self->min_depth_discarded_bases() / $self->total_ref_bases()) * 100 ) );

    # OPTIONAL COLUMNS
    $self->minimum_coverage_depth($coverage_min->sclr);
    $self->maximum_coverage_depth($coverage_max->sclr);

    return $self;
}

sub _round {
    my $value = shift;
    return sprintf( "%.2f", $value );
}

sub stats_array_ref {
    my $self = shift;
    my @stats;
    for my $header ($self->headers) {
        if (defined($self->$header)) {
            push @stats, $self->$header;
        }
    }
    return \@stats;
}

sub stats_hash_ref {
    my $self = shift;
    my %stats;
    for my $header ($self->headers) {
        if (defined($self->$header)) {
            $stats{$header} = $self->$header;
        }
    }
    return \%stats;
}

sub stats_index {
    my $self = shift;

    my @headers = $self->headers;
    my @descriptions = $self->header_descriptions;
    unless (scalar(@headers) == scalar(@descriptions)) {
        die('Failed to find the correct number of headers and descriptions.');
    }
    my %stats;
    for (my $i = 0; $i < scalar(@headers); $i++) {
        my $header = $headers[$i];
        my $description = $descriptions[$i];
        if (defined($self->$header)) {
            $stats{$i} = { $description => $self->$header};
        }
    }
    return \%stats;
}

sub headers {
    return @HEADERS;
}

sub header_descriptions{
    my @descriptions;
    for my $header (@HEADERS) {
        push @descriptions, $HEADER_DESCRIPTIONS{$header};
    }
    return @descriptions;
}

sub save_stats {
    my ($self, $file) = @_;

    # Require a file path.
    if (!$file) { croak (__PACKAGE__ . ' method save_stats requires a "file" argument.'); }

    # Order stats by field keys and print to STDOUT, attach units
    open (OUT, ">$file") or die 'could not open save file for stats';
    foreach my $field_key (sort {$a <=> $b} keys %{ $self->stats_index() }) {
        foreach my $field_name (keys %{ $self->stats_index()->{$field_key} }) {
            if (defined($self->stats_index()->{$field_key}->{$field_name})) {
                my $line = sprintf( "%-45s %-15s\n", $field_name . ':', $self->stats_index()->{$field_key}->{$field_name});
                print OUT $line;
            }
        }
    }
    close (OUT);

    return $self;
}

sub print_stats {
    my $self = shift;

    # Order stats by field keys and print to STDOUT, attach units
    foreach my $field_key (sort {$a <=> $b} keys %{ $self->stats_index() }) {
        foreach my $field_name (keys %{ $self->stats_index()->{$field_key} }) {
            if (defined($self->stats_index()->{$field_key}->{$field_name})) {
                printf( "%-45s %-15s\n", $field_name . ':', $self->stats_index()->{$field_key}->{$field_name});
            }
        }
    }

    return $self;
}

sub as_string {
    my $self = shift;
    my $string = join ("\t", @{ $self->stats_array_ref() });
    return $string;
}

# DEPRECATED: see generate_clusters method below
sub clusters {
    my $self = shift;

    my $chr_pdl = $self->{_coverage_pdl};

    my ($quantity,$value) = rle($chr_pdl >= $self->min_depth);
    my ($padded_quantity,$padded_value) = $quantity->where($value,$quantity!=0);
    my $padded_quantity_cumsum = $padded_quantity->cumusumover;

    my ($chr_gt_min_depth_idx,$chr_lt_min_depth_idx) = which_both($padded_value);

    my $gt_first = $chr_gt_min_depth_idx(0)->sclr;
    my $lt_first = $chr_lt_min_depth_idx(0)->sclr;

    my $gt_last = $chr_gt_min_depth_idx(-1)->sclr;
    my $lt_last = $chr_lt_min_depth_idx(-1)->sclr;

    my $start = $padded_quantity_cumsum->index($chr_lt_min_depth_idx);
    my $end = $padded_quantity_cumsum->index($chr_gt_min_depth_idx);

    if ($gt_first < $lt_first) {
        # Coverage first because the first index is zero
        # Add the initial start position(ie. zero)
        $start = append(zeroes(1),$start);
    }
    if ($gt_last < $lt_last) {
        # Gap on the end, remove the last value
        $start = $start(0:-2);
    }
    my @clusters;
    for (my $i = 0; $i < $start->getdim(0); $i++) {
        # Start is zero-based and End is one-based, just like BED format
        push @clusters, [$start($i)->sclr,$end($i)->sclr];
    }
    return \@clusters;
}

sub generate_clusters {
    my $self = shift;
    my $offset = shift || 0;

    my $chr_pdl = $self->{_coverage_pdl};
    unless (defined($chr_pdl)) { return; }
    my ($quantity,$value) = rle($chr_pdl >= $self->min_depth);
    my ($padded_quantity,$padded_value) = $quantity->where($value,$quantity!=0);
    my $padded_quantity_cumsum = $padded_quantity->cumusumover;

    my ($chr_gt_min_depth_idx,$chr_lt_min_depth_idx) = which_both($padded_value);
    unless ($chr_gt_min_depth_idx->nelem) {
        return;
    }
    my $gt_first = $chr_gt_min_depth_idx(0)->sclr;
    my $lt_first = $chr_lt_min_depth_idx(0)->sclr;

    my $gt_last = $chr_gt_min_depth_idx(-1)->sclr;
    my $lt_last = $chr_lt_min_depth_idx(-1)->sclr;

    my $start = $padded_quantity_cumsum->index($chr_lt_min_depth_idx);
    my $end = $padded_quantity_cumsum->index($chr_gt_min_depth_idx);

    if ($gt_first < $lt_first) {
        # Coverage first because the first index is zero
        # Add the initial start position(ie. zero)
        $start = append(zeroes(1),$start);
    }
    if ($gt_last < $lt_last) {
        # Gap on the end, remove the last value
        $start = $start(0:-2);
    }
    my @clusters;
    for (my $i = 0; $i < $start->getdim(0); $i++) {
        # Start is zero-based and End is one-based, just like BED format
        my $start_coordinate = $start($i)->sclr;
        # The end position is really the start of the gap
        my $end_coordinate = $end($i)->sclr - 1;
        my $cluster_pdl = $chr_pdl($start_coordinate:$end_coordinate)->sever;
        my ($gt_idx,$zero_idx) = which_both($cluster_pdl >= $self->min_depth);
        # TODO: This can be removed once everything is validated
        if ($zero_idx->getdim(0)) {
            die('Unexpected number in PDL: '. $start_coordinate ."\t". $end_coordinate ."\t". $cluster_pdl);
        }
        my $start_pos = $start_coordinate + $offset;
        my $end_pos = $end_coordinate + $offset + 1;
        push @clusters, [ $start_pos, $end_pos, $cluster_pdl];
    }
    return @clusters;
}


sub generate_gaps {
    my $self = shift;
    my $offset = shift || 0;

    my $chr_pdl = $self->{_coverage_pdl};
    unless (defined($chr_pdl)) { return; }
    my ($quantity,$value) = rle($chr_pdl < $self->min_depth);
    my ($padded_quantity,$padded_value) = $quantity->where($value,$quantity!=0);
    my $padded_quantity_cumsum = $padded_quantity->cumusumover;

    my ($chr_lt_min_depth_idx,$chr_gt_min_depth_idx) = which_both($padded_value);
    # No gaps across entire region
    unless ($chr_lt_min_depth_idx->nelem) {
        return;
    }
    my $lt_first = $chr_lt_min_depth_idx(0)->sclr;
    # The entire region is a gap
    unless ($chr_gt_min_depth_idx->nelem) {
        my $end = $padded_quantity_cumsum->index($chr_lt_min_depth_idx);
        my @gaps = [$lt_first, ( ($end(0)->sclr - 1) + $offset)];
        return @gaps;
    }

    my $gt_first = $chr_gt_min_depth_idx(0)->sclr;

    my $gt_last = $chr_gt_min_depth_idx(-1)->sclr;
    my $lt_last = $chr_lt_min_depth_idx(-1)->sclr;

    my $start = $padded_quantity_cumsum->index($chr_gt_min_depth_idx);
    my $end = $padded_quantity_cumsum->index($chr_lt_min_depth_idx);

    if ($gt_first > $lt_first) {
        # Coverage first because the first index is zero
        # Add the initial start position(ie. zero)
        $start = append(zeroes(1),$start);
    }
    if ($gt_last > $lt_last) {
        # Gap on the end, remove the last value
        $start = $start(0:-2);
    }
    my @gaps;
    for (my $i = 0; $i < $start->getdim(0); $i++) {
        # Start is zero-based and End is one-based, just like BED format
        my $start_coordinate = $start($i)->sclr;
        # The end position is really the start of the gap
        my $end_coordinate = $end($i)->sclr - 1;
        my $start_pos = $start_coordinate + $offset;
        my $end_pos = $end_coordinate + $offset;
        push @gaps, [ $start_pos, $end_pos ];
    }
    return @gaps;
}

1;
