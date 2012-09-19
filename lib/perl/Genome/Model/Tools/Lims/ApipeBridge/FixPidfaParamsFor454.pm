package Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsFor454;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsFor454 { 
    is => 'Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForBase', 
};

sub instrument_data_type { return '454'; }
sub valid_prior_processes { return ( 'analyze 454 output', 'analyze 454 run', 'analyze 454 region', 'demux 454 region' ); }

sub _get_sequence_item_from_prior {
    my ($self, $prior) = @_;

    my @run_regions = $prior->get_454_run_regions;
    if ( not @run_regions ) {
        $self->error_message('No 454 run regions for prior PSE!');
        return;
    }
    my @region_ids = map {$_->region_id() } @run_regions;
    my @region_indexes = GSC::RegionIndex454->get(region_id => \@region_ids);
    if ( not @region_indexes ) {
        $self->error_message('No 454 region indexes for prior PSE!');
        return;
    }
    elsif ( @region_indexes > 1 ) {
        $self->status_message('More than one 454 region index found for prior PSE. If you know the instrument data id, use that as a staring point.');
        return;
    }

    return $region_indexes[0];
}

sub _additional_params_to_fix {
    my ($self, $params_to_fix, $sequence_item) = @_;

    my $sff_file;
    if ( $sequence_item->index_sequence ) { # sff from region index
        $sff_file = $sequence_item->get_index_sff;
    }
    else { # get sff from region, it may not exist b/c they used to be in the db
        my $run_region = $sequence_item->get_run_region;
        my $sff_file_path = eval{ $run_region->sff_filesystem_location; };
        $sff_file = $sff_file_path->absolute if $sff_file_path;
    }

    print "SFF $sff_file\n";
    $params_to_fix->{sff_file} = $sff_file if $sff_file and -e $sff_file;

    return 1;
}

1;

