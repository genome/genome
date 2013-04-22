
package Genome::Model::Tools::Picard::IlluminaBasecallsToSam;

# http://picard.sourceforge.net/command-line-overview.shtml#IlluminaBasecallsToSam

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Model::Tools::Picard::IlluminaBasecallsToSam {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        basecalls_directory => {
            is  => 'String',
            doc => 'The basecalls output directory.',
            },
        lane => {
            is  => 'String',
            doc => 'Lane number.',
        },
        run_barcode => {
            is  => 'String',
            doc => 'Prefixed to read names.',
        },
        read_group_id => {
            is  => 'String',
            doc => 'ID used to link RG header record with RG tag in SAM record.',
            is_optional => 1,
        },
        sequencing_center => {
            is  => 'String',
            doc => 'The name of the sequencing center that produced the reads to fill in the RG.CN tag.',
        },
        run_start_date => {
            is  => 'String',
            doc => 'The start date of the run.',
            is_optional => 1,
        },
        platform => {
            is  => 'String',
            doc => 'The name of the sequencing technology that produced the read.',
        },
        read_structure => {
            is  => 'String',
            doc => 'A description of the logical structure of clusters in an Illumina Run',
        },
        library_params => {
            is  => 'String',
            doc => 'Tab-separated file for creating all output BAMs for a run.',
        },
        adapters_to_check => {
            is  => 'String',
            doc => 'Which adapters to look for in the read.',
            is_optional => 1,
        },
        num_processors => {
            is  => 'String',
            doc => 'Run this many threads in parallel.',
            is_optional => 1,
        },
        first_tile => {
            is  => 'String',
            doc => 'If set, this is the first tile to be processed (for debugging).',
            is_optional => 1,
        },
        tile_limit => {
            is  => 'String',
            doc => 'If set, process no more than this many tiles (for debugging).',
            is_optional => 1,
        },
        force_gc => {
            is  => 'String',
            doc => 'If true, call System.gc() periodically.',
            is_optional => 1,
        },
        max_read_in_ram_per_tile => {
            is  => 'String',
            doc => 'Configure SortingCollections to store this many records before spilling to disk.',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Generate a SAM or BAM file from data in an Illumina basecalls output directory.';
}

sub help_detail {
    return <<EOS
    Generate a SAM or BAM file from data in an Illumina basecalls output directory.  For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#IlluminaBasecallsToSam
EOS
}

sub execute {
    my $self = shift;

    my $jar_path = $self->picard_path . '/IlluminaBasecallsToSam.jar';
    unless (-e $jar_path) {
        die('Failed to find '. $jar_path .'!  This command may not be available in version '. $self->use_version);
    }

    my $args = join(' ',
        map {
            my $value = $self->$_;
            defined($value) ? (uc($_) . "='$value'") : ()
        } sort qw(
            basecalls_directory
            lane
            run_barcode
            read_group_id
            sequencing_center
            run_start_date
            platform
            read_structure
            library_params
            adapters_to_check
            num_processors
            first_tile
            tile_limit
            force_gc
            max_read_in_ram_per_tile
        )
    );

    my $cmd = $jar_path . " net.sf.picard.sam.IlluminaBasecallsToSam $args";
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [ $self->basecalls_directory, $self->library_params, ],
#        output_files => [ $self->output ],
#        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
__END__

