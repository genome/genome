package Genome::Model::Tools::DetectVariants2::BamToCna;

use strict;
use warnings;

use Cwd;

use Genome;

class Genome::Model::Tools::DetectVariants2::BamToCna{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has => [
        ratio => {
            type => 'Number',
            is_optional => 1,
            default => 0.25,
            doc => 'Ratio diverged from median, used to find copy number neutral region (default = 0.25).'
        },
        window_size => {
            type => 'Number',
            is_optional => 1,
            default => 10000,
            doc => 'Window size (bp) for counting reads contributing to copy number in that window (resolution, default = 10000 bp).'
        },
        normalize_by_genome => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 1,
            doc => 'enable this flag to normalize by the whole genome median.',
        },
        params => { #calculated for SoftwareResult for now--in the future should support this fully
            is_optional => 1,
            calculate_from => ['ratio', 'window_size', 'normalize_by_genome'],
            calculate => q{
                my $params  = '';
                if(defined $ratio)                  { $params .= ' --ratio ' . $ratio; }
                if(defined $window_size)            { $params .= ' --window-size ' . $window_size; }
                if(defined $normalize_by_genome)    { $params .= ' --normalize-by-genome ' . $normalize_by_genome; }
                return $params;
            },
        }
    ],
};

sub _detect_variants {
    my $self = shift;
    my $b2c_cmd = Genome::Model::Tools::Somatic::BamToCna->create( 
        tumor_bam_file => $self->aligned_reads_input, 
        normal_bam_file => $self->control_aligned_reads_input,
        output_file => $self->_temp_staging_directory."/cnvs.hq",
        ratio => $self->ratio,
        window_size => $self->window_size,
        normalize_by_genome => $self->normalize_by_genome,
    );
    unless($b2c_cmd->execute){
        $self->error_message("Failed to run BamToCna command.");
        die $self->error_message;
    }

    return 1;
}

sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

# Cna current does need to sort its output (header lines being the primary problem, can probably remove them soon)
sub _sort_detector_output {
    return 1;
}

1;
