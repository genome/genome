package Genome::Model::Tools::DetectVariants2::BreakdancerMax;

use warnings;
use strict;

use Genome;
use File::Basename;

class Genome::Model::Tools::DetectVariants2::BreakdancerMax{
    is => 'Genome::Model::Tools::DetectVariants2::Breakdancer',
};


sub chr_list_when_missing_idxstats {
    my $self = shift;
    die $self->error_message("chr list from samtools idxstats is empty, giving up");
}

sub _get_chr_list_from_idx_file {
    my ($self, $idxstats, $idx_file) = @_;

    if ($self->chromosome eq 'all') {
        # this is the _real_ <all>, since <all> is a legacy flag meaning only 1-22,X,Y,MT
        return @{$idxstats->map_ref_list($idx_file)};
    }
    else {
        my %map_chr_list = map { $_ => 1 } @{$idxstats->map_ref_list($idx_file)};
        if ($map_chr_list{$self->chromosome}) {
            return @{$self->chromosome};
        }
        else {
            die sprintf("Chromosome %s not found in BAM idxstats!",
                $self->chromosome);
        }
    }
}

1;
