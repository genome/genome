package Genome::Db::Cosmic::Command::ImportVcf;

use strict;
use warnings;
use Genome;

class Genome::Db::Cosmic::Command::ImportVcf {
    is => 'Command',
    has => [
        urls => {
            is => 'Text',
            is_many => 1,
        },
        output_file => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;
    my $temp_dir = Genome::Sys->create_temp_directory;
    my $count = 0;
    for my $url ($self->urls) {
        `wget -O $temp_dir/$count.vcf.gz --no-check-certificate $url`;
        $count++;
    }
    `gunzip $temp_dir/*.gz`;

    Genome::Model::Tools::Joinx::Sort->execute(
        input_files => [map{"$temp_dir/$_.vcf"} (0..$count-1)],
        output_file => $self->output_file,
    );
    return 1;
}

1;

