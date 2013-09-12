package Genome::File::Vcf;

use strict;
use warnings;
use Genome;

use Genome::Utility::Vcf qw(get_vcf_header get_samples_from_header parse_vcf_line);

class Genome::File::Vcf {
    is => 'Genome::File::Base',
    has_transient_optional => {
        _file_handle => {
            is => 'String',
        },
        sample_names => {
            is => 'String',
            is_many => 1,
            doc => 'Sample names from the header',
        },
        header => {
            is => 'String',
            doc => 'The entire header',
        },
    },
};

sub open {
    my $self = shift;
    
    my ($header, $fh) = get_vcf_header($self->path);
    $self->header($header);

    my @samples = get_samples_from_header($header);
    $self->sample_names(\@samples);

    $self->_file_handle($fh);
    return 1;
}

sub next {
    my $self = shift;

    my $line = $self->_file_handle->getline;
    return unless $line;
    return parse_vcf_line($line, [$self->sample_names]);
}

1;

