package Genome::Model::Tools::DetectVariants2::Combine::UnionuniqueIndel;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Combine::UnionuniqueIndel{
    is => 'Genome::Model::Tools::DetectVariants2::Combine',
    doc => 'Union indels into one file',
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants combine unionunique-indel --variant-file-a samtools.hq.v1.bed --variant-file-b varscan.hq.v1.bed --output-file 
EOS
}

sub _variant_type { 'indels' };

1;
