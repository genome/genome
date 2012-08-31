package Genome::Model::Tools::DetectVariants2::Combine::UnionIndel;


use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Combine::UnionIndel{
    is => 'Genome::Model::Tools::DetectVariants2::Combine',
};


sub help_brief {
    "Union two indel variant bed files",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants combine union indel --variant-file-a samtools.hq.v1.bed --variant-file-b varscan.hq.v1.bed --output-file 
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

sub _variant_type { 'indels' };

1;
