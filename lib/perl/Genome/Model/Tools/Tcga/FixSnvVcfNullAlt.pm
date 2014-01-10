package Genome::Model::Tools::Tcga::FixSnvVcfNullAlt;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;

class Genome::Model::Tools::Tcga::FixSnvVcfNullAlt {
	is  => 'Genome::Model::Tools::Tcga',
	has => [
		input_file  => {
			type => 'String',
			doc  => 'input vcf file',
		},
		output_file => {
			type => 'String',
			doc  => 'output vcf file',
		},
	],
};


sub execute {
    my $self = shift;
    my $input_file  = $self->input_file;
    my $output_file = $self->output_file;

    my $reader = Genome::File::Vcf::Reader->new($input_file);
    my $header = $reader->header;
    my $writer = Genome::File::Vcf::Writer->new($output_file, $header);

    while (my $entry = $reader->next) {
        unless (@{$entry->{alternate_alleles}}) {
            $entry->{alternate_alleles} = ['N'];
        }
        $writer->write($entry)
    }

    return 1;
}

1;


