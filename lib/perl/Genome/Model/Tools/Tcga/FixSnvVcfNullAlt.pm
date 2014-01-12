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

#TODO This should be replaced by Genome::File::Vcf::Reader/Writer, but
#currently a bug in sample metainfo fails to properly parse sample
#info if sample field contains more than one sample info like:
###SAMPLE <XXXX>,<xxxxx>,<xxxx>
sub execute {
    my $self = shift;
    my $input_file  = $self->input_file;
    my $output_file = $self->output_file;

    my $in_fh  = Genome::Sys->open_file_for_reading($input_file)  or die "Failed to open $input_file for reading";
    my $out_fh = Genome::Sys->open_file_for_writing($output_file) or die "Failed to open $output_file for writing";

    while (my $line = $in_fh->getline) {
        if ($line =~ /^#/) {
            $out_fh->print($line);
        }
        else {
            my @columns = split /\t/, $line;
            if ($columns[4] eq '.') {
                $columns[4] = 'N';
                $line = join "\t", @columns;
            }
            $out_fh->print($line);
        }
    }
    $in_fh->close;
    $out_fh->close;

    return 1;
}

1;


