package Genome::Model::Tools::454::Sfffile;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::454::Sfffile {
    is => ['Genome::Model::Tools::454'],
    has => [
            in_sff_files => {
                            doc => 'The sff file to operate',
                     },
            out_sff_file => {
                            is => 'string',
                            doc => 'The output file path',
                        },
            params => {
                       is => 'string',
                       doc => 'The params to pass to sfffile',
                       is_optional => 1,
		   },
        ],
};

sub help_brief {
    "constructs a single SFF file containing the reads from a list of SFF files and/or 454 runs"
}

sub help_detail {
    return <<EOS
see 'sfffile' usage for valid params
EOS
}

sub execute {
    my $self = shift;

    my @in_sff_files = @{$self->in_sff_files};
    my $out_sff_file = $self->out_sff_file;
    my $params = $self->params || '';
    $params .= ' -o '. $out_sff_file;
    my $cmd = $self->bin_path .'/sfffile '. $params .' '. join(' ',@in_sff_files);
    Genome::Sys->shellcmd(
                                         cmd => $cmd,
                                         input_files => \@in_sff_files,
                                         output_files => [$out_sff_file],
                                     );
    return 1;
}

1;


