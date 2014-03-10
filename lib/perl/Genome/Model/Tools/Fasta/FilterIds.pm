package Genome::Model::Tools::Fasta::FilterIds;
use strict;
use warnings;
use Genome;
use Bio::SeqIO;

class Genome::Model::Tools::Fasta::FilterIds {
    is => 'Command::V2',
    has_input => [
        input_filename => {
            is => 'FilesystemPath',
            shell_args_position => 1,
            doc => 'the input file',
        },
        output_filename =>  {
            is => 'FilesystemPath',
            shell_args_position => 2,
            doc => 'the path to the file that will be created',
        },
        whitelist_regex => {
            is => 'Text',
            is_optional => 1,
            doc => 'include only IDs that match this pattern',
        },
        blacklist_regex => {
            is => 'Text',
            is_optional => 1,
            doc => 'exclude any IDs that match this pattern',
        },
    ],
    has_param => [
        verbose => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'more messages'
        },          
    ],
    doc => "filter sequences from a fasta file based on patterns applied to the IDs",
};

sub execute {
    my $self = shift;
    my $input_filename = $self->input_filename;
    my $output_filename = $self->output_filename;
    my $verbose = $self->verbose;

    my $blacklist_regex = $self->blacklist_regex;
    my $whitelist_regex = $self->whitelist_regex;

    my $reader = Bio::SeqIO->new( '-file' => '< '.$input_filename, '-format' => 'fasta');
    my $writer = Bio::SeqIO->new( '-file' => '> '.$output_filename, '-format' => 'fasta');
    
    while (my $seq = $reader->next_seq) { 
        my $id = $seq->id;
        if ($blacklist_regex and $id =~ $blacklist_regex) {
            $self->debug_message("skipping $id because it matches the blacklist pattern");
            next;
        }
        elsif ($whitelist_regex and not $id =~ $whitelist_regex) {
            $self->debug_message("skipping $id because it does not match the whitelist pattern");
            next;
        }
        elsif ($verbose) {
            $self->debug_message("keeping $id");
        }
        $writer->write_seq($seq);
    }
    return 1;
}

sub help_synopsis {
return <<'EOS'
    gmt fasta filter-ids in.fa out.fa --whitelist '^(\d+|X,Y)$' --blacklist '6'
EOS
}

sub help_detail {
    return <<EOS
This tool filters a FASTA sequence file, removing entries based on the ID in the FASTA header.

If the "whitelist regex" (-w) option is supplied, only IDs that match this regular expression will be included.

If the "blacklist regex" (-b) option is supplied, only IDs that do NOT match this regular expression will be included.  The blacklist takes precedence.
EOS
}

1;

