package Genome::Model::Tools::Blat::MergeAmplicons;

use strict;
use warnings;

use Genome;
use Genome::Utility::FastaStreamIn;
use Genome::Utility::FastaStreamOut;

class Genome::Model::Tools::Blat::MergeAmplicons {
    is => 'Command',
    has => [
            amplicon_files => {
                               doc => 'the amplicon files to merge',
                               is => 'String',
                               is_input => 1,
                               is_many => 1,
                       },
            output_file => {
                             doc => 'the output file for the merged amplicons',
                             is => 'String',
                             is_output => 1,
                         },
        ],
};

sub help_brief {
    "",
}

sub help_synopsis {
    return <<"EOS"

EOS
}

sub help_detail {
    return <<EOS

EOS
}

sub execute {
    my $self = shift;

    my %amplicons;
    for my $amplicon_file ($self->amplicon_files) {
        my $fh = IO::File->new($amplicon_file,'r');
        my $reader = Genome::Utility::FastaStreamIn->new($fh);
        while (my $header = $reader->next_header) {
           $header =~ /^>(\S+)/;
           my $amplicon_name = $1;
           unless (defined $amplicons{$amplicon_name}) {
               $amplicons{$amplicon_name} = $header;
           } else {
               if ($header ne $amplicons{$amplicon_name}) {
                   my $error_message = 'Different headers found for '. $amplicon_name .":\n";
                   $error_message .= "\t". $header ."\n";
                   $error_message .= "\t". $amplicons{$amplicon_name} ."\n";
                   $self->error_message($error_message);
                   return;
               }
           }
       }
    }
    my $fh = IO::File->new($self->output_file,'w');
    my $writer = Genome::Utility::FastaStreamOut->new($fh);
    for my $amplicon (keys %amplicons) {
        $writer->print_header($amplicons{$amplicon});
    }
    return 1;
}


1;
