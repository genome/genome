package Genome::Model::Tools::Vcf::CleanupVcf;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Vcf::CleanupVcf {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            is_optional => 0,
            doc => "The file to output",
        },
        input_file => {
            is => 'Text',
            is_optional => 0,
            is_input => 1,
            doc => 'VCF file to remove tagged columns from',
        },
    ],
};

sub help_synopsis {
    <<'HELP';
Remove per-caller columns from a vcf.  The per-caller columns should be tagged with square brackets around the caller name (example: F_AKE_SAMPLENAME000-[Samtools])
HELP
}

sub execute {
    my $self = shift;
    unless(-s $self->input_file){
        $self->error_message("Input VCF: " . $self->input_file . " does not exist or has no size"); 
    }

    my $ifh = Genome::Sys->open_gzip_file_for_reading($self->input_file);
    my $ofh = Genome::Sys->open_gzip_file_for_writing($self->output_file);

    my @columns_to_print = ();
    while (my $line = $ifh->getline){
        chomp $line;
        if($line =~ /^##/){
            print $ofh $line."\n";
        }elsif($line =~ /#\w+/){
            my @headers = split("\t", $line);
            for (my $i = 0; $i < scalar @headers; $i++){
                push(@columns_to_print, $i) unless $headers[$i] =~ /\[.*\]/;
            }
            print $ofh $self->_format_output_line($line, @columns_to_print) . "\n";
        }else{
            print $ofh $self->_format_output_line($line, @columns_to_print) . "\n";
        }
    }

    $ifh->close;
    $ofh->close;

    return 1;
}

sub _format_output_line {
    my ($self, $line, @columns_to_print) = @_;
    my @columns = split("\t", $line);
    return join("\t", map($columns[$_], @columns_to_print));
}

1;
