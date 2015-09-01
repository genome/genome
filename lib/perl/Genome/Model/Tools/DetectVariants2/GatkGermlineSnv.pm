package Genome::Model::Tools::DetectVariants2::GatkGermlineSnv;

use strict;
use warnings;

use Cwd;

use Genome;

class Genome::Model::Tools::DetectVariants2::GatkGermlineSnv{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has_constant => [
        detect_snvs => { value => 1 },
        detect_svs => {},
        detect_indels => {},
    ],
    has => [
        mb_of_ram => {
            is => 'Text',
            doc => 'Amount of memory to allow GATK to use',
            default => 5000,
        },
    ],
    has_param => [
         lsf_resource => {
             default_value => Genome::Config::get('lsf_resource_dv2_gatk'),
         },
     ],
};

sub _detect_variants {
    my $self = shift;

    die $self->error_message('The GATK GermlineSnv detector is no longer available');
}

sub has_version {
    my $self = shift;

    return;
}

#TODO clean all of this up. It is usually/should be based on logic from Genome::Model::Tools::Bed::Convert logic in process_source... 
# this should be smarter about using that work ... perhaps process_source should call a method that just parses one line, and this method can be replaced by a call to that instead
sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    # Skip header lines
    if ($line =~ /^#/) {
        return;
    }

    my ($chromosome, $start, $id, $reference, $variant) = split("\t", $line);
    my $stop;
    if(length($reference) == 1 and length($variant) == 1) {
        #SNV case
        $stop = $start;
        $start -= 1; #convert to 0-based coordinate
    } else {
        die $class->error_message("Unhandled variant type encountered for line: $line");
    }

    unless (defined $chromosome && defined $stop && defined $reference && defined $variant) {
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return [$chromosome, $stop, $reference, $variant];
}

1;
