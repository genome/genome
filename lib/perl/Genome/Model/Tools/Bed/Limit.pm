package Genome::Model::Tools::Bed::Limit;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bed::Limit {
    is => ['Command'],
    has => [
        gene_list => {
            is => 'Text',
            is_optional => 1,
        },
        input_bed_file => {
            is => 'Text',
        },
        output_bed_file => {
            is => 'Text',
        },
        feature_types => {
            is => 'Text',
            doc => 'A comma delimited list of features that can include: intron,cds_exon,utr_exon,rna',
            default_value => 'cds_exon',
        },
    ],
};

sub execute {
    my $self = shift;

    my %feature_types;
    if ($self->feature_types =~ /,/) {
        %feature_types = map { $_ => 1 } split(',',$self->feature_types);
    } else {
        $feature_types{$self->feature_types} = 1;
    }
    my $input_fh = Genome::Sys->open_file_for_reading($self->input_bed_file);
    my $output_fh = Genome::Sys->open_file_for_writing($self->output_bed_file);

    my %include;
    if ($self->gene_list) {
        my $list_fh = Genome::Sys->open_file_for_reading($self->gene_list);
        while (my $line = $list_fh->getline) {
            chomp($line);
            if ($line =~ /^(\S+)/) {
                $include{$1} = 1;
            }
        }
        $list_fh->close;
    }
    while (my $line = $input_fh->getline) {
        chomp($line);
        my @entry = split("\t",$line);
        my $name = $entry[3];
        my ($gene,$transcript,$feature_type,$ordinal) = split(':',$name);
        if ($feature_types{$feature_type}) {
            if ($self->gene_list) {
                if ($include{$gene}) {
                    print $output_fh $line ."\n";
                }
            } else {
                print $output_fh $line ."\n";
            }
        }
    }
    $input_fh->close;
    $output_fh->close;
    return 1;
}

1;
