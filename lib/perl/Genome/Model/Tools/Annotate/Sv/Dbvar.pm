package Genome::Model::Tools::Annotate::Sv::Dbvar;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Dbvar {
    is => "Genome::Model::Tools::Annotate::Sv::Base",
    doc => "Annotate overlaps with dbvar SVs",
    has_input => [
        annotation_file => {
            is => 'Text',
            doc => 'File containing UCSC table',
            example_values => ["/gsc/scripts/share/BreakAnnot_file/human_build37/GRCh37.remap.all.germline.ucsc.gff"],
        },
        overlap_fraction => {
            is => 'Number',
            doc => 'Fraction of overlap (reciprocal) required to hit',
            default => 0.5,
        },
    ],
};

sub help_detail {
    return "Determines whether the SV breakpoints match a dbVar SV within some distance.  It also checks to see that the SV and the dbVar SV reciprocally overlap each other by a given fraction.";
}

sub process_breakpoint_list {
    my ($self, $breakpoints_list) = @_;
    my %output;
    my $dbvar_annotation = $self->read_dbvar_annotation($self->annotation_file);
    $self->annotate_interval_overlaps($breakpoints_list, $dbvar_annotation, "dbvar_annotation");
    foreach my $chr (keys %{$breakpoints_list}) {
        foreach my $item (@{$breakpoints_list->{$chr}}) {
            my $key = $self->get_key_from_item($item);
            $output{$key} = [$self->get_var_annotation($item, $item->{dbvar_annotation}->{bpB})];
        }
    }
    return \%output;
}

sub column_names {
    return ('dbvar');
}

sub read_dbvar_annotation {
    my $self = shift;
    my $file = shift;

    my %annotation;
    my $in = Genome::Sys->open_file_for_reading($file) || die "Unable to open annotation: $file\n";

    while (my $line = <$in>) {
       chomp $line;
       next if ($line =~ /^\#/);
       my $p;
       my @fields = split /\t+/, $line;
       $p->{chrom} = $fields[0];
       $p->{bpA} = $fields[3];
       $p->{bpB} = $fields[4];
       $p->{name} = $fields[8];
       $p->{chrom} =~ s/chr//;
       push @{$annotation{$p->{chrom}}}, $p;
    }
    return \%annotation;
}

1;

