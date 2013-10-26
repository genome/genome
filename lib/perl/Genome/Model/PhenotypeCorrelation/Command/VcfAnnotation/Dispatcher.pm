package Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Dispatcher;

use Carp qw/confess/;
use Data::Dumper;
use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Dispatcher {
    is => "Genome::Command::Base",
    has_input => [
        species_name => {
            is => "Text",
            default_value => "human",
            doc => "Species name",
        },
        strategy => {
            is => "Text",
            doc => "The annotation strategy to apply",
        },
        input_file => {
            is => "Text",
            is_optional => 1,
            doc => "Input vcf file",
        },
        output_file => {
            is => "Text",
            is_optional => 1,
            is_output => 1,
            doc => "Output vcf file",
        },
    ]
};

sub help_synopsis {
    return <<EOS

    genome model phenotype-correlation vcf-annotation dispatcher \\
        --input-file in.vcf \\
        --output-file out.vcf \\
        --strategy 'vep {condel: b, sift: b, polyphen: b, hgnc: 1} | joinx [{source_name: dbsnp, source_version: 137, info_fields: GMAF}, {source_name: 1kg-wgs, source_version: 20101123}]'
EOS
}

sub help_detail {
    return <<EOS
Annotate a vcf file with vep and/or joinx vcf-annotate.
EOS
}


sub _tool_class_name {
    my $name = shift;
    my $leaf = join("", map {ucfirst} split("-", $name));
    my @package = split("::", __PACKAGE__);
    pop @package;
    my $package_base = join("::", @package);
    return sprintf("%s::%s", $package_base, $leaf);
}

sub _create_command {
    my ($self, $cmd, $input_file, $output_file) = @_;

    my $pkg = _tool_class_name($cmd->{type});

    my %params = (
        params => $cmd->{params},
        input_file => $input_file,
        output_file => $output_file,
        species_name => $self->species_name,
        );

    my $obj = $pkg->create(%params);
    return $obj;
}

sub execute {
    my $self = shift;
    my $strategy = Genome::Model::PhenotypeCorrelation::Command::VcfAnnotation::Strategy->create(
            strategy => $self->strategy);

    my $tree = $strategy->execute();
    my $input_file = $self->input_file;

    for my $idx (0..$#$tree) {
        my $item = $tree->[$idx];
        my $output_file;
        if ($idx == $#$tree) {
            $output_file = $self->output_file;
        }
        else {
            $output_file = Genome::Sys->create_temp_file_path;
        }

        $self->status_message("Executing $item->{type}, $input_file -> $output_file");
        my $cmd = $self->_create_command($item, $input_file, $output_file);
        if (!$cmd->execute) {
            confess "Failed to execute command: " . Dumper($item);
        }
        $input_file = $output_file;
    }

    $self->status_message("Annotation complete.");
    return 1;
}
