package Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory { 
    is => 'UR::Singleton',
};

sub default_format {
    return 'vcf';
}

sub build_writer {
    my ($class, $output_params_string) = @_;

    my $output_params = $class->parse_output_params_string($output_params_string);
    return if not $output_params;

    my $writer_class = $class->_resolve_writer_class($output_params);
    return if not $writer_class;

    my $writer = eval{ $writer_class->create(%$output_params); };
    if ( not $writer ) {
        $class->error_message('Failed to open writer! Params: '.Data::Dumper::Dumper($output_params));
        return;
    }

    return $writer;
}

sub parse_output_params_string {
    my ($class, $output_params_string) = @_;

    my %output_params;
    if ( not $output_params_string ) { # do the default thing
        return \%output_params;
    }

    my @output_config_tokens = split(':', $output_params_string);
    if ( @output_config_tokens == 1 and $output_config_tokens[0] !~ /=/ ) {
        $output_params{output} = $output_config_tokens[0];
    }
    else { 
        for my $output_config_token ( @output_config_tokens ) {
            my ($key, $value) = split('=', $output_config_token);
            if ( not defined $value ) {
                $output_params{output} = $key;
                next;
            }
            if ( exists $output_params{$key} ) {
                $class->error_message('Duplicate attribute in output config: '.$key);
                return;
            }
            $output_params{$key} = $value;
        }
    }

    return \%output_params;
}

sub _resolve_writer_class {
    my ($class, $output_params) = @_;

    my $format = delete $output_params->{format} || default_format();
    my %writer_classes = (
        csv => 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv',
        vcf => 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf',
    );
    my $writer_class = $writer_classes{$format};
    if ( not $writer_class ) {
        $class->error_message('Unknown output format! '.$format);
        return;
    }

    return $writer_class;
}

1;

