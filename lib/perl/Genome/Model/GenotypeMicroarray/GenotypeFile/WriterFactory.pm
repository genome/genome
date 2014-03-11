package Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Writer;
use Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory { 
    is => 'UR::Singleton',
};

sub build_writer {
    my ($class, %params) = @_;

    my $header = delete $params{header};
    if ( not $header or not $header->isa('Genome::File::Vcf::Header') ) {
        $class->error_message('No header given to create writer!');
        return;
    }

    my $string = delete $params{string};
    if ( %params ) {
        $class->error_message('Unknwon params sent to build writer! '.Data::Dumper::Dumper(\%params));
        return;
    }

    my $writer_params = $class->_parse_params_string($string);
    return if not $writer_params;

    $writer_params->{header} = $header;
    my $writer = $class->_build_writer($writer_params);
    if ( not $writer ) {
        $class->error_message('Failed to build writer! Params: '.Data::Dumper::Dumper($writer_params));
        return;
    }

    return $writer;
}

sub _parse_params_string {
    my ($class, $string) = @_;

    my @writer_config_tokens = split(':', $string);
    my %writer_params;
    if ( @writer_config_tokens == 1 and $writer_config_tokens[0] !~ /=/ ) {
        $writer_params{output} = $writer_config_tokens[0];
    }
    else { 
        for my $writer_config_token ( @writer_config_tokens ) {
            my ($key, $value) = split('=', $writer_config_token);
            if ( not defined $value ) {
                $writer_params{output} = $key;
                next;
            }
            if ( exists $writer_params{$key} ) {
                $class->error_message('Duplicate attribute in writer config: '.$key);
                return;
            }
            $writer_params{$key} = $value;
        }
    }

    # Use STDOUT if no output
    $writer_params{output} ||= '-' if not $writer_params{file};

    # Resolve format
    if ( not $writer_params{format} ) { 
        # use suffix
        my ($suffix) = $writer_params{output} =~ /\.(\w+)$/;
        if ( $writer_params{separator} ) {
            $writer_params{format} = 'csv';
        }
        elsif ( not $suffix ) {
            $writer_params{format} = 'vcf';
        }
        elsif ( $suffix =~ /tsv/g ) {
            $writer_params{format} = 'csv';
        }
        else {
            $writer_params{format} = $suffix;
        }
    }

    if ( not grep { $writer_params{format} eq  $_ } (qw/ csv vcf /) ) {
        $class->error_message('Unknown format to write! '. $writer_params{format});
        return;
    }

    return \%writer_params;
}

sub _build_writer {
    my ($class, $writer_params) = @_;

    my $format = delete $writer_params->{format};

    my $method = '_build_'.$format.'_writer';
    return $class->$method($writer_params);
}

sub _build_vcf_writer {
    my ($class, $writer_params) = @_;
    return Genome::File::Vcf::Writer->new($writer_params->{output}, $writer_params->{header});
}

sub _build_csv_writer {
    my ($class, $writer_params) = @_;
    return Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv->create(%$writer_params);
}

1;

