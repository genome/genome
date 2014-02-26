package Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory { 
    is => 'UR::Singleton',
};

sub build_writer {
    my ($class, $writer_params_string) = @_;

    my $writer_params = $class->_parse_writer_params_string($writer_params_string);
    return if not $writer_params;

    my $writer = $class->_build_writer($writer_params);
    if ( not $writer ) {
        $class->error_message('Failed to build writer! Params: '.Data::Dumper::Dumper($writer_params));
        return;
    }

    return $writer;
}

sub _parse_writer_params_string {
    my ($class, $writer_params_string) = @_;

    my %writer_params;
    if ( not $writer_params_string ) { # do the default thing
        return \%writer_params;
    }

    my @writer_config_tokens = split(':', $writer_params_string);
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

    return \%writer_params;
}

sub _build_writer {
    my ($class, $writer_params) = @_;

    my $format = delete $writer_params->{format} || 'vcf';
    if ( not grep { $format eq  $_ } (qw/ csv vcf /) ) {
        $class->error_message('Unknown format to write! '. $format);
        return;
    }

    $writer_params->{output} ||= '-' if not $writer_params->{file};

    my $method = '_build_'.$format.'_writer';
    return $class->$method($writer_params);
}

sub _build_vcf_writer {
    my ($class, $writer_params) = @_;
    my $header = Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->header;
    return Genome::File::Vcf::Writer->new( $writer_params->{output}, $header );
}

sub _build_csv_writer {
    my ($class, $writer_params) = @_;

    $writer_params->{separator} = "\t" if not $writer_params->{separator} or $writer_params->{separator} =~ /^tab$/i;
    $writer_params->{print_headers} = delete $writer_params->{headers} if exists $writer_params->{headers};
    $writer_params->{ignore_extra_columns} = 1;

    my $fields = delete $writer_params->{fields};
    if ( $fields ) {
        $writer_params->{headers} = [ split(',', $fields) ];
    }
    else {
        $writer_params->{headers} = [qw/ chromosome position alleles reference id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /];
    }

    return Genome::Utility::IO::SeparatedValueWriter->create(%$writer_params);
}

1;

