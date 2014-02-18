package Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory;

use strict;
use warnings;

use Genome;

use Genome::Utility::IO::SeparatedValueWriter;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriterFactory { 
    is => 'UR::Singleton',
};

sub availble_fields {
    return (qw/ chromosome position alleles id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /);
}

sub build_writer {
    my ($class, $output_params_string) = @_;

    my %defaults = (
        output => '-',
        format => 'csv', # vcf
    );
    my %csv_defaults = (
        separator => "\t",
        headers => 'chromosome,position,alleles',
        print_headers => 1,
        in_place_of_null_value => 'NA',
        ignore_extra_columns => 1,
    );

    my %output_params;
    if ( not $output_params_string ) {
        $output_params{output} = $defaults{output};
    }
    else {
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
    }

    my $format = delete $output_params{format};
    $format = $defaults{format} if not $format;
    my $writer_class;
    if ( $format eq 'vcf' ) {
        $class->error_message("VCF NOT IMPLEMENTED"); return;
        $writer_class = 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf';
    }
    elsif ( $format eq 'csv' ) {
        $writer_class = 'Genome::Utility::IO::SeparatedValueWriter';
        for my $default ( keys %csv_defaults ) {
            if ( not exists $output_params{$default} ) {
                $output_params{$default} = $csv_defaults{$default};
            }
        }
        $output_params{separator} = "\t" if $output_params{separator} =~ /^tab$/i;
        $output_params{headers} = [ split(',', $output_params{headers}) ];
    }
    else {
        $class->error_message('Unknown output format! '.$output_params{format});
        return;
    }

    my $writer = $writer_class->create(%output_params);
    if ( not $writer ) {
        $class->error_message('Failed to open writer! Params: '.Data::Dumper::Dumper(\%output_params));
        return;
    }

    return $writer;
}

1;

