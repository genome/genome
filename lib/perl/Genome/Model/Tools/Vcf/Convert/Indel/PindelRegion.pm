package Genome::Model::Tools::Vcf::Convert::Indel::PindelRegion;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Writer;
use Genome::File::Vcf::Reader;


class Genome::Model::Tools::Vcf::Convert::Indel::PindelRegion {
    is =>  'Genome::Model::Tools::Vcf::Convert::Indel::Pindel',
};


sub execute {
    my $self = shift;
    $self->SUPER::_execute_body(@_);

    my $tmp_output = Genome::Sys->create_temp_file_path;
    Genome::Sys->move_file($self->output_file, $tmp_output);

    my $in_vcf = Genome::File::Vcf::Reader->new($tmp_output);
    my $header = $in_vcf->header;

    my $FT_header_str = '<ID=FT,Number=1,Type=String,Description="Sample genotype filter">';
    $header->add_format_str($FT_header_str);

    my $out_vcf = Genome::File::Vcf::Writer->new($self->output_file, $header);

    while (my $entry = $in_vcf->next) {
        $entry->add_format_field('FT');
        if ($entry->sample_field(0, 'GT') eq '0/0' and $entry->sample_field(1, 'GT') =~ /1/) { #Somatic
            map{$entry->set_sample_field($_, 'FT', 'PASS')}qw(0 1);
        }
        else {
            map{$entry->set_sample_field($_, 'FT', '.')}qw(0 1);
        }
        $out_vcf->write($entry);
    }
    $out_vcf->close;

    return 1;
}

1;

