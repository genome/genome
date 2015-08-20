package Genome::VariantReporting::Suite::Vep::TestHelper;

use Exporter 'import';
@EXPORT_OK = qw(create_vcf_header create_entry_with_vep);

use Genome::File::Vcf::Reader;

sub create_vcf_header {
    my $csq_format = shift;

    my $header_txt = <<EOS;
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="Passed all filters">
##FILTER=<ID=BAD,Description="This entry is bad and it should feel bad">
##INFO=<ID=A,Number=1,Type=String,Description="Info field A">
##INFO=<ID=C,Number=A,Type=String,Description="Info field C (per-alt)">
##INFO=<ID=E,Number=0,Type=Flag,Description="Info field E">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: $csq_format">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Filter">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
EOS
    my @lines = split("\n", $header_txt);
    my $header = Genome::File::Vcf::Header->create(lines => \@lines);
    return $header
}

sub create_entry_with_vep {
    my $vep_string = shift;
    my $csq_format = shift;

    my $info = join(';', 'A=B', 'C=8,9', 'E', $vep_string);
    my @fields = (
        '1',            # CHROM
        10,             # POS
        '.',            # ID
        'A',            # REF
        'C,G',          # ALT
        '10.3',         # QUAL
        'PASS',         # FILTER
        $info,          # INFO
        "GT:DP",        # FORMAT
        "0/1:12"        # S1
    );
    my $entry_txt = join("\t", @fields);

    my $entry = Genome::File::Vcf::Entry->new(create_vcf_header($csq_format), $entry_txt);
    return $entry;
}

1;
