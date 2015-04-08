package Genome::Model::Tools::Vcflib;

use strict;
use warnings;

use Genome;

my $DEFAULT = '1.0';

class Genome::Model::Tools::Vcflib {
    is  => 'Command::V2',
    has => [
        use_version => {
            is            => 'Version',
            is_optional   => 1,
            default_value => $DEFAULT,
            doc           => "Version of vcflib to use, default is $DEFAULT"
        },
    ],
};

sub help_brief {
    return "interface to a set of VCF-related tools that are part of vcflib";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools vcflib ...
EOS
}

sub help_detail {
    return <<EOS
    See the programs that are part of the vcflib project at

        https://github.com/ekg/vcflib.
    
    Version 1.0 represents the state of the github repository at commit
    "b1e9b31d2d95f957f86cdfd1e7c9ec25b4950ee8".
EOS
}

sub vcflib_tool_path {
    my ($self, $tool, $version) = @_;

    my $v = $version || $self->use_version;
    unless ($v) {
        die $self->error_message("no vcflib tool version specified!");
    }

    my $toolset = $self->get_version_toolset($v);
    unless (exists $toolset->{$tool}) {
        die $self->error_message(
            "no vcflib tool named '%s' exists in version '%s'!",
            $tool, $v
        );
    }

    my $path = $toolset->{$tool};
    return $path;
}

sub get_version_toolset {
    my ($self, $version) = @_;
    my %toolset = (
        '1.0' => {
            # core programs
            'vcf2dag'                    => '/usr/bin/vcf2dag1.0',
            'vcf2fasta'                  => '/usr/bin/vcf2fasta1.0',
            'vcf2tsv'                    => '/usr/bin/vcf2tsv1.0',
            'vcfaddinfo'                 => '/usr/bin/vcfaddinfo1.0',
            'vcfafpath'                  => '/usr/bin/vcfafpath1.0',
            'vcfallelicprimitives'       => '/usr/bin/vcfallelicprimitives1.0',
            'vcfaltcount'                => '/usr/bin/vcfaltcount1.0',
            'vcfannotate'                => '/usr/bin/vcfannotate1.0',
            'vcfannotategenotypes'       => '/usr/bin/vcfannotategenotypes1.0',
            'vcfbreakmulti'              => '/usr/bin/vcfbreakmulti1.0',
            'vcfcat'                     => '/usr/bin/vcfcat1.0',
            'vcfcheck'                   => '/usr/bin/vcfcheck1.0',
            'vcfclassify'                => '/usr/bin/vcfclassify1.0',
            'vcfcleancomplex'            => '/usr/bin/vcfcleancomplex1.0',
            'vcfcombine'                 => '/usr/bin/vcfcombine1.0',
            'vcfcommonsamples'           => '/usr/bin/vcfcommonsamples1.0',
            'vcfcountalleles'            => '/usr/bin/vcfcountalleles1.0',
            'vcfcreatemulti'             => '/usr/bin/vcfcreatemulti1.0',
            'vcfdistance'                => '/usr/bin/vcfdistance1.0',
            'vcfecho'                    => '/usr/bin/vcfecho1.0',
            'vcfentropy'                 => '/usr/bin/vcfentropy1.0',
            'vcfevenregions'             => '/usr/bin/vcfevenregions1.0',
            'vcffilter'                  => '/usr/bin/vcffilter1.0',
            'vcffixup'                   => '/usr/bin/vcffixup1.0',
            'vcfflatten'                 => '/usr/bin/vcfflatten1.0',
            'vcfgeno2alleles'            => '/usr/bin/vcfgeno2alleles1.0',
            'vcfgeno2haplo'              => '/usr/bin/vcfgeno2haplo1.0',
            'vcfgenosamplenames'         => '/usr/bin/vcfgenosamplenames1.0',
            'vcfgenosummarize'           => '/usr/bin/vcfgenosummarize1.0',
            'vcfgenotypecompare'         => '/usr/bin/vcfgenotypecompare1.0',
            'vcfgenotypes'               => '/usr/bin/vcfgenotypes1.0',
            'vcfglbound'                 => '/usr/bin/vcfglbound1.0',
            'vcfglxgt'                   => '/usr/bin/vcfglxgt1.0',
            'vcfhetcount'                => '/usr/bin/vcfhetcount1.0',
            'vcfhethomratio'             => '/usr/bin/vcfhethomratio1.0',
            'vcfindex'                   => '/usr/bin/vcfindex1.0',
            'vcfinfo2qual'               => '/usr/bin/vcfinfo2qual1.0',
            'vcfinfosummarize'           => '/usr/bin/vcfinfosummarize1.0',
            'vcfintersect'               => '/usr/bin/vcfintersect1.0',
            'vcfkeepgeno'                => '/usr/bin/vcfkeepgeno1.0',
            'vcfkeepinfo'                => '/usr/bin/vcfkeepinfo1.0',
            'vcfkeepsamples'             => '/usr/bin/vcfkeepsamples1.0',
            'vcfleftalign'               => '/usr/bin/vcfleftalign1.0',
            'vcflength'                  => '/usr/bin/vcflength1.0',
            'vcfnumalt'                  => '/usr/bin/vcfnumalt1.0',
            'vcfoverlay'                 => '/usr/bin/vcfoverlay1.0',
            'vcfparsealts'               => '/usr/bin/vcfparsealts1.0',
            'vcfprimers'                 => '/usr/bin/vcfprimers1.0',
            'vcfqual2info'               => '/usr/bin/vcfqual2info1.0',
            'vcfrandom'                  => '/usr/bin/vcfrandom1.0',
            'vcfrandomsample'            => '/usr/bin/vcfrandomsample1.0',
            'vcfremap'                   => '/usr/bin/vcfremap1.0',
            'vcfremoveaberrantgenotypes' => '/usr/bin/vcfremoveaberrantgenotypes',
            'vcfremovesamples'           => '/usr/bin/vcfremovesamples1.0',
            'vcfroc'                     => '/usr/bin/vcfroc1.0',
            'vcfsample2info'             => '/usr/bin/vcfsample2info1.0',
            'vcfsamplediff'              => '/usr/bin/vcfsamplediff1.0',
            'vcfsamplenames'             => '/usr/bin/vcfsamplenames1.0',
            'vcfsitesummarize'           => '/usr/bin/vcfsitesummarize1.0',
            'vcfstats'                   => '/usr/bin/vcfstats1.0',
            'vcfstreamsort'              => '/usr/bin/vcfstreamsort1.0',
            'vcfuniq'                    => '/usr/bin/vcfuniq1.0',
            'vcfuniqalleles'             => '/usr/bin/vcfuniqalleles1.0',

            # script tools
            bed2region                   => '/usr/share/vcflib-tgi1.0/bed2region',
            'plot_roc.r'                 => '/usr/share/vcflib-tgi1.0/plot_roc.r',
            'vcf2bed.py'                 => '/usr/share/vcflib-tgi1.0/vcf2bed.py',
            'vcf2sqlite.py'              => '/usr/share/vcflib-tgi1.0/vcf2sqlite.py',
            vcf_strip_extra_headers      => '/usr/share/vcflib-tgi1.0/vcf_strip_extra_headers',
            vcfbiallelic                 => '/usr/share/vcflib-tgi1.0/vcfbiallelic',
            vcfclearid                   => '/usr/share/vcflib-tgi1.0/vcfclearid',
            vcfclearinfo                 => '/usr/share/vcflib-tgi1.0/vcfclearinfo',
            vcfcomplex                   => '/usr/share/vcflib-tgi1.0/vcfcomplex',
            vcffirstheader               => '/usr/share/vcflib-tgi1.0/vcffirstheader',
            'vcfgtcompare.sh'            => '/usr/share/vcflib-tgi1.0/vcfgtcompare.sh',
            vcfindelproximity            => '/usr/share/vcflib-tgi1.0/vcfindelproximity',
            vcfindels                    => '/usr/share/vcflib-tgi1.0/vcfindels',
            vcfmultiallelic              => '/usr/share/vcflib-tgi1.0/vcfmultiallelic',
            vcfmultiway                  => '/usr/share/vcflib-tgi1.0/vcfmultiway',
            vcfmultiwayscripts           => '/usr/share/vcflib-tgi1.0/vcfmultiwayscripts',
            vcfnobiallelicsnps           => '/usr/share/vcflib-tgi1.0/vcfnobiallelicsnps',
            vcfnoindels                  => '/usr/share/vcflib-tgi1.0/vcfnoindels',
            vcfnosnps                    => '/usr/share/vcflib-tgi1.0/vcfnosnps',
            vcfnulldotslashdot           => '/usr/share/vcflib-tgi1.0/vcfnulldotslashdot',
            'vcfplotaltdiscrepancy.r'    => '/usr/share/vcflib-tgi1.0/vcfplotaltdiscrepancy.r',
            'vcfplotaltdiscrepancy.sh'   => '/usr/share/vcflib-tgi1.0/vcfplotaltdiscrepancy.sh',
            'vcfplotsitediscrepancy.r'   => '/usr/share/vcflib-tgi1.0/vcfplotsitediscrepancy.r',
            'vcfplottstv.sh'             => '/usr/share/vcflib-tgi1.0/vcfplottstv.sh',
            'vcfprintaltdiscrepancy.r'   => '/usr/share/vcflib-tgi1.0/vcfprintaltdiscrepancy.r',
            'vcfprintaltdiscrepancy.sh'  => '/usr/share/vcflib-tgi1.0/vcfprintaltdiscrepancy.sh',
            vcfqualfilter                => '/usr/share/vcflib-tgi1.0/vcfqualfilter',
            vcfregionreduce              => '/usr/share/vcflib-tgi1.0/vcfregionreduce',
            vcfregionreduce_and_cut      => '/usr/share/vcflib-tgi1.0/vcfregionreduce_and_cut',
            vcfregionreduce_pipe         => '/usr/share/vcflib-tgi1.0/vcfregionreduce_pipe',
            vcfregionreduce_uncompressed => '/usr/share/vcflib-tgi1.0/vcfregionreduce_uncompressed',
            vcfremovenonATGC             => '/usr/share/vcflib-tgi1.0/vcfremovenonATGC',
            vcfsnps                      => '/usr/share/vcflib-tgi1.0/vcfsnps',
            vcfsort                      => '/usr/share/vcflib-tgi1.0/vcfsort',
            vcfvarstats                  => '/usr/share/vcflib-tgi1.0/vcfvarstats',
        }
    );

    if (exists $toolset{$version}) {
        return $toolset{$version}
    }

    die $self->error_message(
        "Couldn't find the vcflib tool set for version: '%s'",
        $version
    );
}

1;
