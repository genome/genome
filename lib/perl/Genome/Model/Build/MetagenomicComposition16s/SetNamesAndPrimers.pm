package Genome::Model::Build::MetagenomicComposition16s::SetNamesAndPrimers;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::MetagenomicComposition16s::SetNamesAndPrimers {};

sub set_names_and_primers_for {
    my ($self, $primer_set_name) = @_;
    Carp::confess('No primer set name to get amplicon set names and primers!') if not $primer_set_name;
    my $metohod = '_'.$primer_set_name;
    return $self->$metohod;
}

sub _sanger {
    return ( '' => [], );
}

sub _454 {
    return (
        V1_V3 => [qw/
            ATTACCGCGGCTGCTGG 
        /],
        V3_V5 => [qw/ 
            CCGTCAATTCATTTAAGT
            CCGTCAATTCATTTGAGT
            CCGTCAATTCCTTTAAGT
            CCGTCAATTCCTTTGAGT
        /],
        V6_V9 => [qw/
            TACGGCTACCTTGTTACGACTT
            TACGGCTACCTTGTTATGACTT
            TACGGTTACCTTGTTACGACTT
            TACGGTTACCTTGTTATGACTT
        /],
    );
}

sub X_454 { # New primers for 454. Not sure if we will use these regularly.
    return (
        V1_V3 => [qw/
            ATTACCGCGGCTGCTGG
            ATTTACCGCGGCTGCTGG
            ATTTTACCGCGGCTGCTGG
        /],
        V3_V5 => [qw/
            CCGTCAATTCATTTAAGT
            CCGTCAATTTCATTTAAGT
            CCGTCAATTCATTTTAAGT
            CCGTCAATTTCATTTTAAGT
            CCGTCAATTCATTTGAGT
            CCGTCAATTTCATTTGAGT
            CCGTCAATTCATTTTGAGT
            CCGTCAATTTCATTTTGAGT
            CCGTCAATTCCTTTAAGT
            CCGTCAATTTCCTTTAAGT
            CCGTCAATTCCTTTTAAGT
            CCGTCAATTTCCTTTTAAGT
            CCGTCAATTCCTTTGAGT
            CCGTCAATTTCCTTTGAGT
            CCGTCAATTCCTTTTGAGT
            CCGTCAATTTCCTTTTGAGT
        /],
    );
}

sub _solexa { #not yet implemented
    return ( '' => [], ); # FIXME rm this when implemented
    return (
        #GA[AG]TTTGATC[ACT]TGGCTCAG
        'V1.F' => [qw/
        GAATTTGATCATGGCTCAG
        GAATTTGATCCTGGCTCAG
        GAATTTGATCTTGGCTCAG
        GAGTTTGATCATGGCTCAG
        GAGTTTGATCCTGGCTCAG
        GAGTTTGATCTTGGCTCAG
        /],
        #TTACT[AC]ACCC[GT]T[CGT]CGCC
        'V1.R' => [qw/
        TTACTAACCCGTCCGCC
        TTACTCACCCGTCCGCC
            TTACTAACCCTTCCGCC
            TTACTCACCCTTCCGCC
            TTACTAACCCGTGCGCC
            TTACTCACCCGTGCGCC
            TTACTAACCCTTGCGCC
            TTACTCACCCTTGCGCC
            TTACTAACCCGTTCGCC
            TTACTCACCCGTTCGCC
            TTACTAACCCTTTCGCC
            TTACTCACCCTTTCGCC
        /],
        #GGCG[GAC]A[CA]GGGT[TG]AGTAA
        'V2.F' => [qw/
            GGCGGACGGGTTAGTAA
            GGCGAACGGGTTAGTAA
            GGCGCACGGGTTAGTAA
            GGCGGAAGGGTTAGTAA
            GGCGAAAGGGTTAGTAA
            GGCGCAAGGGTTAGTAA
            GGCGGACGGGTGAGTAA
            GGCGAACGGGTGAGTAA
            GGCGCACGGGTGAGTAA
            GGCGGAAGGGTGAGTAA
            GGCGAAAGGGTGAGTAA
            GGCGCAAGGGTGAGTAA
        /],
        #CC[AGT]TTACC[CT]CACC[AT]ACTA
        'V2.R' => [qw/
            CCATTACCCCACCAACTA
            CCGTTACCCCACCAACTA
            CCTTTACCCCACCAACTA
            CCATTACCTCACCAACTA
            CCGTTACCTCACCAACTA
            CCTTTACCTCACCAACTA
            CCATTACCCCACCTACTA
            CCGTTACCCCACCTACTA
            CCTTTACCCCACCTACTA
            CCATTACCTCACCTACTA
            CCGTTACCTCACCTACTA
            CCTTTACCTCACCTACTA
        /],
        #CCTACGG[AG][AGT]GGC[ACT]GCAG
        'V3.F' => [qw/
            CCTACGGAAGGCAGCAG
            CCTACGGGAGGCAGCAG
            CCTACGGAGGGCAGCAG
            CCTACGGGGGGCAGCAG
            CCTACGGATGGCAGCAG
            CCTACGGGTGGCAGCAG
            CCTACGGAAGGCCGCAG
            CCTACGGGAGGCCGCAG
            CCTACGGAGGGCCGCAG
            CCTACGGGGGGCCGCAG
            CCTACGGATGGCCGCAG
            CCTACGGGTGGCCGCAG
            CCTACGGAAGGCTGCAG
            CCTACGGGAGGCTGCAG
            CCTACGGAGGGCTGCAG
            CCTACGGGGGGCTGCAG
            CCTACGGATGGCTGCAG
            CCTACGGGTGGCTGCAG
        /],
        #T[GT]ACCGC[AG]GCTGCTGGCAC
        'V3.R' => [qw/
            TGACCGCAGCTGCTGGCAC
            TTACCGCAGCTGCTGGCAC
            TGACCGCGGCTGCTGGCAC
            TTACCGCGGCTGCTGGCAC
        /],
        #GTGCCAGCAGC[CT]GCGGT[AC]A
        'V4.F' => [qw/
            GTGCCAGCAGCCGCGGTAA
            GTGCCAGCAGCTGCGGTAA
            GTGCCAGCAGCCGCGGTTA
            GTGCCAGCAGCTGCGGTTA
        /],
        #CG[GC]ATTTCAC[CT][GC]CTAC
        'V4.R' => [qw/
            CGGATTTCACCGCTAC
            CGCATTTCACCGCTAC
            CGGATTTCACTGCTAC
            CGCATTTCACTGCTAC
            CGGATTTCACCCCTAC
            CGCATTTCACCCCTAC
            CGGATTTCACTCCTAC
            CGCATTTCACTCCTAC
        /],
        #C[AG]AAC[ACG]GGATTAGATACCC
        'V5.F' => [qw/
            CAAACAGGATTAGATACCC
            CGAACAGGATTAGATACCC
            CAAACCGGATTAGATACCC
            CGAACCGGATTAGATACCC
            CAAACGGGATTAGATACCC
            CGAACGGGATTAGATACCC
        /],
        #CCCGTCAATT[CT][ACT]TTT[AG]AGT
        'V5.R' => [qw/
            CCCGTCAATTCATTTAAGT
            CCCGTCAATTTATTTAAGT
            CCCGTCAATTCCTTTAAGT
            CCCGTCAATTTCTTTAAGT
            CCCGTCAATTCTTTTAAGT
            CCCGTCAATTTTTTTAAGT
            CCCGTCAATTCATTTGAGT
            CCCGTCAATTTATTTGAGT
            CCCGTCAATTCCTTTGAGT
            CCCGTCAATTTCTTTGAGT
            CCCGTCAATTCTTTTGAGT
            CCCGTCAATTTTTTTGAGT
        /],
        #A[CA][GA]CGA[AG]GAACCTTACC
        'V6.F' => [qw/
            ACGCGAAGAACCTTACC
            AAGCGAAGAACCTTACC
            ACACGAAGAACCTTACC
            AAACGAAGAACCTTACC
            ACGCGAGGAACCTTACC
            AAGCGAGGAACCTTACC
            ACACGAGGAACCTTACC
            AAACGAGGAACCTTACC
        /],
        #AC[AG][AG]CACGAGCTG[AT]CGAC
        'V6.R' => [qw/
            ACAACACGAGCTGACGAC
            ACGACACGAGCTGACGAC
            ACAGCACGAGCTGACGAC
            ACGGCACGAGCTGACGAC
            ACAACACGAGCTGGCGAC
            ACGACACGAGCTGGCGAC
            ACAGCACGAGCTGGCGAC
            ACGGCACGAGCTGGCGAC
        /],
        #[GT][CT]AACGAGCGCAACCCTT
        'V7.F' => [qw/
            GCAACGAGCGCAACCCTT
            TCAACGAGCGCAACCCTT
            GTAACGAGCGCAACCCTT
            TTAACGAGCGCAACCCTT
        /],
        #CGTC[AG]TCC[CT]C[AT]CCTTCC
        'V7.R' => [qw/
            CGTCATCCCCACCTTCC
            CGTCGTCCCCACCTTCC
            CGTCATCCTCACCTTCC
            CGTCGTCCTCACCTTCC
            CGTCATCCCCTCCTTCC
            CGTCGTCCCCTCCTTCC
            CGTCATCCTCTCCTTCC
            CGTCGTCCTCTCCTTCC
        /],
        #CTA[Cg]A[Ca]ACG[Tc]GCTACAATG
        'V8.F' => [qw/
            CTACACACGTGCTACAATG
            CTAGACACGTGCTACAATG
            CTACAAACGTGCTACAATG
            CTAGAAACGTGCTACAATG
            CTACACACGCGCTACAATG
            CTAGACACGCGCTACAATG
            CTACAAACGCGCTACAATG
            CTAGAAACGCGCTACAATG
        /],
        #CCG[AG]GAACGTATTCAC[GC]
        'V8.R' => [qw/
            CCGAGAACGTATTCACG
            CCGGGAACGTATTCACG
            CCGAGAACGTATTCACC
            CCGGGAACGTATTCACC
        /],
        #CGTTC[CT]CGGG[CT]CTTGTAC
        #[GC]GTGAATACGTTC[CT]CGG
        'V9.F' => [qw/
            CGTTCCCGGGCCTTGTAC
            CGTTCTCGGGCCTTGTAC
            CGTTCCCGGGTCTTGTAC
            CGTTCTCGGGTCTTGTAC

            GGTGAATACGTTCCCGG
            CGTGAATACGTTCCCGG
            GGTGAATACGTTCTCGG
            CGTGAATACGTTCTCGG
        /],
        #CGG[CT]TACCTTGTTA[CT]GACTT
        'V9.R' => [qw/
            CGGCTACCTTGTTACGACTT
            CGGTTACCTTGTTACGACTT
            CGGCTACCTTGTTATGACTT
            CGGTTACCTTGTTATGACTT
        /],
    );
}

1;

