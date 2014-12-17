
package Genome::File::Vcf::BamReadcountParser;

use strict;
use warnings;

sub encode {
    my $line = shift;

    $line =~ s/([\t])/?/g;
    $line =~ s/([:])/;/g;
    return $line;
}

sub decode {
    my $line = shift;

    $line =~ s/([?])/\t/g;
    $line =~ s/([;])/:/g;
    return $line;
}
1;

