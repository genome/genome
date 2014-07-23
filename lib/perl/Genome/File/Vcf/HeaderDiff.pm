package Genome::File::Vcf::HeaderDiff;

use strict;
use warnings FATAL => 'all';
use Params::Validate qw(validate validate_pos :types);

=head1 NAME

Genome::File::Vcf::HeaderDiff - Representation of differences between
two VCF file headers.

=head1 SYNOPSIS

    use Genome::File::Vcf::Differ;

    my $differ = Genome::File::Vcf::Differ->new('filea.vcf', 'fileb.vcf');
    my $diff = $differ->header;
    if ($diff) {
       $diff->print;
    }

=head1 DESCRIPTION

    The HeaderDiff print method displays the differences between two VCF file headers.

=cut

sub new {
    my ($class, $path_a, $diffs_a, $path_b, $diffs_b) = validate_pos(@_,
        {type => SCALAR},
        {type => SCALAR},
        {type => ARRAYREF},
        {type => SCALAR},
        {type => ARRAYREF},
    );
    my $self = {
        _a => $path_a,
        _b => $path_b,
        _diffs_a => $diffs_a,
        _diffs_b => $diffs_b,
    };

    bless $self, $class;
    return $self;
}

sub print {
    my $self = shift;

    print $self->to_string . "\n";
}

sub to_string {
    my $self = shift;

    return _to_string($self->{_a}, @{$self->{_diffs_a}}) .
           _to_string($self->{_b}, @{$self->{_diffs_b}});
}

sub _to_string {
    my $file_name = shift;
    my @lines = @_;

    my $indent = '    ';
    printf "Lines unique to %s are:\n%s%s\n", $file_name, $indent, join("\n$indent", @lines);
}

1;
