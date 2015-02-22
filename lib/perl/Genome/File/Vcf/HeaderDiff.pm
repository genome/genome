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
    my ($class, $path_a, $line_diffs_a, $sample_diffs_a,
                $path_b, $line_diffs_b, $sample_diffs_b) = validate_pos(@_,
        {type => SCALAR},
        {type => SCALAR},
        {type => OBJECT},
        {type => OBJECT},
        {type => SCALAR},
        {type => OBJECT},
        {type => OBJECT},
    );
    my $self = {
        _a => $path_a,
        _b => $path_b,
        _line_diffs_a => [sort $line_diffs_a->members()],
        _line_diffs_b => [sort $line_diffs_b->members()],
        _sample_diffs_a => [sort $sample_diffs_a->members()],
        _sample_diffs_b => [sort $sample_diffs_b->members()],
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

    return _samples_to_string($self->{_a}, $self->{_sample_diffs_a}) .
           _samples_to_string($self->{_b}, $self->{_sample_diffs_b}) .
           _lines_to_string($self->{_a}, $self->{_line_diffs_a}) .
           _lines_to_string($self->{_b}, $self->{_line_diffs_b});
}

sub _lines_to_string {
    my ($file_name, $lines) = @_;

    if (@$lines) {
        my $indent = '    ';
        return sprintf "Lines unique to %s are:\n%s%s\n", $file_name, $indent,
            join("\n$indent", @$lines);
    } else {
        return '';
    }
}

sub _samples_to_string {
    my ($file_name, $names) = @_;

    if (@$names) {
        my $indent = '    ';
        return sprintf "Sample names unique to %s are:\n%s%s\n", $file_name,
            $indent, join("\n$indent", @$names);
    } else {
        return '';
    }

}

1;
