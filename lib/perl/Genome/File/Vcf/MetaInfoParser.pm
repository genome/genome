package Genome::File::Vcf::MetaInfoParser;

use strict;
use warnings;
use Genome;

use Parse::RecDescent;

#$::RD_TRACE = 1;
my $grammar = q{
    evaluate_metainfo: map
      | arbitraryString 

    map: "<" pairs ">"
      { $return = $item{pairs} }

    pairs: pair "," pairs
      { $return = {%{$item{pairs}}, %{$item{pair}}}}
      | pair

    pair: key "=" value
      { $return = {$item{key} => $item{value}} }
      | flag
      { $return = {$item{flag} => undef} }

    flag: key

    key: /[\w\d]+/

    value: map
      | string
      | list 

    list: string "," list
      { unshift @{$item{list}}, $item{string}; $return = $item{list} }
      | string
      { $return = [$item{string}] }

    arbitraryString:  /.+/
      { $return = Genome::File::Vcf::Header::String->new(content => $item[1], is_quoted => 0)}

    string: /[^<>=,"]+/
      { $return = Genome::File::Vcf::Header::String->new(content => $item[1], is_quoted => 0)}
      | <perl_quotelike>
      { $return = Genome::File::Vcf::Header::String->new(content => $item[1]->[2],
                                                            is_quoted => 1)}

};

my $parser;
sub parse {
    my $self = shift;
    my $str = shift;

    $parser = Parse::RecDescent->new($grammar) unless $parser;
    unless ($parser) {
        $self->error_message("Failed to create parser");
        return;
    }
    return $parser->evaluate_metainfo($str);
}

package Genome::File::Vcf::Header::String;
sub new {
    my ($class, %args) = @_;
    return bless { %args }, $class;
}


1;

