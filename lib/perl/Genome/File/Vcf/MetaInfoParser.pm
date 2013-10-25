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
      { $return = {$item{flag} => "IS_VCF_FLAG"} }

    flag: string

    key: /[\w\d]+/

    value: map
      | string
      | list 

    list: string "," list
      { unshift @{$item{list}}, $item{string}; $return = $item{list} }
      | string
      { $return = [$item{string}] }

    arbitraryString:  /.+/

    string: /[\d\w\s-_\/\:\.]+/
      | <perl_quotelike>
      { $return = $item[1]->[2] }

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

1;

