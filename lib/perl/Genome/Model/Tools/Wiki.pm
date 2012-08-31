package Genome::Model::Tools::Wiki;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Wiki {
    is => 'Command',
    has => [
        _wiki => { is => 'GSCApp::Wikibot', is_optional => 1 }
    ],
};

sub post_page_to_wiki {
    my $self = shift;
    
    unless($self->_wiki) {
        $self->_wiki(GSCApp::Wikibot->new());    
    }
    
    return $self->_wiki->post_page_to_wiki(@_);
}

sub create_section {
    my $self = shift;
    my $heading = shift;
    my $body = shift;
    
    #Use level one headings by default; wrap $heading in more '='s for lower-level headings
    return "\n=" . $heading . "=\n" . $body . "\n";
}

sub create_item {
    my $self = shift;
    my $name = shift;
    my @values = @_;
    
    return ';' . $name . join("\n", map(': ' . $_, @values) ) . "\n";
}

sub help_brief {
    "Tools to work with the wiki.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt wiki ...    
EOS
}

sub help_detail {                           
    return <<EOS 
Tools to work with the wiki.
EOS
}

1;
