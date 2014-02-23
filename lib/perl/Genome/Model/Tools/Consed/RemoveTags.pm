package Genome::Model::Tools::Consed::RemoveTags;

use strict;
use warnings;

use Genome;
use Data::Dumper 'Dumper';

class Genome::Model::Tools::Consed::RemoveTags {
    is => 'Command',
    has => [
        ace_in => {
            is => 'Text',
            doc => 'Input ace file',
        },
        ace_out => {
            is => 'Text',
            doc => 'Output ace file',
            is_optional => 1,
            is_mutable => 1,
        },
        tags => {
            is => 'Text',
            doc => 'List of tag types to remove',
            is_many => 1,
        },
    ],
};

sub help_brief {
    'Tool to remove tags from ace file',
}

sub help_synopsis {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    #check ace in
    if ( not -s $self->ace_in ) {
        $self->error_message("Failed to find input ace or file is zero size: ".$self->ace_in);
        return;
    }

    #set ace out
    if ( not $self->ace_out ) {
        $self->ace_out( $self->ace_in.'.tags_removed' );
    } 

    #ace reader
    my $reader = Genome::Model::Tools::Consed::AceReader->create(
        file => $self->ace_in,
    );
    if ( not $reader ) {
        $self->error_message("Failed to create ace reader for file: ".$self->ace_in);
        return;
    }

    #ace writer
    my $writer = Genome::Model::Tools::Consed::AceWriter->create(
        file => $self->ace_out,
    );
    if ( not $writer ) {
        $self->error_message("Failed to create ace writer for file: ".$self->ace_out);
        return;
    }

    #write contigs to new ace file
    while ( my $contig = $reader->next_contig ) {
        $self->debug_message("Exporting contig: ".$contig->{name});
        $writer->add_contig( contig => $contig );
    }

    #write tags
    my @tags;
    for my $tag ( @{$reader->contig_tags} ) {
        if ( grep { $tag->{tag_type} eq lc $_ } $self->tags ) {
            $self->debug_message("Excluding tag with type: ".$tag->{tag_type});
            next;
        }
        push @tags, $tag;
    }

    $writer->add_contig_tags( \@tags ) if @tags;
    $writer->add_assembly_tags( $reader->assembly_tags ) if $reader->assembly_tags;

    $writer->close;

    return 1;
}

1;
