package Genome::Model::Tools::Sample::NameConversion;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Sample::NameConversion {
    is  => 'Command',
    has => [
        input => {
            type => 'String',
            doc  => "provide a file of sample names",
        },
        output => {
            type => 'String',
            doc  => "provide a file name to write the converted names to.",
            is_optional => 1,
        },
        short_to_long => {
            is => 'Boolean',
            doc => "a conversion method that should be used to convert your old TCGA in house short name to the new in house long name",
            is_optional => 1,
        },
    ],
    doc => "A tool to convert sample names."
};

sub help_synopsis {
    return <<EOS

gmt sample name-coversion -h

EOS
}

sub help_detail {
    return <<EOS

	This Tool will convert a sample name.

	short-to-long will run the query
	      sqlrun "select * from genomic_dna where genomic_name like 'sample'"
	      and parse the new sample name from the results 
	      your sample name will be returned with _not_found attached if the query fails
EOS
}

sub execute {

    my $self = shift;

    unless ( $self->short_to_long ) {
        $self->error_message( "\n\nPlease indicate the conversion method."
              . " See gmt sample name-conversion --help for options.\n\n" );
        return;
    }

    my $output = $self->output;
    if ( $self->output ) {
        open( OUT, ">$output" )
          || $self->error_message(
            "\n\ncouldn't open the output file $output\n\n")
          && return;
    }

    my $input = $self->input;

    #unless (-f $input) { print "\n\nCould not open the input $input.\n\n";return; }
    open( IN, $input )
      || $self->error_message("\n\ncouldn't open the input file $input\n\n")
      && return;
    while (<IN>) {
        chomp;
        my $name = $_;

        my $converted_name;
        if ( $self->short_to_long ) {
            ($converted_name) = &sql_query($name);
        }
        if ( $self->output ) {
            print OUT qq($converted_name\n);
        }
        else {
            print qq($converted_name\n);
        }
    }
    close IN;
    close OUT;
}

sub sql_query {

    my ($name) = @_;

    my ( $prefix, $id ) = $name =~ /([\S]+)\-([\S]+)/;
    unless ( $prefix && $id ) {
        my $return_name = "$name\_id_not_found";
        return ($return_name);
    }
    my ($type) = $id =~ /(\D)$/;
    unless ($type) {
        my $return_name = "$name\_type_not_found";
        return ($return_name);
    }

    my $match;
    if ( $type eq "t" ) {
        $match = "0" . "\\d";
    }
    elsif ( $type eq "n" ) {
        $match = "[1-9]" . "\\d";
    }
    else {
        my $return_name = "$name\_match_not_found";
        return ($return_name);
    }

    $id =~ s/\D$//;
    my $like = $prefix . "%" . $id . "%";

    my @run =
      `sqlrun "select * from genomic_dna where genomic_name like '$like'"`;

    my ($return_name);
    for my $line (@run) {
        if ( $line =~ /(\S+$id\-$match\S+)/ ) {
            $return_name = $1;

            #print qq(\t$return_name\n);
        }
        if ( $line =~ /$name/ ) {
            unless ($return_name) { $return_name = $name; }
        }
    }
    unless ($return_name) {
        $return_name = "$name\_query_not_found";
        return ($return_name);
    }
    return ($return_name);
}

1;
