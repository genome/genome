package Genome::Model::Tools::ImportAnnotation::Watson;

use strict;
use warnings;

use Genome;
use Carp;
use Text::CSV_XS;
use IO::File;
use File::Slurp;
use IPC::Run;

class Genome::Model::Tools::ImportAnnotation::Watson {
    is  => 'Command',
    has => [
        outputdir => {
            is  => "Text",
            doc => "output directory",
        },
    ],
    has_many => [
        files => {
            is  => "Text",
            doc => "input files",
        },
    ],
};

sub sub_command_sort_position {12}

sub help_brief
{
    "Import watson snp/variation data";
}

sub help_synopsis
{
    return <<EOS

gmt import-annotation watson --outputdir <directory to dump annotation data> --file <input file>
EOS
}

sub help_detail
{
    return <<EOS
Fill this out.
EOS
}

sub execute
{
    my $self = shift;
    $self->process_watson_snp_file();
    my ($sort_stdout,$sort_stderr);
    IPC::Run::run(['sort',
                   '-n',
                   '-k5,6',
                   '-o','variations.csv.sorted',
                   'variations.csv'],
                  \undef,
                  '2>',
                  \$sort_stderr,
                  '>',
                  \$sort_stdout,

        ) or croak "problem sorting:\n$sort_stderr";


    return 1;
}

sub determine_file_type
{
    my $self = shift;
    my $c    = Text::CSV_XS->new( { sep_char => "\t" } );
    my $file = $self->file;
    my $fh   = IO::File->new($file);
    my $line = $fh->readline();
    $c->parse($line);
    my @f = $c->fields();

    return $#f;
}

sub process_watson_snp_file
{
    my $self = shift;

    my $c         = Text::CSV_XS->new( { sep_char => "\t" } );
    my $var_id    = 1;
    my $outputdir = $self->outputdir;
    #$DB::single = 1;
    my @files = $self->files;
    foreach my $f (@files)
    {
        my $fh = IO::File->new($f);
        while ( my $line = <$fh> )
        {
            $c->parse($line);
            my @f = $c->fields();

            my $chromosome_name = $f[0];
            $chromosome_name =~ s/chr//;
            my ( $start, $stop ) = ( $f[3], $f[4] );
            my $meta = $f[8];
            my ( $extvar_id, $alleles ) = split( /\s{0,1};\s{0,1}/, $meta );
            $alleles   =~ s/alleles //;
            $extvar_id =~ s/SNP //;
            my @variation = (
                $var_id, $extvar_id, $alleles, 'SNP', $chromosome_name,
                $start, $stop, undef
            );
            $c->combine(@variation);
            write_file(
                $outputdir . '/variations.csv',
                { append => 1 },
                $c->string() . "\n"
            );
            my @var_inst = ( $var_id, 'watson', 2, undef );
            $c->combine(@var_inst);
            write_file(
                $outputdir . '/variation_instances.csv',
                { append => 1 },
                $c->string() . "\n"
            );
            $var_id += 1;
        }
        $fh->close;
    }
    return 1;
}

1;
