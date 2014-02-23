package Genome::Model::Tools::ImportAnnotation::DbSnp;

use strict;
use warnings;

use Genome;
use IO::File;
use Text::CSV_XS;
use File::Slurp;
use IPC::Run;
use Carp;

class Genome::Model::Tools::ImportAnnotation::DbSnp {
    is => 'Command',
    has => [
#        file => { is => "Array",
#                  doc => "input files",
#                },
        outputdir => { is => "Text",
                       doc => "output directory for data sources",
                     },
        start_id => { is => "Number",
                      doc => "variation id to start at",
                      default => 1,
                    },
        ],
    has_many => [
        files => { is => "Text",
                   doc => "input files", },
                ],
    
};


sub sub_command_sort_position {12}

sub help_brief
{
    "Import dbsnp snps to the file based data source";
}

sub help_synopsis
{
    return <<EOS
change this
gmt import-annotation db-snp.....
EOS
}

sub help_detail
{
    return <<EOS
fix me.
EOS
}


sub execute
{
    my $self = shift;
    $self->process_dbsnp();
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


sub process_dbsnp
{
    my $self   = shift;
    my @files  = $self->files;
    my $var_id = $self->start_id;
    my $outdir = $self->outputdir;
    my $c      = Text::CSV_XS->new( { sep_char => "\t" } );

    foreach my $f (@files)
    {
        my $fh = IO::File->new($f);

        { # yay! we need to eff with $/ !
            local $/ = "\n\n";
            while(my $line = <$fh>)
            {
                my $rec = $self->parse_record($line);
                my @variations = ($var_id, $rec->{extvarid}, $rec->{alleles},
                                  $rec->{class}, $rec->{chromosome} , 
                                  $rec->{start}, $rec->{stop} );
                $c->combine(@variations);
                write_file( $outdir."/variations.csv",
                            {append => 1},
                            $c->string()."\n");

                #variation instance stuff
                my @var_inst = ($var_id,"dbSNP",2,undef);
                $c->combine(@var_inst);
                write_file( $outdir."/variation_instance.csv",
                            {append => 1},
                            $c->string()."\n");
                $var_id += 1; 
            }

        } # end of the block where we mess with $/
        $fh->close;
    }
    $self->debug_message("variation id at $var_id\n");
    return 1;
}

sub parse_record
{
    my $self = shift;
    my $rec = shift;

    my @lines = split(/\n/,$rec);

    my $rsid    = undef;
    my $alleles = undef;
    my $chrom   = undef;
    my $start   = undef;
    my $stop    = undef;
    my $class   = 'SNP';
    foreach my $line (@lines)
    {
        if($line =~ /^rs\d+/)
        {
            my @fields = split(/ \| /,$line);
            $rsid = $fields[0];
        }
        elsif($line =~ /^SNP/)
        {
            my @fields = split(/ \| /,$line);
            $alleles = $fields[1];
            $alleles =~ s/alleles='//;
            $alleles =~ s/'$//;
            if($alleles =~ /-/)
            {
                $class = 'indel';
            }
        }
        elsif($line =~ /^VAL/) # do we need to worry about validated snps?
        {
            next;
        }
        elsif($line =~ /^CTG | assembly=reference/)
        {
            my @fields = split(/ \| /,$line);
            $chrom = $fields[2];
            $start = $fields[3];
            $stop = $fields[3];
            $chrom =~ s/chr=(\S+)/$1/;
            $start =~ s/chr-pos=(\S+)/$1/;
            $stop =~ s/chr-pos=(\S+)/$1/;
        }
    }
    my %struct;
    $struct{chromosome} = $chrom;
    $struct{start} = $start;
    $struct{stop} = $stop;
    $struct{alleles} = $alleles;
    $struct{extvarid} = $rsid;
    $struct{class} = $class;
    return \%struct;
}


1;

# $Id$
