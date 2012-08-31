package Genome::Model::Tools::ImportAnnotation::Venter;

use strict;
use warnings;

use Genome;
use IO::File;
use Text::CSV_XS;
use File::Slurp;
use Carp;

class Genome::Model::Tools::ImportAnnotation::Venter {
    is => 'Command',
    has => [
#        file => { is => "Array",
#                  doc => "",
#                },
        outputdir => { is  => "Text",
                       doc => "output directory",
                     },
    ],
    has_many => [
        files =>{ is  => "Text",
                  doc => "files",
        },

    ],
};


sub sub_command_sort_position {12}

sub help_brief
{
    "Import Venter snps to the file based data source";
}

sub help_synopsis
{
    return <<EOS
change this
gmt import-annotation venter --files <file1,file2> \
 --outputdir <output directory>
EOS
}

sub help_detail
{
    return <<EOS
This tool is used for importing the Venter snp data files to a file based datasource.
EOS
}


sub execute
{
    my $self = shift;
    #$DB::single = 1;
    $self->process_venter_snps();
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


sub process_venter_snps
{
    my $self = shift;
    my @files = $self->files;
    my $var_id = 1;
    my $outputdir = $self->outputdir;
    my $c = Text::CSV_XS->new({sep_char => "\t"});

    foreach my $f (@files)
    {
        my $fh = IO::File->new($f);

        while(my $line = <$fh>)
        {
            $c->parse($line);
            my @f = $c->fields();
            if($#f == 8)
            {
                my $chromosome_name = $f[0];
                my $extvar_id = $f[1];
                my ($start,$stop) = ($f[3],$f[4]);
                
                my $meta = $f[7];
                my ($alleles,$rmr,$tr) = split(/;/,$meta);
                my @variations = ($var_id,$extvar_id,$alleles,"indel",$chromosome_name,
                                  $start,$stop,undef);
                $c->combine(@variations);
                write_file($outputdir."/variations.csv",
                           {append => 1},
                           $c->string()."\n");
            
                my @var_inst = ($var_id,"venter",2,undef);
                $c->combine(@var_inst);
                write_file($outputdir."/variation_instances.csv",
                           {append => 1},
                           $c->string()."\n");
            }
            elsif($#f == 10)
            {
                my $chromosome_name = $f[0];
                my $extvar_id = $f[1];
                my ($start,$stop) = ($f[3],$f[4]);

                my $alleles = $f[9];
                if($f[10] =~ /Deletion/)
                {
                    $alleles = $alleles."/-";
                }
                elsif($f[10] =~ /Insertion/)
                {
                    $alleles = "-/".$alleles;
                }
                my @variations = ($var_id,$extvar_id,$alleles,"indel",$chromosome_name,
                                  $start,$stop,undef);
                $c->combine(@variations);
                write_file($outputdir."/variations.csv",
                           {append => 1},
                           $c->string()."\n");
            
                my @var_inst = ($var_id,"venter",2,undef);
                $c->combine(@var_inst);
                write_file($outputdir."/variation_instances.csv",
                           {append => 1},
                           $c->string()."\n");
            
            }
            $var_id += 1;
        }
    }
    return 1;
}

1;

# $Id$
