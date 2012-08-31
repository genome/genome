package Genome::Model::Tools::Assembly::RemoveReads;

use strict;
use warnings;

use Genome;
use Cwd;
use IO::File;

class Genome::Model::Tools::Assembly::RemoveReads
{
    is => 'Command',
    has => 
    [
        ace_file => {
            type => "String",
            optional => 0,
            doc => "This is the input ace file"
        }, 
        contig => {
            type => "String",
            optional => 0,
            doc => "This is the contig that is going to have reads removed",
        },
	    read_list => {
            type => "String",
            optional => 0,
		    doc => "This is the comma delimited list of reads to remove from the contig",
	    },
	    out_file_name => {
            type => "String",
            optional => 0,
		    doc => "This is the name of the output file",
	    },	    
    ]
};

sub help_brief {
    ""
}

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    remove-reads --ace-file=in.ace --contig=Contig0.1 --read-list=read1.g1,read2.b1 --out-file-name=out.ace
EOS
}

sub execute
{
    my $self = shift;
    my $ace_file = $self->ace_file;
    my $contig_name = $self->contig;
    my $read_list = $self->read_list;
    my $out_file_name = $self->out_file_name;
    
    
    my $ct = Genome::Model::Tools::Pcap::ContigTools->new;
    my $phd_object;
    $self->error_message("Ace file does not exist\n") and return unless (-e $ace_file);
    my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file, using_db => 1);
    
    if(-e "../phdball_dir/phd.ball.1")
    {
        $phd_object = Genome::Model::Tools::Pcap::Phd->new(input_file => "../phdball_dir/phd.ball.1");
    }
    elsif(-e "../phd_dir/")
    {
        $phd_object = Genome::Model::Tools::Pcap::Phd->new(input_directory => "../phd_dir/");
    }
    else
    {
        $self->error_message("Need to either have a ../phd_dir or a phdball file named ../phdball_dir/phd.ball.1") and return;
    }   
     
    my %params;


	my $clip_contig = $ao->get_contig($contig_name);
    
	$params{remove_reads} = [ split /,/, $read_list ];
    
	if($ct->remove_reads($clip_contig,$phd_object,%params))
	{
		$ao->add_contig($clip_contig);
	}
	else
	{
		$ao->remove_contig($clip_contig->name);
	}
    
    $ao->write_file(output_file => $out_file_name);
    
    return 1;	
}

1;



