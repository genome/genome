package Genome::Model::Tools::Assembly::SplitContig;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::Tools::Assembly::SplitContig
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
            doc => "This is the contig that is going to be split",
        },
	    split_position => {
            type => "String",
            optional => 0,
		    doc => "This is the split position, in padded units",
	    },
	    no_gui => {
            type => "String",
            optional => 1,
		    doc => "This allows the tool to run without a GUI, using the defaults to determine which contigs a read is moved to",
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
    split-contig --ace-file=in.ace --contig=Contig0.1 --split-position=500 --out-file-name=out.ace
EOS
}

sub execute
{
    my $self = shift;
    my $ace_file = $self->ace_file;
    my $split_contig_name = $self->contig;
    my $split_position = $self->split_position;
    my $out_file_name = $self->out_file_name;
    my $no_gui = $self->no_gui;

    $self->error_mesage("Input ace file $ace_file does not exist\n") and return unless (-e $ace_file);
    
    my $ace_object = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file, using_db => 1);
    
    my $phd_object;    
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
    
    my $contig = $ace_object->get_contig($split_contig_name,1);
    my $ct = Genome::Model::Tools::Pcap::ContigTools->new;
    my ($left_contig, $right_contig) = $ct->split($contig, $phd_object, split_position => $split_position, no_gui => $no_gui);

    $ace_object->remove_contig($contig->name);
    $ace_object->add_contig($left_contig);
    $ace_object->add_contig($right_contig);

    $ace_object->write_file(output_file => $out_file_name);
    

    return 1;
}

sub read_data_and_split_contig
{
    my ($self, $ace_file, $split_contig_name, $split_position, $no_gui, $out_file_name) = @_;

    my $ace_object = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file, using_db => 1);
    
    my $phd_object;    
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
    
    my $contig = $ace_object->get_contig($split_contig_name,1);
    my $ct = Genome::Model::Tools::Pcap::ContigTools->new;
    my ($left_contig, $right_contig) = $ct->split($contig, $phd_object, split_position => $split_position, no_gui => $no_gui);

    $ace_object->remove_contig($contig->name);
    $ace_object->add_contig($left_contig);
    $ace_object->add_contig($right_contig);

    $ace_object->write_file(output_file => $out_file_name);

}

1;





