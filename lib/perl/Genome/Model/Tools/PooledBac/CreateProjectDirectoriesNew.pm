package Genome::Model::Tools::PooledBac::CreateProjectDirectoriesNew;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Pcap::Ace;
use Genome::Model::Tools::Pcap::Phd;
use Genome::Model::Tools::PooledBac::Utils;

class Genome::Model::Tools::PooledBac::CreateProjectDirectoriesNew {
    is => 'Command',
    has => 
    [        
        pooled_bac_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "Pooled BAC Assembly Directory",    
        },
        ace_file_name =>
        {
            type => 'String',
            is_optional => 0,
            doc => "Ace file containing pooled bac contigs"
        },
        phd_file_name_or_dir =>
        {
            type => 'Sring',
            is_optional => 1,
            doc => "Phd file or dir containing read bases and quals"       
        },
        project_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "output dir for separate pooled bac projects"        
        },
        contig_map_file =>
        {
            type => 'String',
            is_optional => 1,
            doc => "this file contains a list of contigs and where they map to",
        }
    ]
};

sub help_brief {
    "Move Pooled BAC assembly into separate projects"
}   

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Move Pooled BAC Assembly into separate projects
EOS
}

############################################################
sub execute { 
    my $self = shift;
    print "Creating Project Directories...\n";
    my $pooled_bac_dir = $self->pooled_bac_dir;
    my $project_dir = $self->project_dir;
    my $phd_dir_or_ball = $self->phd_file_name_or_dir;
    $phd_dir_or_ball = $pooled_bac_dir.'/consed/phdball_dir/phd.ball.1' unless $phd_dir_or_ball;
    my $blastfile = $project_dir."/bac_region_db.blast";
    $self->error_message("$blastfile does not exist") and die unless (-e $blastfile);
    my $out = Genome::Model::Tools::WuBlast::Parse->execute(blast_outfile => $blastfile);   
    $self->error_message("Failed to parse $blastfile") and die unless defined $out;

    my $ace_file = $pooled_bac_dir.'/consed/edit_dir/'.$self->ace_file_name;
    $self->error_message("Ace file $ace_file does not exist") and die unless (-e $ace_file);
    my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => $ace_file, using_db => 1);
    $self->error_message("Failed to open ace file") and die unless defined $ao;
    my $po;

    
    my $contig_map_file = $self->contig_map_file || "CONTIG_MAP";
    $contig_map_file = $project_dir.'/'.$contig_map_file;
    $self->error_message("Contig map file, $contig_map_file, does not exist.\n") and die unless (-e $contig_map_file);
    
    my $ut = Genome::Model::Tools::PooledBac::Utils->create;
    $self->error_message("Genome::Model::Tools::PooledBac::Utils->create failed.\n") unless defined $ut;

    my $contig_map = $ut->open_contig_map($contig_map_file);
    my ($match_list, $orphan_list) = $ut->create_match_and_orphan_lists($contig_map);
    my $list = $ut->get_matching_contigs_list($out->{result});$out=undef;
    my %bac_contigs;
    foreach my $contig_name (keys %{$contig_map})
    {
        my $bac_name = $contig_map->{$contig_name}{maps_to};
        
        push @{$bac_contigs{$bac_name}},$contig_name;
    }
    
    foreach my $bac_name (keys %bac_contigs)
    {
        my $bac_dir = $project_dir."/$bac_name/";
        my @contig_names = @{$bac_contigs{$bac_name}};
        $self->error_message("Error creating directory $bac_dir") and die unless Genome::Sys->create_directory($bac_dir);
        my $old_dir = `pwd`;
        chdir($bac_dir);
        $ut->create_project_from_contig_names($ao,$bac_dir."/edit_dir/$bac_name.fasta.screen.ace", \@contig_names, $pooled_bac_dir);    
        chdir($old_dir);
    }
}


1;
