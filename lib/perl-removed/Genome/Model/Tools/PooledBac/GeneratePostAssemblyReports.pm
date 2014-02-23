package Genome::Model::Tools::PooledBac::GeneratePostAssemblyReports;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::Pcap::Ace;
use Genome::Model::Tools::Pcap::Phd;
use Genome::Model::Tools::PooledBac::Utils;
use Genome::Sys;
use Cwd;
use List::Util qw(max min);

class Genome::Model::Tools::PooledBac::GeneratePostAssemblyReports {
    is => 'Command',
    has => 
    [        
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

sub print_assembly_size_report
{
    my ($self, $project_names, $report_name) = @_;
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $name (@{$project_names})
    {
        my $size = 0;        
        if (-e "$name/edit_dir/$name.fasta.screen.ace")
        {
            my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => "$name/edit_dir/$name.fasta.screen.ace", using_db => 1);
            if ( not $ao ) {
                $self->error_message("Could not get ace object for project: $name");
                return;
            }
            my $contig_names = $ao->get_contig_names;
            
            foreach my $contig_name (@{$contig_names})
            {
                my $contig = $ao->get_contig($contig_name);
                my $length = length ($contig->unpadded_base_string);
                $size += $length;            
            }            
        }
        $fh->print ("$name $size\n");
    }
}

sub print_contig_size_report
{
    my ($self, $project_names, $report_name) = @_;
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $name (@{$project_names})
    {
        my $size = 0;
        $fh->print ("$name:\n");
        if (-e "$name/edit_dir/$name.fasta.screen.ace")
        {
            my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file => "$name/edit_dir/$name.fasta.screen.ace", using_db => 1);
            my $contig_names = $ao->get_contig_names;
            
            foreach my $contig_name (@{$contig_names})
            {
                my $contig = $ao->get_contig($contig_name);
                my $size = length ($contig->unpadded_base_string);
                $fh->print("$contig_name $size\n");            
            }
        }
        else
        {
            $fh->print( "NO CONTIGS\n");
        }
        $fh->print("\n");        
    }
}

sub create_fasta_and_qual_files
{
    my ($self, $project_names) = @_;
    foreach my $name (@{$project_names})
    {             
        if (-e "$name/edit_dir/$name.fasta.screen.ace")
        {
            $self->ace2fastaqual("$name/edit_dir/$name.fasta.screen.ace", "$name/edit_dir/contigs.fasta");
        }      
    }
}

sub print_contigs_only_consensus_report
{
    my ($self, $project_names, $report_name) = @_;
    my $fh = IO::File->new('>'.$report_name);
    $self->error_message("Failed to open $report_name for writing.")  and die unless defined $fh;
    foreach my $name (@{$project_names})
    {
        my $size = 0;
        $fh->print ("$name:\n");
        if (-e "$name/edit_dir/$name.fasta.screen.ace")
        {
            my $ao = Genome::Model::Tools::Pcap::Ace->new(input_file=>"$name/edit_dir/$name.fasta.screen.ace", using_db => 1);
            my $contig_names = $ao->get_contig_names;
            foreach my $contig_name (@{$contig_names})
            {
                my $contig = $ao->get_contig($contig_name);
                my $read_hash = $contig->reads;
                my $consensus_only = 1;
                foreach my $read_name (keys %{$read_hash})
                {
                    unless($read_name =~ /.*\.c1$/)
                    {
                        $consensus_only = 0;
                        last;
                    }                
                }  
                if($consensus_only)
                {
                    $fh->print ($contig_name,"\n");
                }          
            }
            system "/bin/rm $name/edit_dir/$name.fasta.screen.ace.db" if -e "$name/edit_dir/$name.fasta.screen.ace.db";
        }
        else
        {
            $fh->print( "NO CONTIGS\n");
        }
        $fh->print("\n");        
    }
}

sub get_bac_names
{
    my ($self) = @_;
    my $fh = IO::File->new('ref_seq.fasta');
    $self->error_message("Failed to open ref_seq.fasta for reading.")  and die unless defined $fh;
    my @names;
    while (my $line = <$fh>)
    {
        chomp $line;
        my ($name) = $line =~ /^\>(.*)/;
        push @names,$name if defined $name;    
    }
    return \@names;
}

sub ace2fastaqual
{
    my ($self,$infile, $outfile) = @_;

    my $reader = Genome::Model::Tools::Consed::AceReader->create(file => $infile);
    $self->error_message("Error creating ace reader for $infile.") and die unless defined $reader;
    my $outfh = IO::File->new(">$outfile");
    $self->error_message("Error opening $outfile.") and die unless defined $outfh;
    my $outfh2 = IO::File->new(">$outfile.qual");
    $self->error_message("Error opening $outfile.qual.") and die unless defined $outfh2;
    while(my $line = $reader->_fh->getline)
    {
        if($line =~ /^CO/)
        {
            $reader->_fh->seek(-length($line),1);
            my $item = $reader->next;    
            if($item->{type} eq 'contig')
            {
                $outfh->print(">",$item->{name},"\n");
                $outfh2->print(">",$item->{name},"\n");
                $item->{consensus} =~ tr/Xx/Nn/;
                $item->{consensus} =~ s/\*//g;

                $outfh->print($item->{consensus},"\n");
                $outfh2->print($item->{base_qualities},"\n");
            }
        }
    }
}

sub _remove_contigs_idx_directories {
    my ( $self, $bac_names ) = @_;

    foreach my $bac_name ( @{$bac_names} ) {
	my $idx_dir = $self->project_dir."/$bac_name/edit_dir/$bac_name".'.fasta.screen.ace.idx';

	if ( -d $idx_dir ) {
	    unless ( File::Path::rmtree( $idx_dir ) ){
		$self->debug_message("Attempted to remove contigs index directory but failed: $idx_dir");
	    }
	}
    }
    return 1;
}
############################################################
sub execute { 
    my $self = shift;
    print "Generating Post Assembly Reports...\n";
    my $orig_dir = cwd();
    my $project_dir = $self->project_dir;
    chdir($project_dir);
    my $reports_dir = $project_dir."/reports/";
    $self->error_message("Failed to create directory $reports_dir")  and die unless Genome::Sys->create_directory($reports_dir);
    
    my $ut = Genome::Model::Tools::PooledBac::Utils->create;
    $self->error_message("Genome::Model::Tools::PooledBac::Utils->create failed.\n") unless defined $ut;

    my $contig_map_file = $self->contig_map_file || "CONTIG_MAP";
    $contig_map_file = $project_dir.'/'.$contig_map_file;    
    my $contig_map = $ut->open_contig_map($contig_map_file);
    
    my $names = $self->get_bac_names;
    $self->print_assembly_size_report($names,$reports_dir."assembly_size_report");
    $self->print_contig_size_report($names, $reports_dir."contig_size_report");
    $self->print_contigs_only_consensus_report($names, $reports_dir."contigs_only_consensus");
    $self->create_fasta_and_qual_files($names);
    $self->_remove_contigs_idx_directories( $names ); 
    chdir( $orig_dir );
    return 1;
}



1;
