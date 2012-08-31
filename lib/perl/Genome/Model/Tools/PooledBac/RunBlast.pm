package Genome::Model::Tools::PooledBac::RunBlast;

use strict;
use warnings;

use Genome;

use Cwd;
use IO::File;

class Genome::Model::Tools::PooledBac::RunBlast {
    is => 'Command',
    has => 
    [ 
        ref_seq_file => 
        {
            type => 'String',
            is_optional => 0,
            doc => "File pointing to ref seq locations",
        },
        pooled_bac_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "Pooled BAC Assembly Directory",    
        },    
        ace_file_name => 
        {
            type => 'String',
            is_optional => 1,
            doc => "ace file containing pooled bac sequences"
        },
        pooled_bac_fasta_file =>
        {
            type => 'String',
            is_optional => 1,
            doc => "fasta file containing pooled bac sequences"        
        },
        project_dir =>
        {
            type => 'String',
            is_optional => 0,
            doc => "output directory for pooled bac projects"        
        },
        ref_qual_value =>
        {
            type => 'String',
            is_optional => 1,
            doc => "This is the quality value that is used when creating reference .qual files, the default is 37",        
        },
    ]    
};

sub help_brief {
    "The first step of creating Pooled BAC Projects Maps Contigs to a known reference sequence"
}

sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Blasts Pooled BAC contigs against a known reference assembly.
EOS
}

############################################################
sub execute { 
    my $self = shift;
    print "Running Blast...\n";
    my $ref_seq_file = $self->ref_seq_file;
    my $pooled_bac_dir = $self->pooled_bac_dir;
    my $project_dir = $self->project_dir;
    $self->error_message("Error creating directory $project_dir") and die unless Genome::Sys->create_directory($project_dir);
    
    my $ace_file = ''; 
    my $fasta_file = '';
    my $orig_dir = cwd();
    chdir($self->pooled_bac_dir);
    #build db file
    #input that gets converted to query fasta file
    $ace_file = $self->pooled_bac_dir.'/consed/edit_dir/'.$self->ace_file_name if($self->ace_file_name); 
    #if a pooled bac fasta file is provided, we use it instead of creating fasta from the ace file above
    $fasta_file = $self->project_dir.'/'.$self->pooled_bac_fasta_file if($self->pooled_bac_fasta_file); 
    $self->error_message("Need either an ace file or pooled bac fasta file to be specified.\n") and die unless (-e $ace_file || -e $fasta_file);
    my $ref_fasta_file =$self->project_dir.'/ref_seq.fasta';
    #build fasta containing reference sequence regions to be blasted against
    $self->build_fasta_file($self->ref_seq_file) if(!(-e $ref_fasta_file&&-e "$ref_fasta_file.qual"));
    
    chdir($project_dir);
    Genome::Model::Tools::WuBlast::Xdformat::Create->execute(
        database => 'bac_region_db', 
        fasta_files => [$self->project_dir.'/ref_seq.fasta'],        
    ); 
    #build query file
    my $query_fasta_file = $self->project_dir.'/pooled_contigs.fasta';
    if(-e $ace_file )
    {   
        $self->ace2fasta($ace_file,$query_fasta_file);
    }
    else
    {
        system("cp $fasta_file $query_fasta_file");
    }
    my $params = 'M=1 N=-3 R=3 Q=3 W=30 wordmask=seg lcmask hspsepsmax=1000 B=1 V=1 topcomboN=1 -errors -notes -warnings -cpus 4 2>/dev/null';
    
    Genome::Model::Tools::WuBlast::Blastn->execute(
        database => 'bac_region_db',
        query_file => $query_fasta_file,
        params => $params,
    );

    chdir( $orig_dir );
    
    return 1;
}

sub parse_ref_seq_coords_file
{
    my ($self) = @_;
    
    my $ref_coords_file = $self->ref_seq_file;
    $self->error_message("$ref_coords_file does not exist\n") and die unless -e $ref_coords_file;
    my $fh = IO::File->new($ref_coords_file);
    $self->error_message("Error opening $ref_coords_file.\n") and die unless defined $fh;
    
    my %ref_seq_coords;
    my $name;
    my $line;
    my $lc = 0; #line count
    while ($line = $fh->getline)
    {
        $lc++;
        chomp $line;
        last if($line =~ /REF_SEQ_COORDINATES/);
        next unless length($line);
        next if($line=~/^\#/);
        next if($line=~/^\s*$/);
        my ($token, $data) = $line =~ /(.*)\:\s*(.*)/;
        $ref_seq_coords{$token} = $data;
        
    }    
    my $bac_name;
    $self->error_message ("Could not find any coordinates, ref seq coordinates should be preceded by the heading 'REF_SEQ_COORDINATES', followed by a list of 1 or more ref seq coordinates.\n") and die unless ($line =~ /REF_SEQ_COORDINATES/);
    while ($line = $fh->getline) 
    {
        $lc++;
        chomp $line;
        next unless length($line);
        next if($line=~/^\#/);
        next if($line=~/^\s*$/);
        my @tokens = split/\s+/,$line;
        unless ((scalar @tokens) == 4|| (scalar @tokens) == 0)
        {
            $self->error_message("Line number $lc appears to be formatted incorrectly.  Bac coordinate lines should use the following space-delimited format\nBAC_NAME CHROMOSONE_NAME START_POSITION END_POSITION\n") and die;
        }
        
        $bac_name = $tokens[0];
        my %hash;
        @hash{ 'chromosome','start','end'} = @tokens[1..3];
        $self->error_message("line $lc: BAC name, $bac_name, appears to be invalid\n") and die unless ($bac_name=~/^\S+$/);        
        $self->error_message("line $lc: Chromosome name, $tokens[1], appears to be invalid, valid names are 1-22, X, and Y.\n") and die unless ($tokens[1] =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)$/);
        $self->error_message("line $lc: Start position, $tokens[2], is not numeric\n") and die unless ($tokens[2] =~/^\d+$/ );
        $self->error_message("line $lc: End position, $tokens[3], is not numeric\n") and die unless ($tokens[3] =~/^\d+$/ );
        $self->error_message("line $lc: Start and end position need to be positive integers\n") and die unless($tokens[2] >0 && $tokens[3] >0);
        $self->error_message("line $lc: Start position needs to be less than or equal to the end position\n") and die unless($tokens[2] <= $tokens[3]);
        $ref_seq_coords{$bac_name} = \%hash;        
    }
    $self->error_message ("Could not find any coordinates, ref seq coordinates should be preceded by the heading 'REF_SEQ_COORDINATES', followed by a list of 1 or more ref seq coordinates.\n") and die unless(scalar (keys %ref_seq_coords));
    
    return \%ref_seq_coords;
}

sub ace2fasta
{
    my ($self,$infile, $outfile) = @_;

    my $reader = Genome::Model::Tools::Consed::AceReader->create(file => $infile);
    $self->error_message("Error creating ace reader for $infile.") and die unless defined $reader;
    my $outfh = IO::File->new(">$outfile");
    $self->error_message("Error opening $outfile.") and die unless defined $outfh;
    while(my $line = $reader->_fh->getline)
    {
        if($line =~ /^CO/)
        {
            $reader->_fh->seek(-length($line),1);
            my $item = $reader->next;    
            if($item->{type} eq 'contig')
            {
                $outfh->print(">",$item->{name},"\n");
                $item->{consensus} =~ tr/Xx/Nn/;
                $item->{consensus} =~ s/\*//g;

                $outfh->print($item->{consensus},"\n");
            }
        }
    }
}

sub get_seq
{
    my ($self, $ref_seq_dir,$bac_name, $chromosome, $ref_start, $ref_stop) = @_;
    my $ref_seq_fasta = $ref_seq_dir."/$chromosome.fasta";
    my $fh = IO::File->new($ref_seq_fasta);
    $self->error_message("Error opening $ref_seq_fasta") and die unless defined $fh;
    my $line = <$fh>;
    my $seq_string;
    while(my $line = <$fh>)
    {
        chomp $line;
        $seq_string .= $line;
    }
    return substr($seq_string, $ref_start, $ref_stop-$ref_start+1);    
}
sub write_fasta
{
    my ($self, $fh, $qfh,$ref_seq_dir, $bac_name, $chromosome, $ref_start, $ref_stop) = @_;
    my $bac_sequence = $self->get_seq($ref_seq_dir,$bac_name, $chromosome,$ref_start,$ref_stop);
    $fh->print(">$bac_name\n");
    $fh->print("$bac_sequence\n");
    $qfh->print(">$bac_name\n");
    my $ref_qual = $self->ref_qual_value || '37';
    my $qual_string;
    for(my $i=0;$i<length($bac_sequence);$i++)
    {
        $qual_string.="$ref_qual ";
    }
    $qual_string.="\n";
    $qfh->print($qual_string);
}

sub build_fasta_file
{
    my ($self, $ref_seq_file) = @_;
    my $data = $self->parse_ref_seq_coords_file;
    
    my $ref_seq_fasta = $self->project_dir.'/ref_seq.fasta';
    my $ref_seq_dir = delete $data->{REFERENCE_ASSEMBLY};
    
    my $fh = IO::File->new(">$ref_seq_fasta");
    $self->error_message("Failed to open $ref_seq_fasta for writing.") and die unless defined $fh;
    my $qfh = IO::File->new(">$ref_seq_fasta.qual");
    $self->error_message("Failed to open $ref_seq_fasta.qual for writing.") and die unless defined $qfh;
    
    foreach my $bac_name (keys %{$data})
    {
        my ($chromosome, $ref_start, $ref_stop) =  map {  $data->{$bac_name}{$_} } ('chromosome','start','end');#@{{$data}->{$bac_name}};
        $self->write_fasta($fh,$qfh,$ref_seq_dir,$bac_name,$chromosome, $ref_start, $ref_stop);
    }
}


1;
