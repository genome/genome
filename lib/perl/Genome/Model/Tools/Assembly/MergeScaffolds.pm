package Genome::Model::Tools::Assembly::MergeScaffolds;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::Tools::Assembly::MergeScaffolds
{
    is => 'Command',
    has => 
    [
        ace_file => {
            type => "String",
            optional => 0,
            doc => "This is the input ace file"
        }, 
        left_scaffold => {
            type => "String",
            optional => 0,
            doc => "This is the name of the left scaffold to be merged",
        },
        right_scaffold => {
            type => "String",
            optional => 0,
            doc => "This is the name of the right scaffold to be merged",
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
    merge-scaffolds --ace-file=in.ace --left-scaffold=Contig0.1 --right scaffold=Contig1.1 --out-file-name=out.ace
EOS
}

sub get_contig_names
{
    my ($self, $ace_file_name) = @_;
    my $fh = Genome::Sys->open_file_for_reading($ace_file_name);
    $self->error_message("There was an error opening ace file $ace_file_name for reading.") and die unless defined $fh;
    my @contig_names;
    while(my $line = <$fh>)
    {
        if($line =~ /^CO /)
        {
            my @tokens = split /\s+/,$line;
            push @contig_names, $tokens[1];            
        }
    }
    return \@contig_names;
}

sub get_scaffold_contigs
{
    my ($self, $contig_list, $scaffold_contig) = @_;
    
    my ($scaffold_name, $contig_num, $suffix);
    ($scaffold_name, $contig_num) = $scaffold_contig =~ /(Contig\d+)\.(\d+)/;
    ($suffix) = $scaffold_contig =~ /Contig\d+\.\d+(\D+)/;
    $suffix = '' unless defined $suffix;
    my @scaffold_contigs;
    @scaffold_contigs = grep 
    { 
        my $temp; 
        my $result = $_ =~ /$scaffold_name\.\d+$suffix/;
        ($temp) = $_ =~  /$scaffold_name\.\d+$suffix(\D+)/;#check to make sure that we don't count contigs with different extensions as belonging to the same scaffold.
        ($result && !(defined $temp && length $temp));
    } @{$contig_list};    
    
    return \@scaffold_contigs;
}

sub scaffold_comp 
{
    my ($c)= $a =~ /Contig\d+\.(\d+)/; 
    my ($d) = $b =~ /Contig\d+\.(\d+)/; 
    return $c <=>$d;
}

sub execute
{
    my $self = shift;
    my $ace_file = $self->ace_file;
    my $left_scaffold = $self->left_scaffold;
    my $right_scaffold = $self->right_scaffold;
    my $out_file_name = $self->out_file_name;

    $self->error_mesage("Input ace file $ace_file does not exist\n") and return unless (-e $ace_file);
    
    my $contig_names = $self->get_contig_names($ace_file);
    my $left_scaffold_contigs = $self->get_scaffold_contigs($contig_names, $left_scaffold);
    my $right_scaffold_contigs = $self->get_scaffold_contigs($contig_names, $right_scaffold);
    
    unless (@{$left_scaffold_contigs} && @{$right_scaffold_contigs}) 
    {
        print "Couldn't find scaffolds for $left_scaffold and/or $right_scaffold.\n";
        #print "It appears that scaffold $scaffold_name.XX$suffix does not exist.\n";
        return;
    }
    my @left_scaffold_contigs = sort scaffold_comp @$left_scaffold_contigs;
    my @right_scaffold_contigs = sort scaffold_comp @$right_scaffold_contigs;

    my $rightmost_contig = $left_scaffold_contigs[@left_scaffold_contigs-1];
    my ($scaffold_name, $rightmost_contig_num, $suffix);
    ($scaffold_name, $rightmost_contig_num) = $rightmost_contig =~ /(Contig\d+)\.(\d+)/;
    ($suffix) = $rightmost_contig =~ /Contig\d+\.\d+(\D+)/;
    $suffix = '' unless defined $suffix;

    my %right_scaffold_contigs;
    foreach my $contig (@right_scaffold_contigs)
    { 
        $rightmost_contig_num++;
        #print $contig,"\n";
        $right_scaffold_contigs{$contig} = $scaffold_name.'.'.$rightmost_contig_num.$suffix;       
    }
    #do a search and replace for all contig names above
    my $fh = IO::File->new($ace_file);
    my $out_fh = IO::File->new(">$out_file_name");
    while(my $line = <$fh>)
    {   
        if($line =~ /Contig/)
        { 
            foreach my $sub_contig (keys %right_scaffold_contigs)
            {
                if($line =~/$sub_contig\W+/)
                {                
                    $line =~ s/$sub_contig/$right_scaffold_contigs{$sub_contig}/g;
                    last;            
                }            
            }
        }
        $out_fh->print($line);
    }
    return 1;
}


1;





