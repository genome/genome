package Genome::Model::Tools::ContaminationScreen::MegaBlast;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use Math::Round qw(round);

class Genome::Model::Tools::ContaminationScreen::MegaBlast
{
    is => 'Genome::Model::Tools::ContaminationScreen',
    has => [
        header => {
                           doc => 'db header to convert id to species (should correspond to db used)',
                           is => 'String',
                           is_input => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-M 15000000 -R 'select[type==LINUX64] rusage[mem=15000]'",
        }
    ],
};

sub help_brief 
{
    "locate contamination in megablast",
}

sub help_synopsis 
{
    return <<"EOS"
EOS
}

sub create 
{
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute 
{
    my $self = shift;
    my ($input_file,$database,$read_file) = ($self->input_file,$self->database,'megablast.out');
    my ($output_file) = ($self->output_file ? $self->output_file :  $input_file . '.screened');
    my $header = $self->header;
    my ($dir,$base) = (dirname($output_file), basename($output_file)); 

    #Step 0 - exec megablast
    my $raw_file = "$dir/0_$base.raw";
    my $mgb = "megablast -p 90.0 -D3 -F \"m D\" -fT -UT -R -d $database -i $input_file > $raw_file";
print "command:  $mgb\n";        
    $self->_run_cmd($mgb); 

    #Step 1 - exec filter
    my $filtered_file = "$dir/1_$base.filtered";
    my $fmt = "awk '(\$3>=98.0 && \$4>=50)||(\$3>=94.0 && \$4>=100)||(\$3>=90.0 && \$4>=200)' $raw_file > $filtered_file"; 
    $self->debug_message('Running: ' . $fmt);
    $self->_run_cmd($fmt); 

    #Step 2 - sort
    my $sorted_file = "$dir/2_$base.sorted";
    my $srt = "sort $filtered_file > $sorted_file";
    $self->debug_message('Running: ' . $srt);
    $self->_run_cmd($srt); 

    #Step 3 - get max
    #   presumes this column order
    #   0           1           2           3                   4           5               6           7       8           9       10          11
    #   Query id    Subject id  % identity  alignment length    mismatches  gap openings    q. start    q. end  s. start    s. end  e-value     bit score 

    my %maxes;
    my $reader = new FileHandle($sorted_file);
    <$reader>; #skip header

    while (my $line = <$reader>)
    {
        next if ($line=~/#/);
        my @columns = split("\t", $line);
        my $i = 0;
        my ($query, $subject, $e_value) = ($columns[0],$columns[1],$columns[10]);

        unless ($maxes{$query} and $maxes{$query}{e_value} > $e_value)
        {
            $maxes{$query} = {'subject' => $subject, 'e_value' => $e_value}; 
        }
    }    	
    close ($reader);
   
    #step 4 - get hits 
    my %hits;
    my $max_file = "$dir/3_$base.max";
    open (MAX, ">$max_file") or die("could not open $max_file");
    foreach my $key(sort keys %maxes)
    {
        $hits{$maxes{$key}{subject}}++;
        print MAX "$key\t" . $maxes{$key}{subject} . "\t" . $maxes{$key}{e_value} . "\n";
    }
    close(MAX); 

    #step 5 - translate id's to species
    system("rm -f $output_file");
    my ($val, $ind);
    foreach my $key(sort keys %maxes)
    {
        # regex assumes the header reference is of the format 
        #               >gi|147|emb|X05274.1| B.taurus mRNA for pancreatic trypsin inhibitor (BpTI)
        # where the id and species are separated by a space, so concatenate everything to the end of the line following the first space.
        # in the header being used for testing, found some lines like this:
        #               >gi|31341348|ref|NM_181029.2| Bos taurus casein alpha s1 (CSN1S1), mRNA >gi|175|emb|X00564.1| Bovine mRNA for pre-alpha S1-casein B
        # so terminate concatenation if '>' is encountered
        
        $val = $maxes{$key}{e_value};
        $key=~s/\|/\\|/g;
        $ind = "awk '/$key/{ for (i=2; i<=NF && !match(\$i, />.*/); i++) printf(\"%s \", \$i);print(\" = $val\");}' $header >> $output_file";
        system($ind);
    }
    
    $self->output_file($output_file);

    return 1;
}

sub _run_cmd()
{
    my ($self, $cmd) = @_;

    $self->debug_message('Running: '. $cmd);

    my $rv = system($cmd);
    unless ($rv == 0)
    {
        $self->error_message("non-zero return value($rv) from command $cmd");
        return;
    }
} 

sub _query()
{
    my ($self,$key) = @_;
    my $database = $self->database;

    my $cmd = "grep \"$key\" $database";
}

1;
