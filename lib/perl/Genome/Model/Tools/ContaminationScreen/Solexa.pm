package Genome::Model::Tools::ContaminationScreen::Solexa;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Math::Round qw(round);

class Genome::Model::Tools::ContaminationScreen::Solexa
{
    is => 'Genome::Model::Tools::ContaminationScreen',
    has_input => [
                    minscore => {
                                    is => 'Number', doc => 'value for minimum score, based on read length',
                    },
        ],
    has => [
    ],
    has_param => [
            lsf_resource => {
                             default_value => "-M 15000000 -R 'select[type==LINUX64] rusage[mem=15000]'",
            },
    ],
};

sub help_brief 
{
    "locate contamination in GERALD fastq files for Illumina/Solexa",
}

sub help_synopsis 
{
    return <<"EOS"
    gmt cross_match.test --input_file --database -raw -tags -minmatch 16 -minscore -bandwidth 3 -penalty -1 -gap_init -1 -gap_exp -1 >  --output_file 
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
    my $input_file = $self->input_file;
    my ($read_file, $minscore) = ($input_file . '.reads', $self->minscore);
    my ($output_file, $summary_file) = ($self->output_file ? $self->output_file :  $input_file . '.screened', 
                                        $self->summary_file ? $self->summary_file : $input_file . '.summary');
    #create read file
    my $cmd = 'cross_match.test ' . $input_file . ' ' .  $self->database . ' -raw -tags -minmatch 16 -minscore ' . $minscore . ' -bandwidth 3 -penalty -1 -gap_init -1 -gap_ext -1 > ' . $read_file;
    $self->debug_message('Running: '. $cmd);
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return value($rv) from command $cmd");
        return;
    }

    #search for hits and parse out column reads
    my ($grep, $parse, $summary) =  (new IO::File, new IO::File, new IO::File);
    my ($line, @align);
    my ($total_reads, $contaminated_reads) = (0,0);
    my ($hit,%aligns);

    if ($grep->open("< $read_file"))
    {
        if ($parse->open("> $output_file"))
        {
            do
            {
                chomp($line = $grep->getline());
                if ($line=~/ALIGNMENT/)
                {
                    @align = split(/ /, $line);
                    $align[9]=~/(.*\/)(\d)/;
                    
                    #screening paired ends - 
                    #if we get HWI-EAS90:4:1:4:388#0/1 
                    #we also need HWI-EAS90:4:1:4:388#0/2, and vice versa
                    #consume duplicates, so only uniques are written to output
                    $aligns{$1 . 1}++;
                    $aligns{$1 . 2}++;
                    $contaminated_reads++;
                }
                $total_reads++;
            }
            until ($grep->eof());
            print $parse join("\n", sort keys(%aligns));
            $parse->close; 
            $grep->close;
        }
        else
        {
            $self->error_message("error opening $output_file for writing");
        }
    }
    else
    {
        $self->error_message("error opening $read_file for reading");
    }   

    #make summary file
    if ($summary->open("> $summary_file"))
    {
       print $summary "total reads:  $total_reads\n";
       print $summary "contaminated reads $contaminated_reads\n";
       print $summary round($contaminated_reads * 100/$total_reads) . "%";  
    }
    else
    {
       $self->error_message("error opening $summary_file for writing"); 
    } 

    $self->output_file($output_file);
    $self->summary_file($summary_file);
    return 1;
}

1;
