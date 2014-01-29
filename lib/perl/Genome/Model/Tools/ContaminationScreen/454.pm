package Genome::Model::Tools::ContaminationScreen::454;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Math::Round qw(round);

class Genome::Model::Tools::ContaminationScreen::454
{
    is => 'Genome::Model::Tools::ContaminationScreen',
};

sub help_brief 
{
    "locate contamination in quality/vector trimmed data for 454",
}

sub help_synopsis 
{
    return <<"EOS"
    gmt cross_match.test --input_file --database -raw -tags -minmatch 14 -bandwidth 6 -penalty -1 -gap_init -1 -gap_ext -1 >  --read_file
    gmt --parse_script -input --read_file -output --parsed_file -percent 90 -length 50 > --output_file
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
    my ($read_file, $parsed_file, $hits_file) = ($input_file . '.reads', $input_file . '.reads.parsed', $input_file . '.hits.fna');
    my ($output_file, $summary_file) = ($self->output_file ? $self->output_file :  $input_file . '.screened', 
                                        $self->summary_file ? $self->summary_file : $input_file . '.summary');
    my $parse_script = '/gscmnt/233/analysis/sequence_analysis/scripts/parse_crossmatch_results.pl';

    #create read file
    my $cmd = 'cross_match.test ' . $input_file . ' ' .  $self->database . ' -raw -tags -minmatch 14 -bandwidth 6 -penalty -1 -gap_init -1 -gap_ext -1 > ' . $read_file;
    $self->debug_message('Running: '. $cmd);
    my $rv = system($cmd);
    unless ($rv == 0) {
        $self->error_message("non-zero return value($rv) from command $cmd");
        return;
    }

    #run parsing script
    my $parse_cmd = $parse_script . ' -input ' . $read_file . ' -output ' . $parsed_file . ' -percent 90 -length 50'; 
    $self->debug_message('Running: ' . $parse_cmd);
    $rv = system($parse_cmd);
    unless ($rv == 0)
    {
        $self->error_message("non-zero return value($rv) from command $parse_cmd");
        return;
    }

    #search for hits and parse out column reads
    my ($grep, $parse, $summary) =  (new IO::File, new IO::File, new IO::File);
    my ($line, @align);
    my ($total_reads, $contaminated_reads) = (0,0);
    my %aligns;
    
    if ($grep->open("< $parsed_file"))
    {
        if ($parse->open("> $output_file"))
        {
            do
            {
                chomp($line = $grep->getline());
                if ($line=~/====/)
                {
                    @align = split(/ /,$line);
                    $aligns{$align[1]}++; #consume duplicates, so only uniques are written to output
                    $contaminated_reads++;
                }
                $total_reads++;
            }
            until ($grep->eof());
            print $parse join("\n", keys(%aligns));
            $parse->close;
            $grep->close;
        }
        else
        {
            $self->error_message("error opening $output_file for writing:  " . $!);
        }
    }
    else
    {
        $self->error_message("error opening $read_file for reading");
    }
    if ($summary->open("> $summary_file"))
    {
       print $summary "total reads:  $total_reads\n";
       print $summary "contaminated reads $contaminated_reads\n";
       print $summary round($contaminated_reads * 100/$total_reads) . "%" if $total_reads;  
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
