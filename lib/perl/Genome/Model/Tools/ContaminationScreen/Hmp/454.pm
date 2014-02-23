package Genome::Model::Tools::ContaminationScreen::Hmp::454;

use strict;
use warnings;

use Genome;
class Genome::Model::Tools::ContaminationScreen::Hmp::454
{
    is => 'Genome::Model::Tools::ContaminationScreen::Hmp',
    has => [
            parsed_file => {
                             doc => 'file to write parsed data to',
                             is => 'String',
                             is_output => 1,
                             is_optional => 1,
            },
            percent => {
                             doc => 'percent identity of alignment',
                             is => 'Number',
                             is_input => 1,
            },
            query_length => {
                             doc => 'length of query side of alignment',
                             is => 'Number',
                             is_input => 1,
            },
#            filter_list => 
#            {
#                           doc => 'file of id\'s produced from screening, to remove from original fasta',
#                           is => 'String',
#                           is_input => 1,
#                           is_optional => 1,
#            },
    ],
};

sub help_brief 
{
    "locate human contamination in 454",
}

sub help_synopsis 
{
    return <<"EOS"
    gt cross_match.test <input reads> <database> -raw -tags -minmatch 14 -bandwidth 6 -penalty -1 -gap_init -1 -gap_ext -1 
EOS
}

sub create 
{
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    print "454 created \n";
    return $self;
}

sub execute 
{
    #$DB::single = 1;
    my $self = shift;
    my $input_file = $self->input_file;
    my $output_file = $self->output_file ? $self->output_file :  $input_file . '.screened';

    # Step 1 - create read file
    my $cmd = 'cross_match.test ' . $input_file . ' ' .  $self->database . ' -raw -tags -minmatch 14 -bandwidth 6 -penalty -1 -gap_init -1 -gap_ext -1 > ' . $output_file;
    $self->debug_message('Running: '. $cmd);
    print "454 step 1 - cmd:  $cmd\n";
    my $rv = system($cmd);

    $self->output_file($output_file);

    # Step 2 - parse 
    # from: /gscmnt/233/analysis/sequence_analysis/scripts/parse_crossmatch_results.pl
    my $parsed_file = $self->parsed_file ? $self->parsed_file :  $input_file . '.parsed';
    my $filter_list = $self->filter_list ? $self->filter_list :  $input_file . '_FILTERED.txt';
    my $success_marker = $filter_list . ".success";
    my $screened = new IO::File $output_file;
    my $percent = $self->percent;
    my $length = $self->query_length;
    my $data;

    my $tmp_counter = 0;

    print "step 2.0 - starting loop\n";
    while (<$screened>) 
    {
        chomp;
        my $line = $_;
        next if ($line =~ /^\s*$/);
        next unless ($line =~ /^ALIGNMENT/);
        my @line = split(/\s+/,$line);
        my $query_name             = $line[5];
        my $score                  = $line[1];
        my $perc_subs              = $line[2] / 100;
        my $perc_dels              = $line[3] / 100;
        my $perc_ins               = $line[4] / 100;
        my $query_alignment_length = $line[7] - $line[6] + 1;

        # Calculate length of ALIGNMENT based on length of query-side of alignment & info about query-side deletions
        my $num_dels               = int(($query_alignment_length * $perc_dels) + .5);
        my $alignment_length       = $query_alignment_length + $num_dels;

        # percent identity relative to full alignment length
        my $num_subs               = int(($query_alignment_length * $perc_subs) + .5);
        my $num_ins                = int(($query_alignment_length * $perc_ins) + .5);
        my $total_errors           = $num_subs + $num_ins + $num_dels;
        my $total_correct          = $alignment_length - $total_errors;
        my $percent_identity       = int((($total_correct / $alignment_length) * 10000) + .5) / 100; #round to 1e-2
    
        # if this alignment line meets the cutoff criteria
        if (($percent_identity >= $percent) and ($alignment_length >= $length)) 
        {
	    if (defined $data->{$query_name}) 
            {
	        if ($score > $data->{$query_name}->{score}) 
                {
		    $data->{$query_name}->{score}            = $score;
    		    $data->{$query_name}->{percent_identity} = $percent_identity;
		    $data->{$query_name}->{alignment_length} = $alignment_length;
		    $data->{$query_name}->{line}             = $line;
		    $data->{$query_name}->{subs}             = $num_subs;
		    $data->{$query_name}->{ins}              = $num_ins;
		    $data->{$query_name}->{dels}             = $num_dels;
	        }
	    } 
            else 
            {
	        $data->{$query_name}->{score}            = $score;
	        $data->{$query_name}->{percent_identity} = $percent_identity;
	        $data->{$query_name}->{alignment_length} = $alignment_length;
	        $data->{$query_name}->{line}             = $line;
	        $data->{$query_name}->{subs}             = $num_subs;
	        $data->{$query_name}->{ins}              = $num_ins;
	        $data->{$query_name}->{dels}             = $num_dels;
	    }
        }
    }
    print "2.1 - loop completed\n";
    #Dump parsed in order of best score to least of all the alignments meeting the cutoff
    my $parsed = new IO::File ">>$parsed_file";
    my $filtered = new IO::File ">>$filter_list";
    print $parsed "<tag> <score> <%sub> <%del> <%ins> <query> <qstart> <qend> <unaligned length> *<\'C\' or null> <subject> **<sstart> **<send> **<unaligned length>\n";
    print $parsed "(NOTE: *the \'C\' tag may or may not be present (the whole column is deleted...), if present this alignment is to the minus strand)\n";
    print $parsed "(NOTE: **These last 3 columns will be flipped when the alignment is to the minus strand (i.e. the order will be <unaligned length> <send> <sstart>))\n";
    foreach my $query_name (sort{$data->{$b}->{score} <=> $data->{$a}->{score}} keys(%{$data})) 
    {
        print $parsed "========== $query_name (perc_id:$data->{$query_name}->{percent_identity} align_length:$data->{$query_name}->{alignment_length} subs:$data->{$query_name}->{subs} dels:$data->{$query_name}->{dels} ins:$data->{$query_name}->{ins}) ==========\n";
        print $parsed "$data->{$query_name}->{line}\n";

        #print filtered read
        print $filtered "$query_name\n";
    } 
    $parsed->close;
    $filtered->close;
    $self->parsed_file($parsed_file);
    $self->filter_list($filter_list);

    #per Idas - action to confirm completion
    my $mark_status = "touch $success_marker";
    system($mark_status);
    print "contamination screen complete\n";
    return 1;
}

1;
