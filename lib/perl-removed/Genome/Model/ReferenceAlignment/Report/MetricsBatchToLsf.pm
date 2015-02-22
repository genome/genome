package Genome::Model::ReferenceAlignment::Report::MetricsBatchToLsf;

use strict;
use warnings;

use Genome;

#use Genome::Model::Command::Report::Variations;
use Command; 
use PP::LSF;
use File::Basename;

class Genome::Model::ReferenceAlignment::Report::MetricsBatchToLsf 
{
    is => 'Command',                       
    has => 
    [ 
        out_log_file => 
        {
            type => 'String',
            is_optional => 1,
            doc => "Will combine the LSF output files into this file",
        },
        error_log_file => 
        {
            type => 'String',
            is_optional => 1,
            doc => "Will combine the LSF error files into this file",
        },
        input => {
            type => 'String',
            doc => 'File path for input map',
        },  
        snpfile => {
            type => 'String',
            doc => 'File path for input snp file',
        },
        qual_cutoff => {
            type => 'int',
            doc => 'quality cutoff value', 
        },
        output => {
            type => 'String',
            doc => 'File path output metrics file', 
            is_optional => 0,
        },
        chunk_count => {
            type => 'int',
            doc => 'number of jobs to run', 
            is_optional => 1,
            default => 20,
        },   
    ],
};

#########################################################

sub help_brief {
    return "Generates an annotation report for a variation file."
}

sub help_synopsis {
    return;
}

sub help_detail {
    return <<"EOS"
    For a given variant file, this module will split it up into files containing 10000 variants.  Then for each of these split files, a job will be spawned and monitored that will run Genome::Model::Command::Report::Variation.  When the jobs are completed, the resulting report files will be combined into one report.  If requested the output and error logs from lsf can also be combined.
EOS
}

sub create {
    use Carp;
    confess "Moved this module from G:M:C:Report to G:M:RefAlign:Report and needs to be updated";
}

sub execute {
    my $self = shift;
    
    # create child jobs
    my $jobs = $self->_create_jobs;
    $self->error_message("No jobs created")
        and return unless $jobs and @$jobs;

    # run & monitor jobs
    my $success = $self->_run_and_monitor_jobs($jobs);

    # finish
    return $self->_finish($success, $jobs);
}

sub snp_chunk
{
    my ($self, $input_file,$count) = @_;
    my $in = IO::File->new($input_file);
    my $chunk_prefix = $self->snp_chunk_prefix($input_file);
    my @out;

    for(my $i = 0;$i<$count;$i++)
    {
        $out[$i] = IO::File->new(">$chunk_prefix.$i");
    }
    my $i=0;
    #I give line 0 to chunk file 0, line 1 to chunk file 1, etc.
    #I divide it this way because certain clusters of SNP's have quite a few overlapping
    #reads, and I want that work to be distributed
    #The other method of chunking, which would involve give the first 10% of the SNP file
    #to chunk 0, second 10% to chunk file 2, etc., would result in some chunks having a much
    #greater amount of work than others.
    while(my $line = <$in>)
    {
        while($i>=$count){$i-=$count;}#same as $i= $i%$count, but a bit more efficient, as long as $i isn't some huge number
        $out[$i]->print( $line);
        $i++;
    }

    for(my $i = 0;$i<$count;$i++)
    {
        $out[$i]->close;
    }
    return $chunk_prefix;
}

sub snp_chunk_prefix
{
    my ($self, $input_file) = @_;
    my $out_path = dirname($input_file);
    my $out_prefix = basename $input_file;
    my $out_dir = $out_path.'/'.$out_prefix.".chunk";
    `mkdir $out_dir` unless -d $out_dir;
    return "$out_dir/$out_prefix";
}

sub snp_unchunk
{
    my ($self,$prefix,$count,$output) = @_;

    $prefix = dirname($prefix.'.0').'/'.basename($prefix.'.0','.0');
    my @fhs;
    for(my $i=0;$i<$count;$i++)
    {
        push @fhs,IO::File->new($prefix.".$i");
    }
    my $out_fh = IO::File->new(">$output");
    my $line_no = 0;
    my $line;
    while(1)
    {
        if($line_no>$count){$line_no-=$count;}
        last unless($line =$fhs[$line_no%$count]->getline);   
        print $out_fh $line;
        $line_no++;
    }
}


#call chunk function here, and setup chunks
sub _create_jobs
{
    my $self = shift;    
    
    my $chunk_prefix = $self->snp_chunk($self->snpfile,$self->chunk_count);

    my @jobs;
    
    for(my $i=0;$i<$self->chunk_count;$i++)
    {
        unless (push @jobs, $self->_setup_job($i))    
        {
            $self->_finish(0, \@jobs);
            return;
        }
    }
    return \@jobs;    

}

#create jobs here, put params here
sub _setup_job
{
    my ($self, $num) = @_;

    my $report_file_base = sprintf
    (
        '%s.%d', 
        $self->output,
        $num,
    );

    # If logging, get a log file for each job
    my ($out_file, $error_file);
    if ( $self->out_log_file )
    {
        $out_file = sprintf
        (
            '%s.%d', 
            $self->out_log_file, 
            $num,
        );
        unlink $out_file if -e $out_file;
    }

    if ( $self->error_log_file )
    {
        $error_file = sprintf
        (
            '%s.%d', 
            $self->error_log_file, 
            $num,
        );
        unlink $error_file if -e $error_file;
    }
    
    my $prefix = $self->snp_chunk_prefix($self->snpfile);
    my $out_prefix = $self->snp_chunk_prefix($self->output);    
    my $snpfile =  "$prefix.$num";
    my $outfile = "$out_prefix.$num";
    my $hostname = `hostname -s`; 
    chomp $hostname;  
    my %job_params =
    (
        pp_type => 'lsf',        
        q => 'aml',
        #R => "'select[db_dw_prod_runq<10] rusage[db_dw_prod=1]'",
        #R => "'hname!=$hostname && hname!=linuscs50'",#exclude linuscs50
        command => sprintf
        (
            '`which gmt` maq generate-variation-metrics --input "%s" --snpfile %s --qual-cutoff 1 --output %s',
            #'/bin/touch %s',
            $self->input,
            $snpfile,
            $outfile
        ),
    );

    # print $job_params{command},"\n";

    $job_params{o} = $out_file if $out_file;
    $job_params{e} = $error_file if $error_file;

    my $job = PP::LSF->create(%job_params);
    $self->error_message("Can't create job: $!")
        and return unless $job;

    return 
    {
        job => $job,
        snpchunk=> $snpfile,
        outchunk => $outfile,
        out => $out_file,#out log file, needs less confusing name...
        #error => $error_file,
        num => $num,
        tries => 3,
        try_count => 1,
    };
}

sub _run_and_monitor_jobs
{
    my ($self, $jobs) = @_;

    # Start jobs.  To monitor, create hash w/ job ids as keys.
    my %running_jobs;
    for my $num ( 0..(scalar(@$jobs) - 1) )
    {
        # Set local $job for clarity
        my $job = $jobs->[$num]->{job};
        $job->start;
        $running_jobs{ $job->id } = $num;
    }

    no warnings; # these are filling the log files, not sure from where
    # Monitor
    MONITOR: while ( %running_jobs )
    {
        sleep 30;
        for my $job_id ( keys %running_jobs )
        {
            # Set local $job for clarity
            my $job = $jobs->[ $running_jobs{$job_id} ]->{job};
            my $job_hash = $jobs->[ $running_jobs{$job_id} ];
            if ( $job->has_ended )
            {
                if ( $job->is_successful )
                {
                    print "$job_id successful\n";
                    delete $running_jobs{$job_id};
                }
                elsif($job_hash->{try_count} <= $job_hash->{tries} )
                {
                    print "$job_hash->{num} with job_id $job_id failed, retry number $job_hash->{try_count}\n";
                    #unlink $job->{out} if -e $job->{out};
                    #unlink $job->{error} if -e $job->{error};
                    unlink "core" if -e "core";
                    my $new_job = $self->_setup_job($job_hash->{num});
                    if(!$new_job)
                    {   
                        $self->_kill_jobs($jobs);
                        last MONITOR;
                    }
                    $jobs->[$new_job->{num}] = $new_job;
                    $new_job->{try_count} = $job_hash->{try_count} + 1;
                    delete $running_jobs{ $job->id };
                    print "restarting $new_job->{num}\n";
                    $new_job->{job}->start;
                    print "restarted $new_job->{num} with new job_id ".$new_job->{job}->id."\n";
                    
                    #save job in running jobs list...
                    $running_jobs{ $new_job->{job}->id } = $new_job->{num};                 
                    
                }
                else
                {
                
                    print "$job_id failed, killing other jobs\n";
                    $self->_kill_jobs($jobs);
                    last MONITOR;
                }
            }
        }
    }

    return ( %running_jobs ) ? 0 : 1; # success is going thru all running jobs 
}

sub _kill_jobs
{
    my ($self, $jobs) = @_;

    for my $job_ref ( @$jobs )
    {
        my $job = $job_ref->{job};
        next if $job->has_ended;
        $job->kill;
    }
    
    return 1;
}

#run unchunk here
sub _finish
{
    my ($self, $success, $jobs) = @_;

    # Create the main reports, if jobs were successful
    if ( $success )
    {
        my $out_prefix = $self->snp_chunk_prefix($self->output);
 
        $self->snp_unchunk($out_prefix,$self->chunk_count,$self->output);
    }

    JOB: for my $job ( @$jobs )
    {
        LOG_TYPE: for my $log_type (qw/ out error /)
        {
            my $log_file_method = $log_type . '_log_file';
            my $log_file = $self->$log_file_method;
            next LOG_TYPE unless $log_file and -e $job->{$log_type};
            # cat the log file
            system sprintf('cat %s >> %s', $job->{$log_type}, $log_file);
            # remove the job's log file
            unlink $job->{$log_type};
        }
        #remove the chunked input snp file
        unlink $job->{snpchunk} if -e $job->{snpchunk}; 
        
        #remove the chunked output file        
        unlink $job->{outchunk} if -e $job->{outchunk};        
    }
    
    #remove chunked snps files directory
    rmdir $self->snpfile.'.chunk';
    #remove the chunk directory
    rmdir $self->output.'.chunk';

    return $success;
}

1;

=pod 

=head1 Name

Genome::Model::Command::Report::VariationsBatchToLsf

=head1 Synopsis

For a given variant file, this module will split it up into files containing 10000 variants.  Then for each of these split files, a job will be spawned and monitored that will run Genome::Model::Command::Report::Variation.  When the jobs are completed, the resulting report files will be combined into one report.  If requested the output and error logs from lsf can also be combined.

=head1 Usage

 $command = Genome::Model::Command::Report::VariationsBatchToLsf->execute
 (
     variant_type => 'snp', # opt, default snp, valid types: snp, indel
     variant_file => $detail_file, # req
     report_file => sprintf('%s/variant_report_for_chr_%s', $reports_dir, $chromosome), # req
     out_log_file => 'out', # opt, combine the out lsf log files here
     error_log_file => 'err', # opt, combine the error lsf log files here
     flank_range => 10000, # opt, default 50000
     variant_range => 0, # opt, default 0
 );

 if ( $command->result )
 { 
    ...
 }
 else 
 {
    ...
 }
 
=head1 Methods

=head2 execute or create then execute

=over

=item I<Synopsis>   Gets all annotations for a snp

=item I<Arguments>  snp (hash; see 'SNP' below)

=item I<Returns>    annotations (array of hash refs; see 'Annotation' below)

=back

=head1 See Also

B<Genome::SnpAnnotator>, B<Genome::Model::Command::Report::Variations>>

=head1 Disclaimer

Copyright (C) 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
