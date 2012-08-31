package Genome::Model::Tools::Snp::Filters::ScaledBinomialTest;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use Statistics::R;
use Workflow;

class Genome::Model::Tools::Snp::Filters::ScaledBinomialTest
{
    is => 'Command',
    has => [
        input_primary_snp_metrics => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'File of experimental metrics for the model',
        },
        input_other_snp_metrics =>
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'File of experimental metrics for the normal',
        },
        basedir => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Basename for the output file',
        },
        ratio => 
        {
            type => 'Float',
            is_optional => 1,
            default => .13,
            doc => 'Ratio to expect between model and normal reads due to contamination (.13 by default)',
        },
        ref_seq_id => 
        {
            type => 'String',
            is_optional => 0,
            is_input => 1,
            doc => 'Chromosome name or something',
        },
        two_sided_p_threshold => 
        {
            type => 'Float',
            is_optional => 1,
            default => .35,
            doc => 'Threshold below which to bin as non-skin in binomial test',
        },
        greater_p_threshold => 
        {
            type => 'Float',
            is_optional => 1,
            default => .35,
            doc => 'something dave larson thinks is good',
        },
        less_p_threshold => 
        {
            type => 'Float',
            is_optional => 1,
            default => .25,
            doc => 'Threshold below which to bin as non-skin in binomial test',
        },
        tmp_files_to_cleanup =>
        { 
            type => 'List',
            is_optional=> 1,
            doc => 'files to cleanup',
        },
        binomial_output_file =>
        {
            type => 'String',
            is_output => 1,
            doc => 'the name of the output file',
            calculate => q|
                        return $self->basedir . "/binomial.chr" . $self->ref_seq_id . ".nonskin.csv";
                       |
      
        },
    ],
};


#----------------------------------
sub execute {
    my $self = shift;
    $DB::single = $DB::stopper;

    my $in1     = $self->input_primary_snp_metrics;              # lie: the tumor "experimental" metrics
    my $in2     = $self->input_other_snp_metrics;   # lie: not experimental 
    my $out1    = $self->binomial_output_file;              # we make 2 files, the 2nd is derived

    if(-e $out1) {
        $self->warning_message("Binomial file exists. Not regenerating. You didn't see anything here.");
        return 1;
    }

    unless($self->_create_r_files) {
        $self->error_message("Error creating temporary files for R");
        return;
    }
    unless($self->_calculate_p_values) {
        $self->error_message("Error running R binomial test");
        return;
    }
    unless($self->_convert_probabilities_to_snps) {
        $self->error_message("Error splitting metrics file into skin and non-skin bins");
        return;
    }
    return 1;
}
#----------------------------------

sub _create_r_files {
    my $self = shift;
      
          
    #create temporary files for R
    my $tumor_metric_file = new IO::File($self->input_primary_snp_metrics);
    unless(defined($tumor_metric_file)) {
        $self->error_message("Couldn't open " . $self->input_primary_snp_metrics);
        return;
    }
    my $normal_metric_file = new IO::File($self->input_other_snp_metrics);
    unless(defined($normal_metric_file)) {
        $self->error_message("Couldn't open " . $self->input_other_snp_metrics);
        return;
    }
    
    my $chromosome = $self->ref_seq_id;

    my $ratio = $self->ratio;

    my $r_handle = new IO::File "> /tmp/skin_binom_test_chr$chromosome.csv";
    unless(defined($r_handle)) {
        $self->error_message("Couldn't open /tmp/skin_binom_test_chr$chromosome.csv for temporary R file");
        return;
    }

    #print header
    print $r_handle "skin_variant total_variant expected_proportion\n";

    my ($line,$skin_line);
    my $rows=0;
    
    while(($line = $tumor_metric_file->getline) && ($skin_line = $normal_metric_file->getline) && ++$rows) {
        #if (!defined($line) or !defined($skin_line)) {
        #    warn "premature end of data after row $rows!";
        #    last;
        #}
        $rows++;
        chomp $line;
        chomp $skin_line;
        if($line =~ /^chromosome/) {
            $line = $tumor_metric_file->getline;   #don't check additional headers that may have been included by cat
            chomp $line;
        }
        if($skin_line =~ /^chromosome/) {
            $skin_line = $normal_metric_file->getline;
            chomp $skin_line;
        }
        my @data_indices = (0, 1, 2, 3, 5, 24); 
        my @tumor_metrics = split /,\s*/, $line;
        my @skin_metrics = split /,\s*/, $skin_line;
       
        #unless(@tumor_metrics and @skin_metrics) {
        #    die "error parsing lines!\n$line\n$skin_line\n";
        #}


        my ($chr,
            $position,
            $al1,
            $al2,
            $al2_read_hg,
            $ref_read_hg,
        ) = @tumor_metrics[@data_indices];

        my ($skin_chr,
            $skin_position,
            $skin_al1,
            $skin_al2,
            $skin_al2_read_hg,
            $skin_ref_read_hg,
        ) = @skin_metrics[@data_indices];
        if($skin_chr ne $chr || $skin_position ne $position || $skin_al1 ne $al1 || $skin_al2 ne $al2) {
            #FILES ARE OFF
            #
            #Probably the skin site didn't make it through filtering
            if($skin_chr eq $chr && $skin_position > $position) {
                $self->error_message("Files de-synced during R data file creation");
                return;
            }
            else {
                $skin_line = $normal_metric_file->getline;
                redo;
            }
        }
        my $tumor_proportion = 1 - $ratio;
        my $tumor_coverage = $al2_read_hg + $ref_read_hg;
        my $skin_coverage = $skin_al2_read_hg + $skin_ref_read_hg;

        my $coverage_adjusted_proportion = $ratio * $skin_coverage / ($ratio * $skin_coverage + $tumor_proportion * $tumor_coverage);

        print $r_handle "$chr.$position.$al2 $skin_al2_read_hg ", $skin_al2_read_hg + $al2_read_hg, " $coverage_adjusted_proportion\n";  
    }
    return 1;
}

#----------------------------------
sub _calculate_p_values {
    my $self = shift;
    

    ############BEGIN INLINE R CODE##########
    my $R_code = <<'RCODE';

#R script to perform binomial test on a file and produce a new file with the p-values


#function ganked largely from qunyuan's gstat library
aml.binomial_test=function(x,alt="greater") {
    #we expect the x object to contain 3 columns
    #the first is the number of successes
    #the second is the number of trials
    #the third is the expected proportion
    p=numeric(0)
    for (i in c(1:nrow(x))) {
        xi=x[i,]
        pi=binom.test(as.numeric(xi[1]),as.numeric(xi[2]),as.numeric(xi[3]),alternative=alt)$p.value
        p=c(p,pi)
    }
    invisible(p)
}    

aml.epithelial_test=function(infile=NULL,outfile=NULL,alt="greater") {
    if(is.null(infile) || is.null(outfile)) {
        return(NULL)    
    }

    read.table(infile,header=T,row.names=1)->x
    aml.binomial_test(x,alt=alt)->p
    if(alt=="two.sided") {
        #auto convert the p value
        p = 1-p
    }
    if(alt=="greater") {
        #auto convert the p value
        p = 1-p
    }

    cbind(x,p)->x
    write.table(x,file=outfile,quote=F,sep="\t",row.names=T)
}

RCODE
    ############END INLINE R CODE##########

    my $chromosome = $self->ref_seq_id;
    my $snp_file_path = "/tmp/skin_binom_test_chr$chromosome.csv";


    #run the binomial test through R
    my $R_bridge = Statistics::R->new();
    unless(defined($R_bridge)) {
        $self->error_message("Couldn't access R via Statistics::R");
        return;
    }
    $R_bridge->startR();
    $R_bridge->send(qq{$R_code});
    #$R_bridge->send(qq{aml.epithelial_test(infile="$snp_file_path",outfile="/tmp/skin_binom_test_chr$chromosome.less",alt="less")});
    $R_bridge->send(qq{aml.epithelial_test(infile="$snp_file_path",outfile="/tmp/skin_binom_test_chr$chromosome.greater",alt="greater")});
    $R_bridge->stopR();
}

#-----------------------------------
sub _convert_probabilities_to_locations {
    my ($self,$file,$threshold,$hash_ref) = @_;
    if(!defined($hash_ref)) {
        return;
    }
    my $snp_filename = $self->input_primary_snp_metrics;

    my $handle = IO::File->new("$file");
    unless(defined($handle)) {
        $self->error_message("Couldn't open R t-test output file $file: $!");
        return;
    }

    my $header_line = $handle->getline; #ignore header
    while(my $line = $handle->getline) {
        chomp $line;
        my ($pos_str,
            $skin_reads,
            $total_reads,
            $expected_ratio,
            $p,
        ) = split "\t", $line;
        if($p <= $threshold || $skin_reads == 0) { 
            my ($chr,$pos,$var) = split /\./,  $pos_str;
            $hash_ref->{$chr}{$pos}{$var}=1;
        }
    }
    $handle->close;
    return 1;
}

    

#----------------------------------
sub _convert_probabilities_to_snps {
    my ($self) = @_;
    my $snp_filename = $self->input_primary_snp_metrics;
    my $chromosome = $self->ref_seq_id; 
    my $output = $self->basedir . "/binomial.chr" . $self->ref_seq_id;
    #expects both the two-sided and one-sided(less) tests to be present

    my %positions;
    my $result = $self->_convert_probabilities_to_locations("/tmp/skin_binom_test_chr$chromosome.greater", $self->greater_p_threshold,\%positions);    
    unless(defined($result)) {
        $self->error_message("Error creating locations for /tmp/skin_binom_test_chr$chromosome.greater");
        return;
    }
    #at this point, all the non-skin positions should be in %positions

    my $snp_file = new IO::File "< $snp_filename";
    unless(defined($snp_file)) {
        $self->error_message("Couldn't open metrics file $snp_filename for splitting");
        return;
    }

    my $non_skin_output_handle = new IO::File "> $output.nonskin.csv";
    unless(defined($non_skin_output_handle)) {
        $self->error_message("Couldn't open $output.nonskin.csv for writing");
        return;
    }
    
    my $skin_output_handle = new IO::File "> $output.skin.csv";
    unless(defined($skin_output_handle)) {
        $self->error_message("Couldn't open $output.skin.csv for writing");
        return;
    }

    print $non_skin_output_handle $snp_file->getline; #get header and print
    print $skin_output_handle $snp_file->getline; #get header and print

    while(my $line=$snp_file->getline) {
        chomp $line;
        my ($chr,
            $position,
            $al1,
            $al2,
        ) = split /,\s*/, $line;
        #nonskin
        if(exists($positions{$chr}) && exists($positions{$chr}{$position}) && exists($positions{$chr}{$position}{$al2})) {
            print $non_skin_output_handle $line,"\n";
        }
        else {
            print $skin_output_handle $line, "\n";
        }
    }
    return 1;
}

sub DESTROY {
    my $self = shift;
    #unlink all temporary files
    #TODO CHARRIS Add a class variable to track the files to unlink and unlink them here.   
} 
