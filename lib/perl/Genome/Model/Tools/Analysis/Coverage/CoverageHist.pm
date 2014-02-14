package Genome::Model::Tools::Analysis::Coverage::CoverageHist;

use strict;
use Genome;
use IO::File;
use warnings;


class Genome::Model::Tools::Analysis::Coverage::CoverageHist{
    is => 'Command',
    has => [
	refalign_model_id => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'refalign model id to grab coverage stats from',
	},

	somatic_validation_model_id => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'somatic-validation model id to grab coverage stats from',
	},

	somatic_variation_model_id => {
	    is => 'Number',
	    is_optional => 1,
	    doc => 'somatic-variation model id - gets associated refalign models and gets coverage stats from them',
	},

        output_pdf => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'output pdf containing histogram',
        },

        input_coverage_file => {
            is => 'String',
	    is_optional => 1,
	    doc => 'path to refcov file to use for input',
        },

        ]
};

sub help_brief {
    "Given either a model id or a coverage file, produce a histogram of coverage of the target region"
}

sub help_detail {
    "Given either a model id or a coverage file, produce a histogram of coverage of the target region. For paired tumor/normal somval builds, produces two plots."
}



sub execute {
    my $self = shift;
    my $output_pdf = $self->output_pdf;
    my $coverage_file = $self->input_coverage_file;

    my @covs;
    my @headers;

    #get the coverage file from the refalign model
    if(defined($self->refalign_model_id)){
        my $model = Genome::Model->get($self->refalign_model_id) or
            die "Could not find model ($self->refalign_model_id\n";
        my $build = $model->last_succeeded_build or
            die "Could not find last succeeded build from model $self->refalign_model_id.\n";
        my $dir = $build->data_directory;
        print STDERR $dir . "\n";
        my @a = glob("$dir/reference_coverage/wingspan_0/*_STATS.tsv");
        if(@a){
            push(@covs, $a[0]);
            push(@headers, "");
        }


    #get the coverage files from the somval model
    } elsif(defined($self->somatic_validation_model_id)){
        my $model = Genome::Model->get($self->somatic_validation_model_id) or
            die "Could not find model ($self->somatic_validation_model_id\n";
        my $build = $model->last_succeeded_build or
            die "Could not find last succeeded build from model $self->somatic_validation_model_id.\n";
        my $dir = $build->data_directory;
        print STDERR $dir . "\n";
        my @a = glob("$dir/coverage/normal/wingspan_0/*_STATS.tsv");
        if(@a){
            push(@covs, $a[0]);
            push(@headers, "Normal");
        }
        @a = glob("$dir/coverage/tumor/wingspan_0/*_STATS.tsv");
        if(@a){
            push(@covs, $a[0]);
            push(@headers, "Tumor");
        }

    #get the coverage files from the somvar model
    } elsif(defined($self->somatic_variation_model_id)){
        my $model = Genome::Model->get($self->somatic_variation_model_id) or
            die "Could not find somatic model";
        my $tumor_model = $model->tumor_model;
        my $normal_model = $model->normal_model;

        foreach my $m ($tumor_model,$normal_model){
            my $build = $m->last_succeeded_build or
                die "Could not find last succeeded build from model " . $m->id . "\n";
            my $dir = $build->data_directory;
            print STDERR $dir . "\n";
            my @a = glob("$dir/reference_coverage/wingspan_0/*_STATS.tsv");
            if(@a){
                push(@covs, $a[0]);
            }
        }
        push(@headers, "Normal");
        push(@headers, "Tumor");


    #just use the input coverage file
    } else {
        if(defined($coverage_file)){
            push(@covs, $coverage_file);
            push(@headers, "");
        }
    }

    if(@covs < 1){
        die "no coverage files found\n";
    }
    print STDERR "using coverage files:\n" . join("\n",@covs) . "\n";



    #create temporary r file
    my ($rfile,$newfile) = Genome::Sys->create_temp_file;
    unless($rfile) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }

    print $rfile "pdf(\"$output_pdf\", height=6, width=12)\n";

    for my $i (0..(@covs-1)){
        print $rfile 'par(mfrow=c(1,2))' . "\n";
        print $rfile 'b = read.table("' . $covs[$i] . "\")\n";
        print $rfile 'a = b[b$V13==1,]' . "\n";
        print $rfile 'hist(a$V6, breaks=seq(0,max(a$V6)+10,10),' . "col=\"darkgreen\",main=\"$headers[$i] Coverage of Target Regions\", ylab=\"Number of regions\", xlab=\"depth of coverage\")\n";
        print $rfile 'hist(a$V6, breaks=seq(0,max(a$V6)+10,10),' . "col=\"darkgreen\",main=\"$headers[$i] Coverage of Target Regions\", ylab=\"Number of regions\", xlab=\"depth of coverage\"," . 'xlim=c(0, sort(a$V6)[round(length(a$V6)*0.99)]))' . "\n";
        print $rfile "mtext(\"zoom to 99% of regions\")\n";
        print $rfile 'print(paste("Median Coverage:",median(a$V6)))' . "\n";
        print $rfile 'print(paste("Mean Coverage:",mean(a$V6)))' . "\n";

    }
    print $rfile "dev.off()\n";
    close($rfile);

    `cp $newfile /tmp/asdf.R`;

    my $cmd = "Rscript $newfile";
    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }
    return $return;


}
