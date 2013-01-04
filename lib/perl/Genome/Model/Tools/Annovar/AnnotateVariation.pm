package Genome::Model::Tools::Annovar::AnnotateVariation;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Annovar::AnnotateVariation {
    is => 'Genome::Model::Tools::Annovar::Base',
    has_input => [
        buildver => {
            is => "String",
        },
        table_names => {
            is => "String",
            is_many => 1,
        },
        input_file => {
            is => 'File',
        },
        outfile => {
            is => 'Path',
        },
        annotation_type => {
            is => 'String',
            valid_values => ["geneanno", "regionanno"],
        },
    ],
    has_optional_input => [
        scorecolumn => {
            is => 'Number',
        },
        bedfiles => {
            is => 'Text',
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    my $input_file = $self->input_file;
    if (_is_vcf($self->input_file)) {
        $self->status_message("Converting VCF to annovar");
        my $converter_input = $self->input_file;
        if (Genome::Sys->file_is_gzipped($converter_input)) {
            $self->status_message("Uncompressing VCF file");
            $converter_input = Genome::Sys->create_temp_file_path;
            my $in = Genome::Sys->open_gzip_file_for_reading($self->input_file);
            my $out = Genome::Sys->open_file_for_writing($converter_input);
            while(my $line = <$in>) {
                $out->print($line);
            }
            $in->close;
            $out->close;
        }
        $input_file = Genome::Sys->create_temp_file_path;
        my $convert = $self->script_path."convert2annovar.pl $converter_input -format vcf4 -outfile ".$input_file;
        $self->status_message("Running convert command: $convert");
        my $rv = Genome::Sys->shellcmd(cmd => $convert);
        unless ($rv) {
            $self->error_message("Could not convert vcf to annovar input");
            return $rv;
        }
    }

    foreach my $table_name ($self->table_names) {
        my $db;
        my %beds;
        if ($table_name eq "bed") {
            foreach my $bedfile ($self->bedfiles) {
                my $bedname;
                my $bedpath;
                ($bedname, $bedpath) = fileparse($bedfile);
                $beds{$bedname} = $bedpath;
            }
        }
        else {

            $db = Genome::Model::Tools::Annovar::Db->get_or_create(
                buildver => $self->buildver,
                table_name => $table_name,
                test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
            );

            unless ($db) {
                $self->error_message("Couldn't get ANNOVAR db result for ".$table_name." buildver ".$self->buildver);
                return;
            }

            $self->status_message("Annotating with ".$table_name." in ".$db->output_dir);
        }

        if ($table_name eq "bed") {
            foreach my $bedname (keys %beds) {
                my $bedpath = $beds{$bedname};
                my $cmd = $self->script_path."annotate_variation.pl --outfile ".$self->outfile." --buildver ".$self->buildver;
                if ($self->annotation_type eq "geneanno") {
                    $cmd .= " --geneanno ";
                }
                elsif ($self->annotation_type eq "regionanno") {
                    if (defined $self->scorecolumn) {
                        $cmd .= " -scorecolumn ".$self->scorecolumn." ";
                    }
                    $cmd .= " -regionanno -dbtype ".$table_name;
                    $cmd .= " -bedfile $bedname";
                }
                $cmd .= " ".$input_file;
                $cmd .= " $bedpath --colswanted 4";

                $self->status_message("Executing command $cmd");
                my $rv = Genome::Sys->shellcmd(cmd => $cmd);
                unless ($rv) {
                    $self->error_message("Failed to annotate with table $table_name $bedname");
                    return;
                }
                my $short_bed_name = $bedname;
                $short_bed_name =~ s/\.bed$//;
                my $mv_cmd = "mv ".$self->outfile.".".$self->buildver."_bed ".$self->outfile.".".$self->buildver."_bed_".$short_bed_name;
                print "Running mv cmd: $mv_cmd\n";
                `$mv_cmd`;
            }
        }
        else {
            my $cmd = $self->script_path."annotate_variation.pl --outfile ".$self->outfile." --buildver ".$self->buildver;
            if ($self->annotation_type eq "geneanno") {
                $cmd .= " --geneanno ";
            }
            elsif ($self->annotation_type eq "regionanno") {
                if (defined $self->scorecolumn) {
                    $cmd .= " -scorecolumn ".$self->scorecolumn." ";
                }
                $cmd .= " -regionanno -dbtype ".$table_name;
            }
            $cmd .= " ".$input_file;
            $cmd .= " ".$db->output_dir;

            $self->status_message("Executing command $cmd");
            my $rv = Genome::Sys->shellcmd(cmd => $cmd);
            unless ($rv) {
                $self->error_message("Failed to annotate with table $table_name");
                return;
            }
        }
    }

    return 1;
}

sub _is_vcf {
    my $in_file = shift;
    return ($in_file =~ /vcf$/ or $in_file =~ /vcf.gz$/);
}

1;

