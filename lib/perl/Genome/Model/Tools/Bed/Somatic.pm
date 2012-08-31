package Genome::Model::Tools::Bed::Somatic;

use strict;
use warnings;

use Genome;
use File::Temp;

class Genome::Model::Tools::Bed::Somatic {
    is => ['Command'],
    has_input => [
        tumor_bed_file => {
            is => 'Text',
            is_input => 1,
            is_many => 1,
            doc => 'The BED format file of intervals to pull out bases from reference.',
        },
        normal_bed_file => {
            is => 'Text',
            is_input => 1,
            is_many => 1,
            doc => 'The BED format file of intervals to pull out bases from reference.',
        },
        somatic_file => {
            is => 'Text',
            is_output => 1,
            is_input => 1,
            doc => 'The stupid result',
        }
          ],
};


sub execute {
    my $self = shift;
    my @tumor_files = $self->tumor_bed_file;
    my @normal_files= $self->normal_bed_file;

    my $tumor_temp = Genome::Sys->create_temp_file_path;
    my $normal_temp = Genome::Sys->create_temp_file_path;
    my $tumor_out;
    my $normal_out;
    if(scalar(@tumor_files)>1){
        for my $file (@tumor_files){
            system("cat ".$file." >> ".$tumor_temp);
        }
        $tumor_out = $tumor_temp;
    }else{
        $tumor_out = $self->tumor_bed_file;
    }
    if(scalar(@normal_files)>1){
        for my $file (@normal_files){
            system("cat ".$file." >> ".$normal_temp);
        }
        $normal_out = $normal_temp;
    }else{
        $normal_out = $self->normal_bed_file;
    }

    
    my $tier1_cmd = $ENV{GENOME_SW} . "/bedtools/installed-64/intersectBed -wa -v -a " . $tumor_out . " -b " . $normal_out . " > " . $self->somatic_file;  

    my $result = system($tier1_cmd);
    
    return 1;
}
