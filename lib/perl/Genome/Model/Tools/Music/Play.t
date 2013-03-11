#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 13;

BEGIN {
    #It's important this be loaded first before the sub is redefined.
    use_ok('Genome::Model::Tools::Music::Play');
};

#don't run all the commands while testing the glue that binds them--instead just record that we tried
no warnings qw(redefine);
sub Genome::Model::Tools::Music::Play::_run_command {
    my $self = shift;
    my $command = shift;

    pass('Would have run ' . $command->command_name);
    return 1;
}
use warnings qw(redefine);


#fake parameters--since commands not being run (per above)
my $play_cmd = Genome::Model::Tools::Music::Play->create(
    bam_list => 'bam.list',
    roi_file => 'roi.file',
    reference_sequence => 'reference.sequence',
    output_dir => Genome::Sys->create_temp_directory,
    maf_file => 'maf.file',
    genetic_data_type => 'gene',
    pathway_file => 'pathway.file',
);

isa_ok($play_cmd, 'Genome::Model::Tools::Music::Play', 'created play comand');

my $rv = eval { $play_cmd->execute; };
if($@) {
    diag('Error executing: ' . $@);
}
ok($rv, 'Play executed successfully.');
