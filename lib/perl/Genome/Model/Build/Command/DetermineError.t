use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::Command::DetermineError;

use Test::More tests => 5;
use File::Temp;

# constants defined at the bottom
sub PTERO_ERROR_LOG();
sub WORKFLOW_ERROR_LOG();
sub OUTPUT_LOG();

sub _write_to_temp_file {
    my $fh = File::Temp->new();
    $fh->print(@_);
    $fh->close;
    return $fh;
}

subtest 'parse ptero error log' => sub {
    plan tests => 5;

    my $error_log = _write_to_temp_file(PTERO_ERROR_LOG);

    my($error_source_file, $error_source_line, $error_host, $error_date, $error_text)
        = Genome::Model::Build::Command::DetermineError::parse_error_log($error_log->filename);

    is($error_source_file, '/path/to/ptero/error.pm', 'error source file');
    is($error_source_line, 1027, 'error source line');
    is($error_host, 'blade1-2-3.example.com', 'error host');
    is($error_date, '2014/12/04 12:25:42', 'error date');
    is($error_text, 'ERROR: This is a fake error', 'error text');
};

subtest 'workflow error log' => sub {
    plan tests => 5;

    my $error_log = _write_to_temp_file(WORKFLOW_ERROR_LOG);
    
    my($error_source_file, $error_source_line, $error_host, $error_date, $error_text)
        = Genome::Model::Build::Command::DetermineError::parse_error_log($error_log->filename);

    is($error_source_file, '/path/to/workflow/error.pm', 'error source file');
    is($error_source_line, 432, 'error source line');
    is($error_host, 'blade14-2-13', 'error host');
    is($error_date, '2014-12-12 20:28:14', 'error date');
    is($error_text, 'ERROR: This is a fake workflow error', 'error text');
};

subtest 'parse output log' => sub {
    plan tests => 5;

    my $output_log = _write_to_temp_file(OUTPUT_LOG);

    my($output_source_file, $output_source_line, $output_host, $output_date, $output_text)
        = Genome::Model::Build::Command::DetermineError::parse_output_log($output_log->filename);

    is($output_source_file, 'n/a', 'output source file');
    is($output_source_line, 'n/a', 'output source line');
    is($output_host, '4*blade3-3-3.example.com', 'output host');
    is($output_date, 'Tue Nov 25 11:16:46 2014', 'output date');
    is($output_text, 'TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.', 'output_text');
};

subtest 'find ptero die or warn in log' => sub {
    plan tests => 5;

    my $log = _write_to_temp_file(PTERO_ERROR_LOG);

    my($error_source_file, $error_source_line, $error_host, $error_date, $error_text)
        = Genome::Model::Build::Command::DetermineError::find_die_or_warn_in_log($log->filename);

    is($error_source_file, '/path/to/ptero/error.pm', 'error source file');
    is($error_source_line, 1027, 'error source line');
    is($error_host, 'blade1-2-3.example.com', 'error host');
    is($error_date, '2014/12/04 12:25:42', 'error date');
    is($error_text, 'ERROR: This is a fake error', 'error');
};

subtest 'find workflow die or warn in log' => sub {
    plan tests => 5;

    my $error_log = _write_to_temp_file(WORKFLOW_ERROR_LOG);
    my($error_source_file, $error_source_line, $error_host, $error_date, $error_text)
        = Genome::Model::Build::Command::DetermineError::find_die_or_warn_in_log($error_log->filename);

    is($error_source_file, '/path/to/the/workflow/exception.pm', 'error source file');
    is($error_source_line, 255, 'error source line');
    is($error_host, 'blade14-2-13', 'error host');
    is($error_date, '2014-12-12 20:28:21', 'error date');
    is($error_text, 'Workflow did not return correctly.', 'error text');
};

use constant PTERO_ERROR_LOG => <<'PTERO_ERROR';
2014-12-04 11:11:49,483 INFO flow.main naked_main 57: Loading command (workflow-wrapper)
2014-12-04 11:11:49,499 INFO flow_workflow.commands.workflow_wrapper _execute 80: Executing (blade4-3-2.example.com): workflow-wrapper.pl command shortcut Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference /tmp/tmpb8fz8b /tmp/tmp4K29zg
[2014/12/04 11:11:49.499618] Starting log annotation on host: blade4-3-2.example.com
[2014/12/04 11:11:51.583241] =========
[2014/12/04 11:11:51.583241] Attempting to shortcut command Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference...
[2014/12/04 11:11:51.583241] vvvvvvvvv
[2014/12/04 11:11:53.287182] ^^^^^^^^^
[2014/12/04 11:11:53.287182] Failed to shortcut command Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference...
[2014/12/04 11:11:53.287182] =========
2014-12-04 11:11:53,340 INFO flow_workflow.commands.workflow_wrapper _finish 102: Non-zero exit-code: 1 from perl_wrapper.
2014-12-04 11:11:53,341 INFO flow.util.exit exit_process 15: Exitting process: signalling children.
2014-12-04 11:11:53,732 INFO flow.util.exit exit_process 23: Children killed, exiting with code 1
2014-12-04 12:25:36,234 INFO flow.main naked_main 57: Loading command (lsf-pre-exec)
2014-12-04 12:25:36,235 INFO flow.shell_command.lsf.commands.pre_exec _execute 29: Begin LSF pre exec
2014-12-04 12:25:36,245 INFO flow.brokers.amqp.connection_manager _disconnect 147: Closing AMQP connection
2014-12-04 12:25:36,246 INFO flow.shell_command.lsf.commands.pre_exec _teardown 43: End LSF pre exec
2014-12-04 12:25:36,246 INFO flow.util.exit exit_process 15: Exitting process: signalling children.
2014-12-04 12:25:36,333 INFO flow.util.exit exit_process 23: Children killed, exiting with code 0
2014-12-04 12:25:37,258 INFO flow.main naked_main 57: Loading command (workflow-wrapper)
2014-12-04 12:25:37,265 INFO flow_workflow.commands.workflow_wrapper _execute 80: Executing (blade1-2-3.example.com): workflow-wrapper.pl command execute Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference /tmp/tmpVSPp6Y /tmp/tmppORDux
[2014/12/04 12:25:37.265717] Starting log annotation on host: blade1-2-3.example.com
[2014/12/04 12:25:38.593359] =========
[2014/12/04 12:25:38.593359] Attempting to execute command Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference...
[2014/12/04 12:25:38.593359] vvvvvvvvv
[2014/12/04 12:25:42.728920] Generated model name "H_KA-312451-1106998_SV_Contigs-human".
[2014/12/04 12:25:42.821850] ERROR: This is a fake error at /path/to/ptero/error.pm line 1027
[2014/12/04 12:25:42.821850]    UR::Context::create_entity('UR::Context::Process=HASH(0x2e86be8)', 'Genome::Model::ImportedReferenceSequence', 'UR::BoolExpr=HASH(0x7899e98)') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/UR/Object.pm line 20
[2014/12/04 12:25:42.821850]    UR::Object::create('Genome::Model::ImportedReferenceSequence', 'UR::BoolExpr=HASH(0x7899e98)') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Genome/Model.pm line 495
[2014/12/04 12:25:42.821850]    Genome::Model::create('Genome::Model::ImportedReferenceSequence', 'UR::BoolExpr=HASH(0x7899e98)') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Genome/ModelDeprecated.pm line 183
[2014/12/04 12:25:42.821850]    Genome::ModelDeprecated::create('Genome::Model::ImportedReferenceSequence', 'subject_type', 'species_name', 'subject_name', 'human', 'subject_class_name', 'Genome::Taxon', 'subject_id', 1653198737, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Genome/Model/Command/Define/ImportedReferenceSequence.pm line 305
[2014/12/04 12:25:42.821850]    Genome::Model::Command::Define::ImportedReferenceSequence::_get_or_create_model('Genome::Model::Command::Define::ImportedReferenceSequence=HAS...', 'Genome::Taxon=HASH(0x79897a8)') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Genome/Model/Command/Define/ImportedReferenceSequence.pm line 201
[2014/12/04 12:25:42.821850]    Genome::Model::Command::Define::ImportedReferenceSequence::execute('Genome::Model::Command::Define::ImportedReferenceSequence=HAS...') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Command/V2.pm line 214
[2014/12/04 12:25:42.821850]    Command::V2::execute('Genome::Model::Command::Define::ImportedReferenceSequence=HAS...') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Genome/Model/SomaticValidation/Command/ValidateSvs/CreateAssembledContigReference.pm line 116
[2014/12/04 12:25:42.821850]    Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference::execute('Genome::Model::SomaticValidation::Command::ValidateSvs::Creat...') called at /gsc/scripts/opt/genome/snapshots/genome-3546/lib/perl/Command/V2.pm line 214
[2014/12/04 12:25:42.821850]    Command::V2::execute('Genome::Model::SomaticValidation::Command::ValidateSvs::Creat...') called at /usr/bin/workflow-wrapper.pl line 199
[2014/12/04 12:25:42.821850]    eval {...} called at /usr/bin/workflow-wrapper.pl line 199
[2014/12/04 12:25:42.821850]    main::run_command('execute', 'Genome::Model::SomaticValidation::Command::ValidateSvs::Creat...', '/tmp/tmpVSPp6Y', '/tmp/tmppORDux') called at /usr/bin/workflow-wrapper.pl line 225
[2014/12/04 12:25:42.821850]    eval {...} called at /usr/bin/workflow-wrapper.pl line 225
[2014/12/04 12:25:42.821850]    main::safely_wrap('command', 'execute', 'Genome::Model::SomaticValidation::Command::ValidateSvs::Creat...', '/tmp/tmpVSPp6Y', '/tmp/tmppORDux') called at /usr/bin/workflow-wrapper.pl line 264
[2014/12/04 12:25:42.821850] 
[2014/12/04 12:25:42.821850] ^^^^^^^^^
[2014/12/04 12:25:42.821850] Crashed in execute for command Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference.
[2014/12/04 12:25:42.821850] =========
2014-12-04 12:25:42,865 INFO flow_workflow.commands.workflow_wrapper _finish 102: Non-zero exit-code: 1 from perl_wrapper.
2014-12-04 12:25:42,866 INFO flow.util.exit exit_process 15: Exitting process: signalling children.
2014-12-04 12:25:42,945 INFO flow.util.exit exit_process 23: Children killed, exiting with code 1
2014-12-04 12:25:43,251 INFO flow.main naked_main 57: Loading command (lsf-post-exec)
2014-12-04 12:25:43,252 INFO flow.shell_command.lsf.commands.post_exec _execute 33: Begin LSF post exec
2014-12-04 12:25:43,263 INFO flow.brokers.amqp.connection_manager _disconnect 147: Closing AMQP connection
2014-12-04 12:25:43,264 INFO flow.shell_command.lsf.commands.post_exec _teardown 81: End LSF post exec
2014-12-04 12:25:43,264 INFO flow.util.exit exit_process 15: Exitting process: signalling children.
2014-12-04 12:25:43,346 INFO flow.util.exit exit_process 23: Children killed, exiting with code 0
PTERO_ERROR

use constant OUTPUT_LOG => <<'OUTPUT';
[2014/11/23 04:42:13.243745] Starting log annotation on host: blade1-2-2.example.com
[2014/11/23 08:35:03.296285] Starting log annotation on host: blade3-3-3.example.com

------------------------------------------------------------
Sender: LSF System <lsfadmin@blade3-3-3.example.com>
Subject: Job 3327325: <RnaSeq Cufflinks Expression> Exited

Job <RnaSeq Cufflinks Expression> was submitted from host <blade1-1-1.example.com> by user <apipe-builder> in cluster <lsfcluster1>.
Job was executed on host(s) <4*blade3-3-3.example.com>, in queue <apipe>, as user <apipe-builder> in cluster <lsfcluster1>.
</gscuser/apipe-builder> was used as the home directory.
</> was used as the working directory.
Started at Sun Nov 23 08:35:01 2014
Results reported at Tue Nov 25 11:16:46 2014

Your job looked like:

------------------------------------------------------------
# LSBATCH: User input
'/usr/local/bin/flow' 'workflow-wrapper' '--method' 'execute' '--action-type' 'command' '--action-id' 'Genome::Model::RnaSeq::Command::Expression::Cufflinks' '--net-key' 'BOpgKBHUTB+wkgK6Gt5ljQ' '--operation-id' '12' '--parallel-id' '[[12, 2]]'
------------------------------------------------------------

TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.
Exited with signal termination: File size limit exceeded.

Resource usage summary:

    CPU time   : 579768.69 sec.
    Max Memory :     67902 MB
    Max Swap   :     75212 MB

    Max Processes  :         6
    Max Threads    :        11

The output (if any) is above this job summary.



PS:

Read file </gscmnt/gc13030/info/model_data/ab30047671a749e29455cd9d4da796c8/build3bf92295a6044c3a8d7fb9bfd7ddd973/logs/RnaSeq_Cufflinks_Expression.12.12_2.err> for stderr output of this job.
OUTPUT

use constant WORKFLOW_ERROR_LOG => <<'WORKFLOW_ERROR';
2014-12-12 20:23:06-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:06-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:07-0600 blade14-2-13: #################################################
2014-12-12 20:23:07-0600 blade14-2-13: Date Scheduled:  2014-12-12 20:23:06
2014-12-12 20:23:07-0600 blade14-2-13: Type:  execute
2014-12-12 20:23:07-0600 blade14-2-13: LSF Job Id:  4555537
2014-12-12 20:23:07-0600 blade14-2-13: Pid:  16096
2014-12-12 20:23:07-0600 blade14-2-13: HOST:  blade14-2-13.gsc.wustl.edu
2014-12-12 20:23:07-0600 blade14-2-13: USER:  ckang
2014-12-12 20:23:07-0600 blade14-2-13: DEBUG: Executing detect variants step
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: Created directory: /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: Created directory: /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/snv
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: Created directory: /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/indel
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: detector samtools supports SINGLE-sample detection with r963 []
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: detector varscan supports SINGLE-sample detection with 2.2.9 [--min-coverage 3 --min-var-freq 0.20 --p-value 0.10 --strand-filter 1 --map-quality 10]
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: detector varscan supports SINGLE-sample detection with 2.2.9 [--min-coverage 3 --min-var-freq 0.20 --p-value 0.10 --strand-filter 1 --map-quality 10]
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: Disconnecting Genome::DataSource::GMSchema default handle.
2014-12-12 20:23:10-0600 blade14-2-13: DEBUG: Now launching the dispatcher workflow.
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 hub server starting on port 54824
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 lsftail _start
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 lsftail Establishing Wheel on /usr/local/lsf/work/lsfcluster1/logdir/lsb.acct
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 dispatch _start
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 passthru _start
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 passthru Announcing we are at blade14-2-13.gsc.wustl.edu:54824 to /tmp/GG6shAECA6/hub_location
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 dispatch register Client16096
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 dispatch register alias 10.100.5.131:59895
2014-12-12 20:23:11-0600 blade14-2-13: 2014/12/12 20:23:11 passthru start ur (deferring until UR process connects)
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 ur process connecting to hub blade14-2-13.gsc.wustl.edu:54824
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 dispatch register blade14-2-13.gsc.wustl.edu-548ba31000003f17
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 dispatch register alias 10.100.5.131:59896
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 dispatch register alias UR
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 UR server _build hub_port=54824
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 workflow _start
2014-12-12 20:23:12-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:12-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 passthru Registered UR client.
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 pasthru pass_it_on
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 workflow load_and_execute (return_address:HASH(0x357e678))
2014-12-12 20:23:12-0600 blade14-2-13: 2014/12/12 20:23:12 Loading Workflow::Operation from xml
2014-12-12 20:23:13-0600 blade14-2-13: 2014/12/12 20:23:13 Loaded Workflow::Operation(blade14-2-13.gsc.wustl.edu 16151 1418437391 10022) from xml
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: 73135948 input connector
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: 73135953 indel varscan 2.2.9 #3
2014-12-12 20:23:14-0600 blade14-2-13: DEBUG: 73135950 snv samtools r963 #1
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch add_work 73135953
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch start_jobs 0 500
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Varscan /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135953.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135953.err
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Varscan; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch fork_worker started 16155
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch start_jobs submitted P16155 1
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch add_work 73135950
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 dispatch start_jobs 1 500
2014-12-12 20:23:14-0600 blade14-2-13: 2014/12/12 20:23:14 workflow schedule_instance
2014-12-12 20:23:15-0600 blade14-2-13: 2014/12/12 20:23:15 dispatch register blade14-2-13.gsc.wustl.edu-548ba31300003f1c
2014-12-12 20:23:15-0600 blade14-2-13: 2014/12/12 20:23:15 dispatch register alias 10.100.5.131:59904
2014-12-12 20:23:15-0600 blade14-2-13: 2014/12/12 20:23:15 dispatch register alias Worker
2014-12-12 20:23:15-0600 blade14-2-13: 2014/12/12 20:23:15 dispatch get_work P16155
2014-12-12 20:23:15-0600 blade14-2-13: 2014/12/12 20:23:15 workflow begin_instance
2014-12-12 20:23:15-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:15-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch end_work P16155 73135953
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch finalize_work 73135953
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba31300003f1c
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch unregister alias 10.100.5.131:59904
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 workflow end_instance
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch unregister alias Worker
2014-12-12 20:23:18-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch start_jobs 0 500
2014-12-12 20:23:18-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 workflow finalize_instance 73135953   
2014-12-12 20:23:18-0600 blade14-2-13: DEBUG: 73135954 indel_varscan_2.2.9_#3 false-indel v1 #2
2014-12-12 20:23:18-0600 blade14-2-13: DEBUG: 73135955 snv varscan 2.2.9 #3
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 Committing 5 (possible) changes.
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch sig_CHLD 16155 0
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch start_jobs 0 500
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Samtools /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135950.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135950.err
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Samtools; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch fork_worker started 16209
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch start_jobs submitted P16209 1
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 workflow schedule_instance
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch add_work 73135954
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch add_work 73135955
2014-12-12 20:23:18-0600 blade14-2-13: 2014/12/12 20:23:18 dispatch start_jobs 1 500
2014-12-12 20:23:19-0600 blade14-2-13: 2014/12/12 20:23:19 dispatch register blade14-2-13.gsc.wustl.edu-548ba31700003f52
2014-12-12 20:23:19-0600 blade14-2-13: 2014/12/12 20:23:19 dispatch register alias 10.100.5.131:59909
2014-12-12 20:23:19-0600 blade14-2-13: 2014/12/12 20:23:19 dispatch register alias Worker
2014-12-12 20:23:19-0600 blade14-2-13: 2014/12/12 20:23:19 dispatch get_work P16209
2014-12-12 20:23:19-0600 blade14-2-13: 2014/12/12 20:23:19 workflow begin_instance
2014-12-12 20:23:19-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:19-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch end_work P16209 73135950
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch finalize_work 73135950
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba31700003f52
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch unregister alias 10.100.5.131:59909
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch unregister alias Worker
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 workflow end_instance
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch start_jobs 0 500
2014-12-12 20:23:22-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:22-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 workflow finalize_instance 73135950   
2014-12-12 20:23:22-0600 blade14-2-13: DEBUG: 73135951 snv_samtools_r963_#1 snp-filter v1 #1
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 Committing 4 (possible) changes.
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch add_work 73135951
2014-12-12 20:23:22-0600 blade14-2-13: 2014/12/12 20:23:22 dispatch start_jobs 0 500
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch sig_CHLD 16209 0
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch start_jobs 0 500
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Filter::SnpFilter /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135951.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135951.err
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Filter::SnpFilter; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch fork_worker started 16235
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 dispatch start_jobs submitted P16235 1
2014-12-12 20:23:23-0600 blade14-2-13: 2014/12/12 20:23:23 workflow schedule_instance
2014-12-12 20:23:24-0600 blade14-2-13: 2014/12/12 20:23:24 dispatch register blade14-2-13.gsc.wustl.edu-548ba31c00003f6c
2014-12-12 20:23:24-0600 blade14-2-13: 2014/12/12 20:23:24 dispatch register alias 10.100.5.131:59914
2014-12-12 20:23:24-0600 blade14-2-13: 2014/12/12 20:23:24 dispatch register alias Worker
2014-12-12 20:23:24-0600 blade14-2-13: 2014/12/12 20:23:24 dispatch get_work P16235
2014-12-12 20:23:24-0600 blade14-2-13: 2014/12/12 20:23:24 workflow begin_instance
2014-12-12 20:23:24-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:24-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch end_work P16235 73135951
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch finalize_work 73135951
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba31c00003f6c
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch unregister alias 10.100.5.131:59914
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch unregister alias Worker
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 workflow end_instance
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch start_jobs 0 500
2014-12-12 20:23:28-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:28-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 workflow finalize_instance 73135951   
2014-12-12 20:23:28-0600 blade14-2-13: DEBUG: 73135952 snv_samtools_r963_#1 false-positive v1 #2
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 Committing 4 (possible) changes.
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch add_work 73135952
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch start_jobs 0 500
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch sig_CHLD 16235 0
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch start_jobs 0 500
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Filter::FalsePositive /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.err
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch fork_worker started 16261
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 dispatch start_jobs submitted P16261 1
2014-12-12 20:23:28-0600 blade14-2-13: 2014/12/12 20:23:28 workflow schedule_instance
2014-12-12 20:23:29-0600 blade14-2-13: 2014/12/12 20:23:29 dispatch register blade14-2-13.gsc.wustl.edu-548ba32100003f86
2014-12-12 20:23:29-0600 blade14-2-13: 2014/12/12 20:23:29 dispatch register alias 10.100.5.131:59919
2014-12-12 20:23:29-0600 blade14-2-13: 2014/12/12 20:23:29 dispatch register alias Worker
2014-12-12 20:23:29-0600 blade14-2-13: 2014/12/12 20:23:29 dispatch get_work P16261
2014-12-12 20:23:29-0600 blade14-2-13: 2014/12/12 20:23:29 workflow begin_instance
2014-12-12 20:23:29-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:29-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch end_work P16261 73135952
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch add_work 73135952
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch start_jobs 1 500
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 start_jobs calling command: bsub -R 'select[ncpus>=1 && mem>=8000] span[hosts=1] rusage[mem=8000]' -M 8000512 -n 1 -q apipe -o /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.out -e /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.err -g /apipe-workflow-worker -J "snv_samtools_r963_#1 false-positive v1 #2" -P "build14e3d3d9fd8047679e68429824e9137b" annotate-log /usr/bin/perl -e ' use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu", 54824)'
2014-12-12 20:23:32-0600 blade14-2-13: Workflow LSF Dispatcher issuing command: bsub -R 'select[ncpus>=1 && mem>=8000] span[hosts=1] rusage[mem=8000]' -M 8000512 -n 1 -q apipe -o /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.out -e /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135952.err -g /apipe-workflow-worker -J "snv_samtools_r963_#1 false-positive v1 #2" -P "build14e3d3d9fd8047679e68429824e9137b" annotate-log /usr/bin/perl -e ' use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu", 54824)' at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Dispatcher/Lsf.pm line 17.
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch start_jobs submitted 4555637 0
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba32100003f86
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch unregister alias 10.100.5.131:59919
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch unregister alias Worker
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 lsftail add_watcher 4555637
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch start_jobs 1 500
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch sig_CHLD 16261 0
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch start_jobs 1 500
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Filter::FalseIndel /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135954.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135954.err
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 workflow schedule_instance
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalseIndel; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch fork_worker started 16282
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 dispatch start_jobs submitted P16282 1
2014-12-12 20:23:32-0600 blade14-2-13: 2014/12/12 20:23:32 workflow schedule_instance
2014-12-12 20:23:34-0600 blade14-2-13: 2014/12/12 20:23:34 dispatch register blade14-2-13.gsc.wustl.edu-548ba32600003f9b
2014-12-12 20:23:34-0600 blade14-2-13: 2014/12/12 20:23:34 dispatch register alias 10.100.5.131:59922
2014-12-12 20:23:34-0600 blade14-2-13: 2014/12/12 20:23:34 dispatch register alias Worker
2014-12-12 20:23:34-0600 blade14-2-13: 2014/12/12 20:23:34 dispatch get_work P16282
2014-12-12 20:23:34-0600 blade14-2-13: 2014/12/12 20:23:34 workflow begin_instance
2014-12-12 20:23:34-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:34-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch end_work P16282 73135954
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch finalize_work 73135954
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba32600003f9b
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch unregister alias 10.100.5.131:59922
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch unregister alias Worker
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 workflow end_instance
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch start_jobs 1 500
2014-12-12 20:23:37-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:37-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 workflow finalize_instance 73135954   
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch sig_CHLD 16282 0
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch start_jobs 1 500
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Varscan /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135955.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135955.err
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Varscan; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch fork_worker started 16292
2014-12-12 20:23:37-0600 blade14-2-13: 2014/12/12 20:23:37 dispatch start_jobs submitted P16292 1
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 Committing 7 (possible) changes.
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 workflow schedule_instance
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 dispatch register blade14-2-13.gsc.wustl.edu-548ba32b00003fa5
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 dispatch register alias 10.100.5.131:59927
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 dispatch register alias Worker
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 dispatch get_work P16292
2014-12-12 20:23:39-0600 blade14-2-13: 2014/12/12 20:23:39 workflow begin_instance
2014-12-12 20:23:39-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:39-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch end_work P16292 73135955
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch finalize_work 73135955
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba32b00003fa5
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch unregister alias 10.100.5.131:59927
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch unregister alias Worker
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 workflow end_instance
2014-12-12 20:23:41-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch start_jobs 1 500
2014-12-12 20:23:41-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 workflow finalize_instance 73135955   
2014-12-12 20:23:41-0600 blade14-2-13: DEBUG: 73135956 snv_varscan_2.2.9_#3 false-positive v1 #2
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 Committing 4 (possible) changes.
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch add_work 73135956
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch start_jobs 1 500
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch sig_CHLD 16292 0
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch start_jobs 1 500
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch fork_worker Genome::Model::Tools::DetectVariants2::Filter::FalsePositive /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.out /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.err
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch fork_worker annotate-log /usr/bin/perl -e use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu" , 54824, 2)
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch fork_worker started 16312
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 dispatch start_jobs submitted P16312 1
2014-12-12 20:23:41-0600 blade14-2-13: 2014/12/12 20:23:41 workflow schedule_instance
2014-12-12 20:23:43-0600 blade14-2-13: 2014/12/12 20:23:43 dispatch register blade14-2-13.gsc.wustl.edu-548ba32f00003fb9
2014-12-12 20:23:43-0600 blade14-2-13: 2014/12/12 20:23:43 dispatch register alias 10.100.5.131:59932
2014-12-12 20:23:43-0600 blade14-2-13: 2014/12/12 20:23:43 dispatch register alias Worker
2014-12-12 20:23:43-0600 blade14-2-13: 2014/12/12 20:23:43 dispatch get_work P16312
2014-12-12 20:23:43-0600 blade14-2-13: 2014/12/12 20:23:43 workflow begin_instance
2014-12-12 20:23:43-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:23:43-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch end_work P16312 73135956
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch add_work 73135956
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch start_jobs 2 500
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 start_jobs calling command: bsub -R 'select[ncpus>=1 && mem>=8000] span[hosts=1] rusage[mem=8000]' -M 8000512 -n 1 -q apipe -o /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.out -e /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.err -g /apipe-workflow-worker -J "snv_varscan_2.2.9_#3 false-positive v1 #2" -P "build14e3d3d9fd8047679e68429824e9137b" annotate-log /usr/bin/perl -e ' use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu", 54824)'
2014-12-12 20:23:46-0600 blade14-2-13: Workflow LSF Dispatcher issuing command: bsub -R 'select[ncpus>=1 && mem>=8000] span[hosts=1] rusage[mem=8000]' -M 8000512 -n 1 -q apipe -o /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.out -e /gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488b9/build14e3d3d9fd8047679e68429824e9137b/variants/73135947/73135956.err -g /apipe-workflow-worker -J "snv_varscan_2.2.9_#3 false-positive v1 #2" -P "build14e3d3d9fd8047679e68429824e9137b" annotate-log /usr/bin/perl -e ' use Genome; use Genome::Model::Tools::DetectVariants2::Filter::FalsePositive; use Workflow::Server::Worker; Workflow::Server::Worker->start("blade14-2-13.gsc.wustl.edu", 54824)' at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Dispatcher/Lsf.pm line 17.
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch start_jobs submitted 4555678 0
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba32f00003fb9
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch unregister alias 10.100.5.131:59932
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch unregister alias Worker
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 lsftail add_watcher 4555678
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch start_jobs 2 500
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 workflow schedule_instance
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch sig_CHLD 16312 0
2014-12-12 20:23:46-0600 blade14-2-13: 2014/12/12 20:23:46 dispatch start_jobs 2 500
2014-12-12 20:24:11-0600 blade14-2-13: 2014/12/12 20:24:11 Committing 3 (possible) changes.
2014-12-12 20:26:42-0600 blade14-2-13: 2014/12/12 20:26:42 dispatch register blade14-4-14.gsc.wustl.edu-548ba3e20000653d
2014-12-12 20:26:42-0600 blade14-2-13: 2014/12/12 20:26:42 dispatch register alias 10.100.5.164:43935
2014-12-12 20:26:42-0600 blade14-2-13: 2014/12/12 20:26:42 dispatch register alias Worker
2014-12-12 20:26:42-0600 blade14-2-13: 2014/12/12 20:26:42 dispatch get_work 4555637
2014-12-12 20:26:42-0600 blade14-2-13: 2014/12/12 20:26:42 workflow begin_instance
2014-12-12 20:26:42-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:26:42-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 dispatch end_work 4555637 73135952
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 lsftail skip_watcher 4555637 60
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 dispatch unregister blade14-4-14.gsc.wustl.edu-548ba3e20000653d
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 dispatch unregister alias 10.100.5.164:43935
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 dispatch unregister alias Worker
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 workflow end_instance
2014-12-12 20:26:51-0600 blade14-2-13: 2014/12/12 20:26:51 dispatch start_jobs 1 500
2014-12-12 20:26:51-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:26:51-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:26:58-0600 blade14-2-13: 2014/12/12 20:26:58 dispatch register blade9-2-1.gsc.wustl.edu-548ba3f20000718d
2014-12-12 20:26:58-0600 blade14-2-13: 2014/12/12 20:26:58 dispatch register alias 10.100.4.87:50063
2014-12-12 20:26:58-0600 blade14-2-13: 2014/12/12 20:26:58 dispatch register alias Worker
2014-12-12 20:26:58-0600 blade14-2-13: 2014/12/12 20:26:58 dispatch get_work 4555678
2014-12-12 20:26:58-0600 blade14-2-13: 2014/12/12 20:26:58 workflow begin_instance
2014-12-12 20:26:58-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:26:58-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 dispatch end_work 4555678 73135956
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 lsftail skip_watcher 4555678 60
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 dispatch unregister blade9-2-1.gsc.wustl.edu-548ba3f20000718d
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 dispatch unregister alias 10.100.4.87:50063
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 dispatch unregister alias Worker
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 dispatch start_jobs 0 500
2014-12-12 20:27:09-0600 blade14-2-13: 2014/12/12 20:27:09 workflow end_instance
2014-12-12 20:27:09-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:27:09-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:27:11-0600 blade14-2-13: 2014/12/12 20:27:11 Committing 4 (possible) changes.
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 lsftail skip_it 4555637
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 lsftail delete_watcher 4555637
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 Tried to delete non-existent lsftail watcher (job_id=4555637)
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 dispatch finalize_work 73135952
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 workflow finalize_instance 73135952   
2014-12-12 20:27:51-0600 blade14-2-13: 2014/12/12 20:27:51 Committing 1 (possible) changes.
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 lsftail skip_it 4555678
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 lsftail delete_watcher 4555678
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 Tried to delete non-existent lsftail watcher (job_id=4555678)
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 dispatch finalize_work 73135956
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 workflow finalize_instance 73135956   
2014-12-12 20:28:09-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:28:09-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 Committing 1 (possible) changes.
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 workflow error_relay
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 passthru Relaying result from ur process to Client16096
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 dispatch unregister Client16096
2014-12-12 20:28:09-0600 blade14-2-13: 2014/12/12 20:28:09 dispatch unregister alias 10.100.5.131:59895
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch unregister blade14-2-13.gsc.wustl.edu-548ba31000003f17
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch unregister alias 10.100.5.131:59896
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch unregister alias UR
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch sig_HUP_INT_TERM
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch close_out
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch close_out clear queue
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch close_out bkill pending
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 dispatch close_out bkill running
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 Genome::Model::Tools::DetectVariants2::Dispatcher id(blade14-2-13.gsc.wustl.edu 16096 1418437384 10043): DetectVariants2 Dispatcher: Execution halted due to unresolvable dependencies or crashed children.  Status and incomplete inputs:
2014-12-12 20:28:14-0600 blade14-2-13: input connector <73135948> (done)
2014-12-12 20:28:14-0600 blade14-2-13: output connector <73135949> (new)
2014-12-12 20:28:14-0600 blade14-2-13:   snv_output_directory
2014-12-12 20:28:14-0600 blade14-2-13:   snv_result_id
2014-12-12 20:28:14-0600 blade14-2-13:   snv_result_class
2014-12-12 20:28:14-0600 blade14-2-13: snv samtools r963 #1 <73135950> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 snp-filter v1 #1 <73135951> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 false-positive v1 #2 <73135952> (crashed)
2014-12-12 20:28:14-0600 blade14-2-13: indel varscan 2.2.9 #3 <73135953> (done)
2014-12-12 20:28:14-0600 blade14-2-13: indel_varscan_2.2.9_#3 false-indel v1 #2 <73135954> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv varscan 2.2.9 #3 <73135955> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_varscan_2.2.9_#3 false-positive v1 #2 <73135956> (crashed)
2014-12-12 20:28:14-0600 blade14-2-13: unionunique snv_samtools_r963_#1 snv_varscan_2.2.9_#3 <73135957> (new)
2014-12-12 20:28:14-0600 blade14-2-13:   input_a_id
2014-12-12 20:28:14-0600 blade14-2-13:   input_b_id
2014-12-12 20:28:14-0600 blade14-2-13: 
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 false-positive v1 #2: Found allocation at (/gscmnt/gc13026/info/build_merged_alignments/detect-variants--blade9-1-3.gsc.wustl.edu-ckang-27644-b7b71bbdc4b54ad78c365b02567e9f42) but no software result for it's owner ID (b7b71bbdc4b54ad78c365b02567e9f42). This is either because the software result is currently being generated or because the allocation has been orphaned. If it is determined that the allocation has been orphaned then the all
2014-12-12 20:28:14-0600 blade14-2-13: ERROR: This is a fake workflow error at /path/to/workflow/error.pm line 432.
2014-12-12 20:28:14-0600 blade14-2-13: input connector <73135948> (done)
2014-12-12 20:28:14-0600 blade14-2-13: output connector <73135949> (new)
2014-12-12 20:28:14-0600 blade14-2-13:   snv_output_directory
2014-12-12 20:28:14-0600 blade14-2-13:   snv_result_id
2014-12-12 20:28:14-0600 blade14-2-13:   snv_result_class
2014-12-12 20:28:14-0600 blade14-2-13: snv samtools r963 #1 <73135950> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 snp-filter v1 #1 <73135951> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 false-positive v1 #2 <73135952> (crashed)
2014-12-12 20:28:14-0600 blade14-2-13: indel varscan 2.2.9 #3 <73135953> (done)
2014-12-12 20:28:14-0600 blade14-2-13: indel_varscan_2.2.9_#3 false-indel v1 #2 <73135954> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv varscan 2.2.9 #3 <73135955> (done)
2014-12-12 20:28:14-0600 blade14-2-13: snv_varscan_2.2.9_#3 false-positive v1 #2 <73135956> (crashed)
2014-12-12 20:28:14-0600 blade14-2-13: unionunique snv_samtools_r963_#1 snv_varscan_2.2.9_#3 <73135957> (new)
2014-12-12 20:28:14-0600 blade14-2-13:   input_a_id
2014-12-12 20:28:14-0600 blade14-2-13:   input_b_id
2014-12-12 20:28:14-0600 blade14-2-13: 
2014-12-12 20:28:14-0600 blade14-2-13: snv_samtools_r963_#1 false-positive v1 #2: Found allocation at (/gscmnt/gc13026/info/build_merged_alignments/detect-variants--blade9-1-3.gsc.wustl.edu-ckang-27644-b7b71bbdc4b54ad78c365b02567e9f42) but no software result for it's owner ID (b7b71bbdc4b54ad78c365b02567e9f42). This is either because the software result is currently being generated or because the allocation has been orphaned. If it is determined that the allocation has been orphaned then the allocation will need to be removed. at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 190.
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::_check_instance_output('Genome::Model::Tools::DetectVariants2::Result::Filter=HASH(0x...', '/gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 128
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 123
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', 'Genome::Model::Tools::DetectVariants2::Filter::SnpFilter v1', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/SoftwareResult.pm line 252
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::SoftwareResult::get_or_create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', 'Genome::Model::Tools::DetectVariants2::Filter::SnpFilter v1', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 344
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::_summon_filter_result('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 290
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 398
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::call('Workflow::OperationType::Command=HASH(0x3aee850)', 'execute', 'params', 'Workflow::Link::Instance=HASH(0x3aee1d8)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x3aee370)', 'version', 'Workflow::Link::Instance=HASH(0x3aee418)', 'pedigree_file_path', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 284
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::execute('Workflow::OperationType::Command=HASH(0x3aee850)', 'params', 'Workflow::Link::Instance=HASH(0x3aee1d8)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x3aee370)', 'version', 'Workflow::Link::Instance=HASH(0x3aee418)', 'pedigree_file_path', 'Workflow::Link::Instance=HASH(0x3aee4c0)', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x39070a8)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3907000)', 'execute', 'POE::Session=ARRAY(0x3b1b178)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x39070a8)', 'POE::Session=ARRAY(0x3b1b178)', 'execute', 'ARRAY(0x3b24df8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x39070a8)', 'POE::Session=ARRAY(0x39070a8)', 'execute', 2, 'ARRAY(0x3b24df8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x39070a8)', 'execute', 'ARRAY(0x3ae3838)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x3b1b178)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3b1b0e8)', '__thunk', 'POE::Session=ARRAY(0x3aaf338)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3b1b178)', 'POE::Session=ARRAY(0x3aaf338)', '__thunk', 'ARRAY(0x3b24a38)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3b1b178)', 'POE::Session=ARRAY(0x39070a8)', '__thunk', 2, 'ARRAY(0x3b24a38)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3b1b178)', '__thunk', undef, 'ARRAY(0x3b22110)', 'HASH(0x3b21e70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x3b22110)', 'HASH(0x3b21e70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x3a919c0)', 'HASH(0x3aea1f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x3aaf338)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3aaf290)', 'request', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3aaf338)', 'POE::Session=ARRAY(0x3ae5b30)', 'request', 'ARRAY(0x3aec0b0)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3aaf338)', 'POE::Session=ARRAY(0x39070a8)', 'request', 2, 'ARRAY(0x3aec0b0)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'IKC', 'request', 'HASH(0x3aea1f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3ae5aa0)', 'receive', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'receive', 'ARRAY(0x35d4978)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3ae5b30)', 'receive', 'HASH(0x3aea1f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3ae5aa0)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x3ae3b50)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x3ae3b50)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x356b7a8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x356b7a8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade14-2-13.gsc.wustl.edu', 54824) called at -e line 1
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', 'Genome::Model::Tools::DetectVariants2::Filter::SnpFilter v1', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/SoftwareResult.pm line 252
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::SoftwareResult::get_or_create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', 'Genome::Model::Tools::DetectVariants2::Filter::SnpFilter v1', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 344
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::_summon_filter_result('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 290
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 398
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::call('Workflow::OperationType::Command=HASH(0x3aee850)', 'execute', 'params', 'Workflow::Link::Instance=HASH(0x3aee1d8)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x3aee370)', 'version', 'Workflow::Link::Instance=HASH(0x3aee418)', 'pedigree_file_path', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 284
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::execute('Workflow::OperationType::Command=HASH(0x3aee850)', 'params', 'Workflow::Link::Instance=HASH(0x3aee1d8)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x3aee370)', 'version', 'Workflow::Link::Instance=HASH(0x3aee418)', 'pedigree_file_path', 'Workflow::Link::Instance=HASH(0x3aee4c0)', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x39070a8)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3907000)', 'execute', 'POE::Session=ARRAY(0x3b1b178)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x39070a8)', 'POE::Session=ARRAY(0x3b1b178)', 'execute', 'ARRAY(0x3b24df8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x39070a8)', 'POE::Session=ARRAY(0x39070a8)', 'execute', 2, 'ARRAY(0x3b24df8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x39070a8)', 'execute', 'ARRAY(0x3ae3838)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x3b1b178)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3b1b0e8)', '__thunk', 'POE::Session=ARRAY(0x3aaf338)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3b1b178)', 'POE::Session=ARRAY(0x3aaf338)', '__thunk', 'ARRAY(0x3b24a38)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3b1b178)', 'POE::Session=ARRAY(0x39070a8)', '__thunk', 2, 'ARRAY(0x3b24a38)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3b1b178)', '__thunk', undef, 'ARRAY(0x3b22110)', 'HASH(0x3b21e70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x3b22110)', 'HASH(0x3b21e70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x3a919c0)', 'HASH(0x3aea1f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x3aaf338)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3aaf290)', 'request', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3aaf338)', 'POE::Session=ARRAY(0x3ae5b30)', 'request', 'ARRAY(0x3aec0b0)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3aaf338)', 'POE::Session=ARRAY(0x39070a8)', 'request', 2, 'ARRAY(0x3aec0b0)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'IKC', 'request', 'HASH(0x3aea1f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3ae5aa0)', 'receive', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'receive', 'ARRAY(0x35d4978)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3ae5b30)', 'receive', 'HASH(0x3aea1f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Kernel=ARRAY(0x356b7a8)', 'HASH(0x3ae5aa0)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x3ae5b30)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x3ae3b50)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x356b7a8)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Session=ARRAY(0x3ae5b30)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x3ae3b50)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x356b7a8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x356b7a8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade14-2-13.gsc.wustl.edu', 54824) called at -e line 1
2014-12-12 20:28:14-0600 blade14-2-13: 
2014-12-12 20:28:14-0600 blade14-2-13: snv_varscan_2.2.9_#3 false-positive v1 #2: Found allocation at (/gscmnt/gc13027/info/build_merged_alignments/detect-variants--blade11-4-2.gsc.wustl.edu-ckang-15205-e4d0dcabcb0b40b98ac47cd90b179323) but no software result for it's owner ID (e4d0dcabcb0b40b98ac47cd90b179323). This is either because the software result is currently being generated or because the allocation has been orphaned. If it is determined that the allocation has been orphaned then the allocation will need to be removed. at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 190.
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::_check_instance_output('Genome::Model::Tools::DetectVariants2::Result::Filter=HASH(0x...', '/gscmnt/gc9020/info/model_data/078b991a66ee4db2a256ac94fa1488...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 128
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Result/DetectionBase.pm line 123
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/SoftwareResult.pm line 252
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::SoftwareResult::get_or_create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 344
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::_summon_filter_result('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 290
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 398
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::call('Workflow::OperationType::Command=HASH(0x38dd910)', 'execute', 'params', 'Workflow::Link::Instance=HASH(0x38dd298)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x38dd430)', 'version', 'Workflow::Link::Instance=HASH(0x38dd4d8)', 'pedigree_file_path', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 284
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::execute('Workflow::OperationType::Command=HASH(0x38dd910)', 'params', 'Workflow::Link::Instance=HASH(0x38dd298)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x38dd430)', 'version', 'Workflow::Link::Instance=HASH(0x38dd4d8)', 'pedigree_file_path', 'Workflow::Link::Instance=HASH(0x38dd580)', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x36f6178)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x36f60d0)', 'execute', 'POE::Session=ARRAY(0x390a0a8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x36f6178)', 'POE::Session=ARRAY(0x390a0a8)', 'execute', 'ARRAY(0x3913cc8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x36f6178)', 'POE::Session=ARRAY(0x36f6178)', 'execute', 2, 'ARRAY(0x3913cc8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x36f6178)', 'execute', 'ARRAY(0x38d28f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x390a0a8)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x390a018)', '__thunk', 'POE::Session=ARRAY(0x389e468)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x390a0a8)', 'POE::Session=ARRAY(0x389e468)', '__thunk', 'ARRAY(0x3913908)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x390a0a8)', 'POE::Session=ARRAY(0x36f6178)', '__thunk', 2, 'ARRAY(0x3913908)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x390a0a8)', '__thunk', undef, 'ARRAY(0x3911010)', 'HASH(0x3910d70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x3911010)', 'HASH(0x3910d70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x3880a90)', 'HASH(0x38d92b8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x389e468)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x389e3c0)', 'request', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x389e468)', 'POE::Session=ARRAY(0x38d4bf0)', 'request', 'ARRAY(0x38db170)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x389e468)', 'POE::Session=ARRAY(0x36f6178)', 'request', 2, 'ARRAY(0x38db170)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'IKC', 'request', 'HASH(0x38d92b8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x38d4b60)', 'receive', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'receive', 'ARRAY(0x33c3a68)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x38d4bf0)', 'receive', 'HASH(0x38d92b8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x38d4b60)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x38d2c10)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x38d2c10)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x335a778)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x335a778)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade14-2-13.gsc.wustl.edu', 54824) called at -e line 1
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Result::DetectionBase::create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/SoftwareResult.pm line 252
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::SoftwareResult::get_or_create('Genome::Model::Tools::DetectVariants2::Result::Filter', 'chromosome_list', undef, 'reference_build_id', 106942997, 'filter_name', 'Genome::Model::Tools::DetectVariants2::Filter::FalsePositive', 'previous_filter_strategy', undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 344
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::_summon_filter_result('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Filter.pm line 290
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Filter::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Filter::FalsePositive=...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 398
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::call('Workflow::OperationType::Command=HASH(0x38dd910)', 'execute', 'params', 'Workflow::Link::Instance=HASH(0x38dd298)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x38dd430)', 'version', 'Workflow::Link::Instance=HASH(0x38dd4d8)', 'pedigree_file_path', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Command.pm line 284
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Command::execute('Workflow::OperationType::Command=HASH(0x38dd910)', 'params', 'Workflow::Link::Instance=HASH(0x38dd298)', 'previous_result_id', 'Workflow::Link::Instance=HASH(0x38dd430)', 'version', 'Workflow::Link::Instance=HASH(0x38dd4d8)', 'pedigree_file_path', 'Workflow::Link::Instance=HASH(0x38dd580)', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x36f6178)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x36f60d0)', 'execute', 'POE::Session=ARRAY(0x390a0a8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x36f6178)', 'POE::Session=ARRAY(0x390a0a8)', 'execute', 'ARRAY(0x3913cc8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x36f6178)', 'POE::Session=ARRAY(0x36f6178)', 'execute', 2, 'ARRAY(0x3913cc8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x36f6178)', 'execute', 'ARRAY(0x38d28f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x390a0a8)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x390a018)', '__thunk', 'POE::Session=ARRAY(0x389e468)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x390a0a8)', 'POE::Session=ARRAY(0x389e468)', '__thunk', 'ARRAY(0x3913908)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x390a0a8)', 'POE::Session=ARRAY(0x36f6178)', '__thunk', 2, 'ARRAY(0x3913908)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x390a0a8)', '__thunk', undef, 'ARRAY(0x3911010)', 'HASH(0x3910d70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x3911010)', 'HASH(0x3910d70)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x3880a90)', 'HASH(0x38d92b8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x389e468)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x389e3c0)', 'request', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x389e468)', 'POE::Session=ARRAY(0x38d4bf0)', 'request', 'ARRAY(0x38db170)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x389e468)', 'POE::Session=ARRAY(0x36f6178)', 'request', 2, 'ARRAY(0x38db170)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'IKC', 'request', 'HASH(0x38d92b8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x38d4b60)', 'receive', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'receive', 'ARRAY(0x33c3a68)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x38d4bf0)', 'receive', 'HASH(0x38d92b8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Kernel=ARRAY(0x335a778)', 'HASH(0x38d4b60)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x38d4bf0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x38d2c10)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x335a778)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Session=ARRAY(0x38d4bf0)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x38d2c10)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x335a778)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x335a778)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade14-2-13.gsc.wustl.edu', 54824) called at -e line 1
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 Genome::Model::Tools::DetectVariants2::Dispatcher id(blade14-2-13.gsc.wustl.edu 16096 1418437384 10043): Workflow did not return correctly.
2014-12-12 20:28:14-0600 blade14-2-13: ERROR: Workflow did not return correctly.
2014-12-12 20:28:14-0600 blade14-2-13: 2014/12/12 20:28:14 Workflow::OperationType::Event id(6f3ffd815f2342cb86ee887dda0d72d3): Workflow did not return correctly. at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Dispatcher.pm line 255.
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Dispatcher::_detect_variants('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Base.pm line 125
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Base::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Event/Build/ReferenceAlignment/DetectVariants.pm line 54
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Event::Build::ReferenceAlignment::DetectVariants::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V1::execute('Genome::Model::Ev
2014-12-12 20:28:14-0600 blade14-2-13: ERROR: Workflow did not return correctly. at /path/to/other/workflow/error.pm line 255.
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Dispatcher::_detect_variants('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Base.pm line 125
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Base::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Event/Build/ReferenceAlignment/DetectVariants.pm line 54
2014-12-12 20:28:14-0600 blade14-2-13: 	Genome::Model::Event::Build::ReferenceAlignment::DetectVariants::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V1::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 201
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 198
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Event::call('Workflow::OperationType::Event=HASH(0x3b23118)', 'execute', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 106
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Event::execute('Workflow::OperationType::Event=HASH(0x3b23118)', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x552c420)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552c378)', 'execute', 'POE::Session=ARRAY(0x555e0c0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x555e0c0)', 'execute', 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x552c420)', 'execute', 2, 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'execute', 'ARRAY(0x552c6c0)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x555e030)', '__thunk', 'POE::Session=ARRAY(0x54fb578)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x54fb578)', '__thunk', 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x552c420)', '__thunk', 2, 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', '__thunk', undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x5219418)', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x54fb578)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x54fb4d0)', 'request', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552f4c8)', 'request', 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552c420)', 'request', 2, 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'IKC', 'request', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'receive', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'ARRAY(0x5020630)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'HASH(0x55342f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade9-3-10.gsc.wustl.edu', 51228) called at -e line 1
2014-12-12 20:28:14-0600 blade14-2-13: 	Command::V1::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 201
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 198
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Event::call('Workflow::OperationType::Event=HASH(0x3b23118)', 'execute', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 106
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::OperationType::Event::execute('Workflow::OperationType::Event=HASH(0x3b23118)', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x552c420)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552c378)', 'execute', 'POE::Session=ARRAY(0x555e0c0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x555e0c0)', 'execute', 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x552c420)', 'execute', 2, 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'execute', 'ARRAY(0x552c6c0)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x555e030)', '__thunk', 'POE::Session=ARRAY(0x54fb578)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x54fb578)', '__thunk', 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x552c420)', '__thunk', 2, 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', '__thunk', undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x5219418)', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x54fb578)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x54fb4d0)', 'request', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552f4c8)', 'request', 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552c420)', 'request', 2, 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'IKC', 'request', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'receive', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'ARRAY(0x5020630)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'HASH(0x55342f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:14-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:14-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:14-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade9-3-10.gsc.wustl.edu', 51228) called at -e line 1
2014-12-12 20:28:21-0600 blade14-2-13: DEBUG: config key datetime exists
2014-12-12 20:28:21-0600 blade14-2-13: DEBUG: got datetime format %Y-%m-%d %H:%M:%S
2014-12-12 20:28:21-0600 blade14-2-13: Command module died or returned undef.
2014-12-12 20:28:21-0600 blade14-2-13: Workflow did not return correctly. at /path/to/the/workflow/exception.pm line 255.
2014-12-12 20:28:21-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Dispatcher::_detect_variants('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Tools/DetectVariants2/Base.pm line 125
2014-12-12 20:28:21-0600 blade14-2-13: 	Genome::Model::Tools::DetectVariants2::Base::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V2.pm line 214
2014-12-12 20:28:21-0600 blade14-2-13: 	Command::V2::execute('Genome::Model::Tools::DetectVariants2::Dispatcher=HASH(0x8829...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Genome/Model/Event/Build/ReferenceAlignment/DetectVariants.pm line 54
2014-12-12 20:28:21-0600 blade14-2-13: 	Genome::Model::Event::Build::ReferenceAlignment::DetectVariants::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Command/V1.pm line 135
2014-12-12 20:28:21-0600 blade14-2-13: 	Command::V1::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 201
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 198
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::call('Workflow::OperationType::Event=HASH(0x3b23118)', 'execute', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 106
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::execute('Workflow::OperationType::Event=HASH(0x3b23118)', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x552c420)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552c378)', 'execute', 'POE::Session=ARRAY(0x555e0c0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x555e0c0)', 'execute', 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x552c420)', 'execute', 2, 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'execute', 'ARRAY(0x552c6c0)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x555e030)', '__thunk', 'POE::Session=ARRAY(0x54fb578)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x54fb578)', '__thunk', 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x552c420)', '__thunk', 2, 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', '__thunk', undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x5219418)', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x54fb578)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x54fb4d0)', 'request', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552f4c8)', 'request', 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552c420)', 'request', 2, 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'IKC', 'request', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'receive', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'ARRAY(0x5020630)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'HASH(0x55342f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade9-3-10.gsc.wustl.edu', 51228) called at -e line 1
2014-12-12 20:28:21-0600 blade14-2-13: 	Command::V1::execute('Genome::Model::Event::Build::ReferenceAlignment::DetectVarian...') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 201
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 198
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::call('Workflow::OperationType::Event=HASH(0x3b23118)', 'execute', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 106
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::execute('Workflow::OperationType::Event=HASH(0x3b23118)', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x552c420)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552c378)', 'execute', 'POE::Session=ARRAY(0x555e0c0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x555e0c0)', 'execute', 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x552c420)', 'execute', 2, 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'execute', 'ARRAY(0x552c6c0)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x555e030)', '__thunk', 'POE::Session=ARRAY(0x54fb578)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x54fb578)', '__thunk', 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x552c420)', '__thunk', 2, 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', '__thunk', undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x5219418)', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x54fb578)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x54fb4d0)', 'request', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552f4c8)', 'request', 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552c420)', 'request', 2, 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'IKC', 'request', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'receive', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'ARRAY(0x5020630)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'HASH(0x55342f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade9-3-10.gsc.wustl.edu', 51228) called at -e line 1
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::call('Workflow::OperationType::Event=HASH(0x3b23118)', 'execute', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/OperationType/Event.pm line 106
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::OperationType::Event::execute('Workflow::OperationType::Event=HASH(0x3b23118)', 'prior_result', 'Workflow::Link::Instance=HASH(0x3b22e18)', 'prior_result', 1) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 132
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 122
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::_worker_execute(undef, 'POE::Session=ARRAY(0x552c420)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552c378)', 'execute', 'POE::Session=ARRAY(0x555e0c0)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x555e0c0)', 'execute', 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, '__thunk') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'POE::Session=ARRAY(0x552c420)', 'execute', 2, 'ARRAY(0x5561490)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1510, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552c420)', 'execute', 'ARRAY(0x552c6c0)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1510
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::__thunk('POE::Component::IKC::Responder::Thunk', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x555e030)', '__thunk', 'POE::Session=ARRAY(0x54fb578)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x54fb578)', '__thunk', 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'request') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', 'POE::Session=ARRAY(0x552c420)', '__thunk', 2, 'ARRAY(0x555e798)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 1413, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x555e0c0)', '__thunk', undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 1413
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Thunk::thunk(undef, 'ARRAY(0x555dd90)', 'HASH(0x555daf0)', undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 472
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 424
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::Object::request('POE::Component::IKC::Responder::Object=HASH(0x5219418)', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Responder.pm line 108
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Responder::request('POE::Component::IKC::Responder', 'POE::Session=ARRAY(0x54fb578)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x54fb4d0)', 'request', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 482
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552f4c8)', 'request', 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'POE::Wheel::ReadWrite(2) -> select read') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1092
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1078
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x54fb578)', 'POE::Session=ARRAY(0x552c420)', 'request', 2, 'ARRAY(0x55345c8)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Co...', 632, 'execute', ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1789
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'IKC', 'request', 'HASH(0x55342f8)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Component/IKC/Channel.pm line 632
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Component::IKC::Channel::channel_receive(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'receive', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'ARRAY(0x5020630)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wh...', 283, 'execute') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1782
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::call('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'receive', 'HASH(0x55342f8)', 2) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Wheel/ReadWrite.pm line 283
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Wheel::ReadWrite::__ANON__(undef, 'POE::Session=ARRAY(0x552f4c8)', 'POE::Kernel=ARRAY(0x4fb7698)', 'HASH(0x552f438)', 'POE::Wheel::ReadWrite(2) -> select read', 'POE::Session=ARRAY(0x552f4c8)', undef, '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Session.pm line 463
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Session::_invoke_state('POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1100
2014-12-12 20:28:21-0600 blade14-2-13: 	eval {...} called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1099
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_dispatch_event('POE::Kernel=ARRAY(0x4fb7698)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Session=ARRAY(0x552f4c8)', 'POE::Wheel::ReadWrite(2) -> select read', 1024, 'ARRAY(0x552e700)', '/gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Re...', 223, undef, ...) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Resource/FileHandles.pm line 221
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::_data_handle_enqueue_ready(undef, undef, 5) called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 285
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_do_timeslice('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Loop/Select.pm line 323
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::loop_run('POE::Kernel=ARRAY(0x4fb7698)') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/POE/Kernel.pm line 1276
2014-12-12 20:28:21-0600 blade14-2-13: 	POE::Kernel::run('POE::Kernel') called at /gsc/scripts/opt/genome/snapshots/genome-3551/lib/perl/Workflow/Server/Worker.pm line 52
2014-12-12 20:28:21-0600 blade14-2-13: 	Workflow::Server::Worker::start('Workflow::Server::Worker', 'blade9-3-10.gsc.wustl.edu', 51228) called at -e line 1
WORKFLOW_ERROR
