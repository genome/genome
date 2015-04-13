package Genome::Model::Tools::CgHub::Test;

use strict;
use warnings;

use LWP::UserAgent;
use Sub::Install;
use Test::MockObject;

sub overload_lwp_user_agent_request {
    my $class = shift;

    my $request_cnt = 0;
    my $response = Test::MockObject->new();
    $response->mock('content', sub{ return xml_content(); });
    Sub::Install::reinstall_sub({
            code => sub{
                $request_cnt++;
                if ( ${$_[1]->uri} =~ /INVALID/ ) {
                    $response->set_false('is_success');
                }
                else {
                    $response->set_true('is_success');
                }
                return $response;
            },
            into => 'LWP::UserAgent',
            as => 'request',
        });

    return 1;
}

sub xml_content {
    return <<XML;
<?xml version="1.0" encoding="utf-8" standalone="yes"?>
<ResultSet date="2015-04-09 23:51:37" id="32418000">
	<Query>analysis_id:01f22763-6bb2-4edc-b65a-99904f7c6fad</Query>
	<Hits>1</Hits>
	<Result id="1">
		<analysis_id>01f22763-6bb2-4edc-b65a-99904f7c6fad</analysis_id>
		<state>live</state>
		<reason></reason>
		<last_modified>2014-10-24T07:36:34Z</last_modified>
		<upload_date>2014-10-24T17:06:21Z</upload_date>
		<published_date>2014-10-24T07:36:34Z</published_date>
		<center_name>BCCAGSC</center_name>
		<study>phs000178</study>
		<aliquot_id>f43178e9-176f-4c0b-b31e-44a6f485e02e</aliquot_id>
		<files>
			<file>
				<filename>TCGA-L5-A8NH-01A-11R-A37I-31_rnaseq.bam</filename>
				<filesize>14658974566</filesize>
				<checksum type="MD5">530aedb9ae22815fec6e21c68edcdee6</checksum>
			</file>
			<file>
				<filename>TCGA-L5-A8NH-01A-11R-A37I-31_rnaseq.bam.bai</filename>
				<filesize>6368144</filesize>
				<checksum type="MD5">5f476af9bf5d15a508f1e592bc40adc6</checksum>
			</file>
		</files>
		<sample_accession></sample_accession>
		<legacy_sample_id>TCGA-L5-A8NH-01A-11R-A37I-31</legacy_sample_id>
		<disease_abbr>ESCA</disease_abbr>
		<tss_id>L5</tss_id>
		<participant_id>5af0e222-3dc8-400b-ba61-2225921f2fd3</participant_id>
		<sample_id>7d82d2b8-85e1-4959-b440-d4c4cf2fbbb4</sample_id>
		<analyte_code>R</analyte_code>
		<sample_type>01</sample_type>
		<library_strategy>RNA-Seq</library_strategy>
		<platform>ILLUMINA</platform>
		<refassem_short_name>GRCh37-lite</refassem_short_name>
		<analysis_submission_uri>https://cghub.ucsc.edu/cghub/metadata/analysisSubmission/01f22763-6bb2-4edc-b65a-99904f7c6fad</analysis_submission_uri>
		<analysis_full_uri>https://cghub.ucsc.edu/cghub/metadata/analysisFull/01f22763-6bb2-4edc-b65a-99904f7c6fad</analysis_full_uri>
		<analysis_data_uri>https://cghub.ucsc.edu/cghub/data/analysis/download/01f22763-6bb2-4edc-b65a-99904f7c6fad</analysis_data_uri>
	</Result>
	<ResultSummary>
		<downloadable_file_count>2</downloadable_file_count>
		<downloadable_file_size units="GB">13.66</downloadable_file_size>
		<state_count>
			<live>1</live>
		</state_count>
	</ResultSummary>
</ResultSet>
XML
;
}

1;

