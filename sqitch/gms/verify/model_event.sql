-- Verify model_event

BEGIN;

SELECT genome_model_event_id, model_id, run_id, event_type, event_status,
    lsf_job_id, date_scheduled, date_completed, user_name, ref_seq_id,
    retry_count, parent_event_id,prior_event_id, status_detail, build_id,
    instrument_data_id, workflow_instance_id
FROM model.event
WHERE FALSE;

ROLLBACK;
