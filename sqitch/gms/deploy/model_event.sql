-- Deploy model_event
-- requires: model_build
-- requires: model_model

BEGIN;

CREATE TABLE IF NOT EXISTS model.event (
    genome_model_event_id character varying(32) NOT NULL,
    model_id character varying(32),
    run_id character varying(32),
    event_type character varying(255) NOT NULL,
    event_status character varying(255),
    lsf_job_id character varying(64),
    date_scheduled timestamp(6) without time zone,
    date_completed timestamp(6) without time zone,
    user_name character varying(64),
    ref_seq_id character varying(64),
    retry_count integer,
    parent_event_id character varying(32),
    prior_event_id character varying(32),
    status_detail character varying(200),
    build_id character varying(32),
    instrument_data_id character varying(64),
    workflow_instance_id character varying(32),
    CONSTRAINT event_pkey PRIMARY KEY (genome_model_event_id),
    CONSTRAINT fk_m_e_bi FOREIGN KEY (build_id) REFERENCES model.build(build_id),
    CONSTRAINT fk_m_e_mi FOREIGN KEY (model_id) REFERENCES model.model(genome_model_id),
    CONSTRAINT fk_m_e_pei FOREIGN KEY (parent_event_id) REFERENCES model.event(genome_model_event_id),
    CONSTRAINT fk_m_e_prei FOREIGN KEY (prior_event_id) REFERENCES model.event(genome_model_event_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
