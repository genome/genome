-- Deploy model_event_input
-- requires: model_build
-- requires: model_event

BEGIN;

CREATE TABLE IF NOT EXISTS model.event_input (
    event_id character varying(32) NOT NULL,
    param_name character varying(100) NOT NULL,
    param_value character varying(1000) NOT NULL,
    CONSTRAINT event_input_pkey PRIMARY KEY (event_id, param_name, param_value),
    CONSTRAINT event_input_event_id_fkey FOREIGN KEY (event_id) REFERENCES model.event(genome_model_event_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
