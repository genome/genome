-- Deploy model_event_metric
-- requires: model_event

BEGIN;

CREATE TABLE model.event_metric (
    event_id character varying(32) NOT NULL,
    metric_name character varying(100) NOT NULL,
    metric_value character varying(1000),
    CONSTRAINT event_metric_pkey PRIMARY KEY (event_id, metric_name),
    CONSTRAINT event_metric_event_id_fkey FOREIGN KEY (event_id) REFERENCES model.event(genome_model_event_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
