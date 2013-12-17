-- Deploy model_build_metric
-- requires: model_build

BEGIN;

CREATE TABLE IF NOT EXISTS model.build_metric (
    build_id character varying(32) NOT NULL,
    metric_name character varying(100) NOT NULL,
    metric_value character varying(1000) NOT NULL,
    CONSTRAINT build_metric_pkey PRIMARY KEY (build_id, metric_name, metric_value),
    CONSTRAINT build_metric_build_id_fkey FOREIGN KEY (build_id) REFERENCES model.build(build_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
