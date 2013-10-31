-- Deploy result_metric
-- requires: result_software_result

BEGIN;

CREATE TABLE IF NOT EXISTS result.metric (
    software_result_id character varying(32) NOT NULL,
    metric_name character varying(1000) NOT NULL,
    metric_value character varying(1000) NOT NULL,
    CONSTRAINT metric_pkey PRIMARY KEY (software_result_id, metric_name),
    CONSTRAINT metric_software_result_id_fkey FOREIGN KEY (software_result_id) REFERENCES result.software_result(id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
