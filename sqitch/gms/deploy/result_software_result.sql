-- Deploy result_software_result
-- requires: result_schema

BEGIN;

CREATE TABLE IF NOT EXISTS result.software_result (
    id character varying(32) NOT NULL,
    class_name character varying(255) NOT NULL,
    version character varying(64),
    inputs_id character varying(4000),
    params_id character varying(4000),
    outputs_path character varying(1000),
    lookup_hash character varying(32) NOT NULL,
    CONSTRAINT software_result_pkey PRIMARY KEY (id),
    CONSTRAINT sr_cni UNIQUE (id, class_name)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
