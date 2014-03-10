-- Deploy result_param
-- requires: result_software_result

BEGIN;

CREATE TABLE IF NOT EXISTS result.param (
    software_result_id character varying(32) NOT NULL,
    param_name character varying(100) NOT NULL,
    param_value character varying(1000) NOT NULL,
    name character varying(255),
    value_class_name character varying(255),
    value_id character varying(1000),
    CONSTRAINT param_pkey PRIMARY KEY (software_result_id, param_name),
    CONSTRAINT param_software_result_id_fkey FOREIGN KEY (software_result_id) REFERENCES result.software_result(id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
