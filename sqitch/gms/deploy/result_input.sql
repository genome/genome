-- Deploy result_input
-- requires: result_software_result

BEGIN;

CREATE TABLE IF NOT EXISTS result.input (
    software_result_id character varying(32) NOT NULL,
    input_name character varying(100) NOT NULL,
    input_value character varying(1000) NOT NULL,
    name character varying(255),
    value_class_name character varying(255),
    value_id character varying(1000),
    CONSTRAINT input_pkey PRIMARY KEY (software_result_id, input_name),
    CONSTRAINT input_software_result_id_fkey FOREIGN KEY (software_result_id) REFERENCES result.software_result(id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
