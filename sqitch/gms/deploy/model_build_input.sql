-- Deploy model_build_input
-- requires: model_build

BEGIN;

CREATE TABLE IF NOT EXISTS model.build_input (
    build_id character varying(64) NOT NULL,
    value_class_name character varying(255) NOT NULL,
    value_id character varying(1000) NOT NULL,
    name character varying(255) NOT NULL,
    filter_desc character varying(100),
    CONSTRAINT build_input_pkey PRIMARY KEY (build_id, value_class_name, value_id, name),
    CONSTRAINT build_input_build_id_fkey FOREIGN KEY (build_id) REFERENCES model.build(build_id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
