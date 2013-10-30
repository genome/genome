-- Deploy model_model_input
-- requires: model_model

BEGIN;

CREATE TABLE IF NOT EXISTS model.model_input (
    model_id character varying(32) NOT NULL,
    value_class_name character varying(255) NOT NULL,
    value_id character varying(1000) NOT NULL,
    name character varying(255) NOT NULL,
    filter_desc character varying(100),
    CONSTRAINT model_input_pkey PRIMARY KEY (model_id, value_class_name, value_id, name),
    CONSTRAINT model_input_model_id_fkey FOREIGN KEY (model_id) REFERENCES model.model(genome_model_id)
);

COMMIT;
