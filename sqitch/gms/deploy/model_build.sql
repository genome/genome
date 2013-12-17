-- Deploy model_build
-- requires: model_model

BEGIN;

CREATE TABLE IF NOT EXISTS model.build (
    build_id character varying(64) NOT NULL,
    data_directory character varying(1000),
    model_id character varying(64) NOT NULL,
    software_revision character varying(1000),
    subclass_name character varying(255),
    CONSTRAINT build_pkey PRIMARY KEY (build_id),
    CONSTRAINT build_model_id_fkey FOREIGN KEY (model_id) REFERENCES model.model(genome_model_id)
);

COMMIT;
