-- Deploy model_processing_profile_param
-- requires: model_processing_profile

BEGIN;

CREATE TABLE IF NOT EXISTS model.processing_profile_param (
    processing_profile_id character varying(32) NOT NULL,
    param_name character varying(100) NOT NULL,
    param_value character varying(1000) NOT NULL,
    name character varying(255),
    value_class_name character varying(255),
    value_id character varying(1000),
    CONSTRAINT processing_profile_param_pkey PRIMARY KEY (processing_profile_id, param_name, param_value),
    CONSTRAINT processing_profile_param_processing_profile_id_fkey FOREIGN KEY (processing_profile_id) REFERENCES model.processing_profile(id)
);

COMMIT;
