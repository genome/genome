-- Deploy model_processing_profile
-- requires: model_schema

BEGIN;

CREATE TABLE IF NOT EXISTS model.processing_profile (
    id character varying(32) NOT NULL,
    type_name character varying(255),
    name character varying(255),
    subclass_name character varying(255),
    CONSTRAINT processing_profile_pkey PRIMARY KEY (id)
);

COMMIT;
