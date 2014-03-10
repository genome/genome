-- Deploy model_model_group
-- requires: model_schema

BEGIN;

CREATE TABLE IF NOT EXISTS model.model_group (
    id character varying(64) NOT NULL,
    name character varying(255) NOT NULL,
    user_name character varying(64),
    uuid character varying(64),
    CONSTRAINT model_group_pkey PRIMARY KEY (id),
    CONSTRAINT uniq_m_mg_name UNIQUE (name)
);

COMMIT;
