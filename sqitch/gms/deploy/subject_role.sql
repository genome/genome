-- Deploy subject_role
-- requires: subject_schema

BEGIN;

CREATE TABLE IF NOT EXISTS subject.role (
    id character varying(32) NOT NULL,
    name character varying(64) NOT NULL,
    CONSTRAINT role_pkey PRIMARY KEY (id),
    CONSTRAINT role_name_key UNIQUE (name)
);

COMMIT;
