-- Deploy subject_project
-- requires: subject_schema

BEGIN;

CREATE TABLE IF NOT EXISTS subject.project (
    id character varying(64) NOT NULL,
    name character varying(200) NOT NULL,
    CONSTRAINT project_pkey PRIMARY KEY (id)
);

COMMIT;
