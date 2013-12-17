-- Deploy subject_subject
-- requires: subject_schema

BEGIN;

CREATE TABLE IF NOT EXISTS subject.subject (
    subject_id character varying(32) NOT NULL,
    subclass_name character varying(255) NOT NULL,
    name character varying(255),
    CONSTRAINT subject_pkey PRIMARY KEY (subject_id)
);

COMMIT;
