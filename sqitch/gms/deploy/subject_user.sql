-- Deploy subject_user
-- requires: subject_schema

BEGIN;

CREATE TABLE IF NOT EXISTS subject."user" (
    name character varying(64),
    email character varying(255) NOT NULL,
    username character varying(64),
    CONSTRAINT user_pkey PRIMARY KEY (email)
);

COMMIT;
