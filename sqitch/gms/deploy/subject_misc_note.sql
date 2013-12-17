-- Deploy subject_misc_note
-- requires: subject_subject

BEGIN;

CREATE TABLE IF NOT EXISTS subject.misc_note (
    editor_id character varying(200) NOT NULL,
    entry_date timestamp(6) without time zone NOT NULL,
    subject_class_name character varying(255) NOT NULL,
    subject_id character varying(255) NOT NULL,
    header_text character varying(255) NOT NULL,
    body_text character varying(4000),
    id character varying(32) NOT NULL,
    CONSTRAINT misc_note_pkey PRIMARY KEY (id)
)
WITH (autovacuum_enabled=false, toast.autovacuum_enabled=false);

COMMIT;
