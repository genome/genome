-- Deploy subject_subject_attribute
-- requires: subject_subject

BEGIN;

CREATE TABLE IF NOT EXISTS subject.subject_attribute (
    subject_id character varying(32) NOT NULL,
    attribute_label character varying(64) NOT NULL,
    attribute_value character varying(512) NOT NULL,
    nomenclature character varying(64) NOT NULL,
    CONSTRAINT subject_attribute_pkey PRIMARY KEY (subject_id, attribute_label, attribute_value, nomenclature),
    CONSTRAINT subject_attribute_subject_id_fkey FOREIGN KEY (subject_id) REFERENCES subject.subject(subject_id)
);

COMMIT;
