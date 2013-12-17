-- Deploy subject.subject.name
-- requires: subject_subject

BEGIN;

CREATE INDEX subject_name_index on subject.subject using btree (name);

COMMIT;
