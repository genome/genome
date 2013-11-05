-- Deploy subject.misc_note.subject_id
-- requires: subject_misc_note

BEGIN;

CREATE INDEX misc_note_subject_id on subject.misc_note using btree (subject_id);

COMMIT;
