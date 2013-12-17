-- Deploy subject.misc_note.body_text
-- requires: subject_misc_note

BEGIN;

CREATE INDEX misc_note_body_index on subject.misc_note using btree (body_text);

COMMIT;
