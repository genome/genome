-- Deploy subject.misc_note.entry_date_subject_class_name_subject_id
-- requires: subject_misc_note

BEGIN;

CREATE INDEX misc_note_subject_date_index on subject.misc_note using btree (entry_date, subject_class_name, subject_id);

COMMIT;
