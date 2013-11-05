-- Deploy subject.misc_note.subject_class_name_subject_id_header_text
-- requires: subject_misc_note

BEGIN;

CREATE INDEX misc_note_subject_index on subject.misc_note using btree (subject_class_name, subject_id, header_text);

COMMIT;
