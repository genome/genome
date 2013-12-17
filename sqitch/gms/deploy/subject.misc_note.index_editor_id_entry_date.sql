-- Deploy subject.misc_note.editor_id_entry_date
-- requires: subject_misc_note

BEGIN;

CREATE INDEX misc_note_editor_index on subject.misc_note using btree (editor_id, entry_date);

COMMIT;
