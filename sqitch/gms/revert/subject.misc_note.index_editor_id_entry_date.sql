-- Revert subject.misc_note.index_editor_id_entry_date

BEGIN;

DROP INDEX subject.misc_note_editor_index;

COMMIT;
