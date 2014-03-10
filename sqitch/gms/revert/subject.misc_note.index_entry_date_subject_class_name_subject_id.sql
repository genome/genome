-- Revert subject.misc_note.index_entry_date_subject_class_name_subject_id

BEGIN;

DROP INDEX subject.misc_note_subject_date_index;

COMMIT;
