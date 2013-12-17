-- Revert subject.misc_note.index_subject_class_name_subject_id_header_text

BEGIN;

DROP INDEX subject.misc_note_subject_index;

COMMIT;
