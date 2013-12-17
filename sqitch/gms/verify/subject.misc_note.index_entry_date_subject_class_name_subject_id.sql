-- Verify subject.misc_note.index_entry_date_subject_class_name_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_note_subject_date_index';

ROLLBACK;
