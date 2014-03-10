-- Verify subject.misc_note.index_subject_class_name_subject_id_header_text

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_note_subject_index';

ROLLBACK;
