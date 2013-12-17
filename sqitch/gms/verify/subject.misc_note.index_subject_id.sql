-- Verify subject.misc_note.index_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_note_subject_id';

ROLLBACK;
