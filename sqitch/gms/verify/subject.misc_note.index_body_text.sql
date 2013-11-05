-- Verify subject.misc_note.index_body_text

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_note_body_index';

ROLLBACK;
