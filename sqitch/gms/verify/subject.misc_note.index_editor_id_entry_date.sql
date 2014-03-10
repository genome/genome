-- Verify subject.misc_note.index_editor_id_entry_date

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'misc_note_editor_index';

ROLLBACK;
