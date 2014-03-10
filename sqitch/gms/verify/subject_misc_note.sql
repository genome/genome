-- Verify subject_misc_note

BEGIN;

SELECT editor_id, entry_date, subject_class_name, subject_id, header_text, body_text, id
FROM subject.misc_note
WHERE FALSE;

ROLLBACK;
