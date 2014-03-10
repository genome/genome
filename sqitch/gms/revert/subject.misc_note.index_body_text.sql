-- Revert subject.misc_note.index_body_text

BEGIN;

DROP INDEX subject.misc_note_body_index;

COMMIT;
