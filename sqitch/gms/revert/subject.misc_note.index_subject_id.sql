-- Revert subject.misc_note.index_subject_id

BEGIN;

DROP INDEX subject.misc_note_subject_id;

COMMIT;
