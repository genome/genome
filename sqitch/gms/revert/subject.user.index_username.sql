-- Revert subject.user.index_username

BEGIN;

DROP INDEX subject.subject_user_username_index;

COMMIT;
