-- Revert result.user.index_user_id

BEGIN;

DROP INDEX result.sru_uid_i;

COMMIT;
