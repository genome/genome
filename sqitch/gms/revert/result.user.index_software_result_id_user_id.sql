-- Revert result.user.index_software_result_id_user_id

BEGIN;

DROP INDEX result.user_result_id_user_id_index;

COMMIT;
