-- Revert result.user.index_software_result_id

BEGIN;

DROP INDEX result.user_software_result_id_idx;

COMMIT;
