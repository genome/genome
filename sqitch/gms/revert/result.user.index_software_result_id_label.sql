-- Revert result.user.index_software_result_id_label

BEGIN;

DROP INDEX result.user_result_label_index;

COMMIT;
