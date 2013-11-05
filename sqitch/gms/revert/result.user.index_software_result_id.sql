-- Revert result.user.index_software_result_id

BEGIN;

DROP INDEX result.sru_rid_i;

COMMIT;
