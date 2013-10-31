-- Revert result_software_result

BEGIN;

DROP TABLE IF EXISTS result.software_result;

COMMIT;
