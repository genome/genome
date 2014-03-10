-- Revert result.software_result.index_class_name

BEGIN;

DROP INDEX result.sr_cname;

COMMIT;
