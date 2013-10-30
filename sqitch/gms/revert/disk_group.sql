-- Revert disk_group

BEGIN;

DROP TABLE IF EXISTS disk."group";

COMMIT;
