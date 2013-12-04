-- Revert disk_allocation

BEGIN;

DROP TABLE IF EXISTS disk.allocation;

COMMIT;
