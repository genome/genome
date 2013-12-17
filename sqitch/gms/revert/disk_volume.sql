-- Revert disk_volume

BEGIN;

DROP TABLE IF EXISTS disk.volume;

COMMIT;
